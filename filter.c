/********************************************************************************************
 *
 *  Fast local alignment filter for long, noisy reads based on "dumbing down" of my RECOMB 2005
 *     filter with Jens Stoye, and a "smarting up" of the k-mer matching by turning it into
 *     a threaded sort and merge paradigm using a super cache coherent radix sort.  Local
 *     alignment is accomplised with dynamically-banded O(nd) algorithm that terminates when
 *     it fails to find a e-matching patch for a significant distance, and polishes the match
 *     to the last e-prefix-positive 32-mer.
 *
 *  Author :  Gene Myers
 *  First  :  June 2013
 *  Current:  June 1, 2014
 *
 ********************************************************************************************/

//  A complete threaded code for the filter

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

#include "DB.h"
#include "lsd.sort.h"
#include "filter.h"
#include "align.h"

#undef  PROFILE        //  WHen running sensitivity trials, compute histogram of
#define MAXHIT   1000  //    false & true positive hit scores

  //  K-mer selection strategy control:
  //       MODIFER off: Select minimizers in window of size WINDOW
  //       MODIFER on : Select modimizers mod MODULUS < MODTHR (best)

#define MODIFER
#define WINDOW   6

#define MODULUS  101
#define MODTHR    28

  //  Debug Controls

#undef    TEST_KSORT
#undef    TEST_PAIRS
#undef    TEST_CSORT
#define    HOW_MANY   3000   //  Print first HOW_MANY items for each of the TEST options above

#define DO_ALIGNMENT
#define DO_BRIDGING

#undef  TEST_GATHER
#undef  TEST_CONTAIN
#undef  TEST_BRIDGE
#undef  SHOW_OVERLAP          //  Show the cartoon
#undef  SHOW_ALIGNMENT        //  Show the alignment

#define   ALIGN_WIDTH    80   //     Parameters for alignment
#define   ALIGN_INDENT   20
#define   ALIGN_BORDER   10

#ifdef SHOW_OVERLAP
#define NOTHREAD
#endif

#ifdef TEST_GATHER
#define NOTHREAD
#endif

#ifdef TEST_CONTAIN
#define NOTHREAD
#endif

  //  Algorithm constants & global data types

#define THREAD       pthread_t

#define MAX_CODE_16  0xffffu
#define MAX_CODE_32  0xffffffffu
#define MAX_CODE_64  0xffffffffffffffffllu

#define SIGN_BIT     0x1u
#define LONG_BIT     0x80000000u
#define POST_MASK    0x7fffffffu

#define MAXGRAM 10000  //  Cap on k-mer count histogram (in count_thread, merge_thread)

#define PANEL_SIZE     50000   //  Size to break up very long A-reads
#define PANEL_OVERLAP  10000   //  Overlap of A-panels

#define MATCH_CHUNK    100     //  Max initial number of hits between two reads
#define TRACE_CHUNK  20000     //  Max initial trace points in hits between two reads

typedef struct
  { uint32 rpos;
    uint32 read;
    uint64 code;
  } KmerPos;

typedef struct
  { int aread;
    int bread;
    int apos;
    int diag;
  } SeedPair;

/*******************************************************************************************
 *
 *  PARAMETER SETUP
 *
 ********************************************************************************************/

static int    Kmer;
static int    Koff;           //  Kmer + 1;
static int    Kshift;         //  2*Kmer
static uint64 Kmask;          //  2^Kshift - 1

static uint64 LFmask;         //  4^floor(Kmer/2)-1
static uint64 HFmask;         //  Kmask - LFmask;
static uint64 LRmask;         //  4^ceil(Kmer/2)-1
static uint64 HRmask;         //  Kmask - LRmask;

static int Hitmin;
static int Binshift;
static int Suppress;
static int TooFrequent;       //  (Suppress != 0) ? Suppress : INT32_MAX

static int NTHREADS;          //  # of threads to use

void Set_Filter_Params(int kmer, int binshift, int suppress, int hitmin, int nthread)
{ if (kmer > 32)
    { fprintf(stderr,"%s: Kmer length must be <= 32\n",Prog_Name);
      exit (1);
    }

  Kmer     = kmer;
  Koff     = kmer+1;
  Binshift = binshift;
  Suppress = suppress;
  Hitmin   = hitmin;

  Kmer    = kmer;
  if (Kmer >= 32)
    { Kshift = 64;
      Kmask  = MAX_CODE_64;
    }
  else
    { Kshift = 2*Kmer;
      Kmask  = (0x1llu << Kshift) - 1;
    }

  LFmask = (0x1llu<<(Kshift/2))-1;
  HFmask = Kmask - LFmask;

  LRmask = (0x1llu<<((Kshift+1)/2))-1;
  HRmask = Kmask - LRmask;

  if (Suppress == 0)
    TooFrequent = INT32_MAX;
  else
    TooFrequent = Suppress;

  NTHREADS = nthread;
}


/*******************************************************************************************
 *
 *  INDEX BUILD
 *
 ********************************************************************************************/

static DAZZ_DB    *TA_block;
static DAZZ_TRACK *TA_track;

static KmerPos *FR_src;
static KmerPos *FR_trg;

static uint64 Cumber[4];   //  Cumber[i] = (3-i) << (Kshift-2)

typedef struct
  { int    beg;
    int    end;
    int    fill;
  } Tuple_Arg;

  //  for reads [beg,end) computing how many k-tuples are not masked

static void *mask_thread(void *arg)
{ Tuple_Arg *data  = (Tuple_Arg *) arg;
  DAZZ_READ *reads = TA_block->reads;
  int        km1   = Kmer-1;
  int        beg, end, idx;
  int64      a, b, f;
  int        i, p, q;

#ifndef MODIFER
  uint64     min1[WINDOW], min2[WINDOW];
  int        w, ny, m1, m2;
#endif
  int        x;
  uint64     c, u;
  uint64     d, v;
  char      *s;

  beg = data->beg;
  end = data->end;
  idx = 0;

  s = ((char *) (TA_block->bases)) + TA_block->reads[beg].boff;
  if (TA_track != NULL)
    { int64 *anno1 = ((int64 *) (TA_track->anno)) + 1;
      int   *point = (int *) (TA_track->data);

      q = 0;
      f = anno1[beg-1];
      for (i = beg; i < end; i++)
        { b = f;
          f = anno1[i];
          for (a = b; a <= f; a += 2)
            { if (a == b)
                p = 0;
              else
                p = point[a-1];
              if (a == f)
                q = reads[i].rlen;
              else
                q = point[a];
              if (q-p > km1)
                { c = u = 0;
                  for (x = 1; x < Kmer; x++)
                    { c = (c << 2) | s[p];
                      u = (u >> 2) | Cumber[(int) s[p++]];
                    }
#ifndef MODIFER
                  ny = 1;
                  w  = 0;
#endif
                  while (p < q)
                    { x = s[p];

                      d = (c & HFmask);
                      c  = ((c << 2) | x) & Kmask;
                      d = d | (c & LFmask);

                      v = (u & LRmask);
                      u  = (u >> 2) | Cumber[x];
                      v = v | (u & HRmask);

#ifdef MODIFER
                      if (u < c)
                        { if (u % MODULUS < MODTHR)
                            idx += 1;
                        }
                      else
                        { if (c % MODULUS < MODTHR)
                            idx += 1;
                        }

                      if (v < d)
                        { if (v % MODULUS < MODTHR)
                            idx += 1;
                        }
                      else
                        { if (d % MODULUS < MODTHR)
                            idx += 1;
                        }
#else
                      if (u < c)
                        min1[w] = u;
                      else
                        min1[w] = c;
                      if (v < d)
                        min2[w] = v;
                      else
                        min2[w] = d;

                      if (ny)
                        { w += 1;
                          if (w < WINDOW)
                            continue;
                          m1 = m2 = w = w-1;
                          ny = 0;
                        }

                      if (w == m1)
                        { m1 = 0;
                          for (x = 1; x < WINDOW; x++)
                            if (min1[x] < min1[m1])
                              m1 = x;
                          idx += 1;
                        }
                      else if (min1[w] < min1[m1])
                        { m1 = w;
                          idx += 1;
                        } 

                      if (w == m2)
                        { m2 = 0;
                          for (x = 1; x < WINDOW; x++)
                            if (min2[x] < min2[m2])
                              m2 = x;
                          idx += 1;
                        } 
                      else if (min2[w] < min2[m2])
                        { m2 = w;
                          idx += 1;
                        } 
                        
                      w += 1;
                      if (w == WINDOW)
                        w = 0; 
#endif
                      p += 1;
                    }   
                }   
            }
          s += (q+1);
        }
    }
  else
    for (i = beg; i < end; i++)
      { q = reads[i].rlen;
        c = u = 0;
        for (p = 0; p < km1; p++)
          { x = s[p];
            c = (c << 2) | x;
            u = (u >> 2) | Cumber[x];
          }
#ifndef MODIFER
        ny = 1;
        w  = 0;
#endif
        for (p = km1; p < q; p++)
          { x = s[p];

            d = (c & HFmask);
            c  = ((c << 2) | x) & Kmask;
            d = d | (c & LFmask);

            v = (u & LRmask);
            u  = (u >> 2) | Cumber[x];
            v = v | (u & HRmask);

#ifdef MODIFER
            if (u < c)
              { if (u % MODULUS < MODTHR)
                  idx += 1;
              }
            else
              { if (c % MODULUS < MODTHR)
                  idx += 1;
              }

            if (v < d)
              { if (v % MODULUS < MODTHR)
                  idx += 1;
              }
            else
              { if (d % MODULUS < MODTHR)
                  idx += 1;
              }
#else
            if (u < c)
              min1[w] = u;
            else
              min1[w] = c;
            if (v < d)
              min2[w] = v;
            else
              min2[w] = d;

            if (ny)
              { w += 1;
                if (w < WINDOW)
                  continue;
                m1 = m2 = w = w-1;
                ny = 0;
              }

            if (w == m1)
              { m1 = 0;
                for (x = 1; x < WINDOW; x++)
                  if (min1[x] < min1[m1])
                    m1 = x;
                idx += 1;
              }
            else if (min1[w] < min1[m1])
              { m1 = w;
                idx += 1;
              }

            if (w == m2)
              { m2 = 0;
                for (x = 1; x < WINDOW; x++)
                  if (min2[x] < min2[m2])
                    m2 = x;
                idx += 1;
              }
            else if (min2[w] < min2[m2])
              { m2 = w;
                idx += 1;
              }
              
            w += 1;
            if (w == WINDOW)
              w = 0;
#endif
          }    
        s += (q+1);
      }

  data->fill = idx;
  return (NULL);
}

  // for reads [beg,end) place their k-tuples in list starting at index idx

static void *tuple_thread(void *arg)
{ Tuple_Arg *data  = (Tuple_Arg *) arg;
  DAZZ_READ *reads = TA_block->reads;
  int        km1   = Kmer-1;
  KmerPos   *list  = FR_src;
  int        beg, end, idx;
  int64     a, b, f;
  int       i, p, q, x, r;
  uint64    c, u;
  uint64    d, v;
  uint32    lbit;
  char     *s;

#ifndef MODIFER
  uint64     min1[WINDOW], min2[WINDOW];
  int        sgn1[WINDOW], sgn2[WINDOW];
  int        pos[WINDOW];
  int        w, ny, m1, m2;
#endif

  beg = data->beg;
  end = data->end;
  idx = data->fill;

  s = ((char *) (TA_block->bases)) + TA_block->reads[beg].boff;
  if (TA_track != NULL)
    { int64     *anno1 = ((int64 *) (TA_track->anno)) + 1;
      int       *point = (int *) (TA_track->data);

      q = 0;
      f = anno1[beg-1];
      for (i = beg; i < end; i++)
        { r = (i << 1);
          b = f;
          f = anno1[i];
          for (a = b; a <= f; a += 2)
            { if (a == b)
                p = 0;
              else
                p = point[a-1];
              if (a == f)
                q = reads[i].rlen;
              else
                q = point[a];
              if (q-p > km1)
                { c = 0;
                  u = 0;
                  for (x = 1; x < Kmer; x++)
                    { c = ((c << 2) | s[p]);
                      u = (u >> 2) | Cumber[(int) s[p]];
                      p += 1;
                    }
#ifndef MODIFER
                  ny = 1;
                  w  = 0;
#endif

                  lbit = 0;
                  while (p < q)
                    { x = s[p++];

                      d = (c & HFmask);
                      c = ((c << 2) | x) & Kmask;
                      d = d | (c & LFmask);

                      v = (u & LRmask);
                      u = (u >> 2) | Cumber[x];
                      v = v | (u & HRmask);
#ifdef MODIFER
                      if (u < c)
                        { if (u % MODULUS < MODTHR)
                            { list[idx].code = u;
                              list[idx].read = r | SIGN_BIT;
                              list[idx].rpos = p;
		              idx += 1;
                            }
                        }
                      else
                        { if (c % MODULUS < MODTHR)
                            { list[idx].code = c;
                              list[idx].read = r;
                              list[idx].rpos = p;
		              idx += 1;
                            }
                        }

                      if (v < d)
                        { if (v % MODULUS < MODTHR)
                            { list[idx].code = v;
                              list[idx].read = r | SIGN_BIT;
                              list[idx].rpos = p | lbit;
                              idx += 1;
                            }
                        }
                      else
                        { if (d % MODULUS < MODTHR)
                            { list[idx].code = d;
                              list[idx].read = r;
                              list[idx].rpos = p | lbit;
                              idx += 1;
                            }
                        }
#else
                      if (u < c)
                        { min1[w] = u;
                          sgn1[w] = 0x1;
                        }
                      else
                        { min1[w] = c;
                          sgn1[w] = 0x0;
                        }
                      if (v < d)
                        { min2[w] = v;
                          sgn2[w] = 0x1;
                        }
                      else
                        { min2[w] = d;
                          sgn2[w] = 0x0;
                        }
                      pos[w] = p;

                      if (ny)
                        { w += 1;
                          if (w < WINDOW)
                            continue;
                          m1 = m2 = w = w-1;
                          ny = 0;
                        }

                      if (w == m1)
                        { m1 = 0;
                          for (x = 1; x < WINDOW; x++)
                            if (min1[x] < min1[m1])
                              m1 = x;
                          list[idx].read = r | sgn1[m1];
                          list[idx].code = min1[m1];
                          list[idx].rpos = pos[m1];
                          idx += 1;
                        }
                      else if (min1[w] < min1[m1])
                        { m1 = w;
                          list[idx].read = r | sgn1[m1];
                          list[idx].code = min1[m1];
                          list[idx].rpos = pos[m1];
                          idx += 1;
                        }

                      if (w == m2)
                        { m2 = 0;
                          for (x = 1; x < WINDOW; x++)
                            if (min2[x] < min2[m2])
                              m2 = x;
                          list[idx].read = r | sgn2[m2];
                          list[idx].code = min2[m2];
                          list[idx].rpos = pos[m2] | lbit;
                          idx += 1;
                        }
                      else if (min2[w] < min2[m2])
                        { m2 = w;
                          list[idx].read = r | sgn2[m2];
                          list[idx].code = min2[m2];
                          list[idx].rpos = pos[m2] | lbit;
                          idx += 1;
                        }

                      w += 1;
                      if (w == WINDOW)
                        w = 0;
#endif
                      lbit = LONG_BIT;
                    }
                }
            }
          s += (q+1);
        }
    }

  else
    for (i = beg; i < end; i++)
      { q = reads[i].rlen;
        r = (i << 1);
        c = 0;
        u = 0;
        for (p = 0; p < km1; p++)
          { x = s[p];
            c = (c << 2) | x;
            u = (u >> 2) | Cumber[x];
          }
#ifndef MODIFER
        ny = 1;
        w  = 0;
#endif
        lbit = 0;
        while (p < q)
          { x = s[p++];

            d = (c & HFmask);
            c = ((c << 2) | x) & Kmask;
            d = d | (c & LFmask);

            v = (u & LRmask);
            u = (u >> 2) | Cumber[x];
            v = v | (u & HRmask);

#ifdef MODIFER
            if (u < c)
              { if (u % MODULUS < MODTHR)
                  { list[idx].code = u;
                    list[idx].read = r | SIGN_BIT;
                    list[idx].rpos = p;
	            idx += 1;
                  }
              }
            else
              { if (c % MODULUS < MODTHR)
                  { list[idx].code = c;
                    list[idx].read = r;
                    list[idx].rpos = p;
	            idx += 1;
                  }
              }

            if (v < d)
              { if (v % MODULUS < MODTHR)
                  { list[idx].code = v;
                    list[idx].read = r | SIGN_BIT;
                    list[idx].rpos = p | lbit;
                    idx += 1;
                  }
              }
            else
              { if (d % MODULUS < MODTHR)
                  { list[idx].code = d;
                    list[idx].read = r;
                    list[idx].rpos = p | lbit;
                    idx += 1;
                  }
              }
#else
            if (u < c)
              { min1[w] = u;
                sgn1[w] = 0x1;
              }
            else
              { min1[w] = c;
                sgn1[w] = 0x0;
              }
            if (v < d)
              { min2[w] = v;
                sgn2[w] = 0x1;
              }
            else
              { min2[w] = d;
                sgn2[w] = 0x0;
              }
            pos[w] = p;

            if (ny)
              { w += 1;
                if (w < WINDOW)
                  continue;
                m1 = m2 = w = w-1;
                ny = 0;
              }

            if (w == m1)
              { m1 = 0;
                for (x = 1; x < WINDOW; x++)
                  if (min1[x] < min1[m1])
                    m1 = x;
                list[idx].read = r | sgn1[m1];
                list[idx].code = min1[m1];
                list[idx].rpos = pos[m1];
                idx += 1;
              }
            else if (min1[w] < min1[m1])
              { m1 = w;
                list[idx].read = r | sgn1[m1];
                list[idx].code = min1[m1];
                list[idx].rpos = pos[m1];
                idx += 1;
              }

            if (w == m2)
              { m2 = 0;
                for (x = 1; x < WINDOW; x++)
                  if (min2[x] < min2[m2])
                    m2 = x;
                list[idx].read = r | sgn2[m2];
                list[idx].code = min2[m2];
                list[idx].rpos = pos[m2] | lbit;
                idx += 1;
              }
            else if (min2[w] < min2[m2])
              { m2 = w;
                list[idx].read = r | sgn2[m2];
                list[idx].code = min2[m2];
                list[idx].rpos = pos[m2] | lbit;
                idx += 1;
              }

            w += 1;
            if (w == WINDOW)
              w = 0;
#endif
            lbit = LONG_BIT;
          }
        s += (q+1);
      }

  return (NULL);
}

static void *compsize_thread(void *arg)
{ Tuple_Arg  *data  = (Tuple_Arg *) arg;
  int         end   = data->end;
  KmerPos    *src   = FR_src;
  int         n, i, c, p;
  uint64      h, g;

  i = data->beg;
  h = src[i].code;
  n = 0;
  while (i < end)
    { p = i++;
      while (1)
        { g = src[i].code; 
          if (g != h)
            break;
          i += 1;
        }
      if ((c = (i-p)) < TooFrequent)
        n += c;
      h = g;
    }

  data->fill = n;
  return (NULL);
}

static void *compress_thread(void *arg)
{ Tuple_Arg  *data  = (Tuple_Arg *) arg;
  int         end   = data->end;
  KmerPos    *src   = FR_src;
  KmerPos    *trg   = FR_trg;
  int         n, i, p;
  uint64      h, g;

  i = data->beg;
  h = src[i].code;
  n = data->fill;
  while (i < end)
    { p = i++;
      while (1)
        { g = src[i].code; 
          if (g != h)
            break;
          i += 1;
        }
      if (i-p < TooFrequent)
        { while (p < i)
            trg[n++] = src[p++];
        }
      h = g;
    }

  return (NULL);
}

void *Sort_Kmers(DAZZ_DB *block, int *len)
{ THREAD    threads[NTHREADS];
  Tuple_Arg parmt[NTHREADS];

  KmerPos  *src, *trg, *rez;
  int       kmers, nreads;

  nreads = block->nreads;

  if (block->reads[nreads].boff > 0x7fffffffll)
    { fprintf(stderr,"%s: Fatal error, DB blocks are greater than 2Gbp!\n",Prog_Name);
      Clean_Exit(1);
    }

  if (nreads <= 0)
    goto no_mers;

  TA_block = block;
  TA_track = block->tracks;

  Cumber[0] = (0x3llu << (Kshift-2));
  Cumber[1] = (0x2llu << (Kshift-2));
  Cumber[2] = (0x1llu << (Kshift-2));
  Cumber[3] = (0x0llu << (Kshift-2));

  //  Determine how many k-tuples will be listed for each thread
  //    and use that to set up index drop points

  { int i, x, z;

    parmt[0].beg = 0;
    for (i = 1; i < NTHREADS; i++)
      parmt[i].beg = parmt[i-1].end = (((int64) nreads) * i) / NTHREADS;
    parmt[NTHREADS-1].end = nreads;

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,mask_thread,parmt+i);
    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

    x = 0;
    for (i = 0; i < NTHREADS; i++)
      { z = parmt[i].fill;
        parmt[i].fill = x;
        x += z;
      }
    kmers = x;

    if (kmers <= 0)
      goto no_mers;
  }

  //  Allocate k-mer sorting arrays now that # of kmers is known

  if (( (Kshift-1)/8 + (TooFrequent < INT32_MAX) ) & 0x1)
    { src = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+2),"Allocating Sort_Kmers vectors");
      trg = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+2),"Allocating Sort_Kmers vectors");
    }
  else
    { trg = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+2),"Allocating Sort_Kmers vectors");
      src = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+2),"Allocating Sort_Kmers vectors");
    }
  if (src == NULL || trg == NULL)
    Clean_Exit(1);

  if (VERBOSE)
    { printf("\n   Kmer count = ");
      Print_Number((int64) kmers,0,stdout);
      printf("\n   Using %.2fGb of space\n",(1. * kmers) / (0x20000000/sizeof(KmerPos)));
      fflush(stdout);
    }

  //  Build the k-mer list

  { int i;

    FR_src = src;

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,tuple_thread,parmt+i);
    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);
  }

  //  Sort the k-mer list

  { int i;
    int mersort[11];

#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
    for (i = 0; i < (Kmer-1)/4+1; i++)
      mersort[i] = 8+i;
#else
    for (i = 0; i < (Kmer-1)/4+1; i++)
      mersort[i] = 17-i;
#endif
    mersort[i] = -1;

    rez = (KmerPos *) LSD_Sort(kmers,src,trg,16,16,mersort);
  }

  //  Compress frequent tuples if requested

  if (TooFrequent < INT32_MAX && kmers > 0)
    { int    i, x, z;
      uint64 h;

      parmt[0].beg = 0;
      for (i = 1; i < NTHREADS; i++)
        { x = (((int64) i)*kmers) / NTHREADS;
          h = rez[x-1].code;
          while (rez[x].code == h)
            x += 1;
          parmt[i-1].end = parmt[i].beg = x;
        }
      parmt[NTHREADS-1].end = kmers;

      if (rez[kmers-1].code == MAX_CODE_64)
        rez[kmers].code = 0;
      else
        rez[kmers].code = MAX_CODE_64;

      if (src == rez)
        { FR_src = src;
          FR_trg = rez = trg;
        }
      else
        { FR_src = trg;
          FR_trg = rez = src;
        }

      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,compsize_thread,parmt+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);

      x = 0;
      for (i = 0; i < NTHREADS; i++)
        { z = parmt[i].fill;
          parmt[i].fill = x;
          x += z;
        }
      kmers = x;

      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,compress_thread,parmt+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);
    }

  rez[kmers].code   = MAX_CODE_64;
  rez[kmers+1].code = 0;
    
  if (src != rez)
    free(src);
  else
    free(trg);

#ifdef TEST_KSORT
  { int i;

    printf("\nKMER SORT:\n");
    for (i = 0 /*100000000*/; i < 100000000+HOW_MANY && i < kmers; i++)
      { KmerPos *c = rez+i;
        printf(" %9d:  %6d%c / %6d / %016llx\n",i,c->read>>1,
                                               (c->read&0x1)?'c':'n',(c->rpos & POST_MASK),c->code);
      }
    fflush(stdout);
  }
#endif

#ifdef HISTOGRAM_KSORT
  { int    hist[100];
    uint64 ca;
    int    i, j;

    for (i = 0; i < 100; i++)
      hist[i] = 0;

    i = 0;
    while (i < kmers)
      { ca = rez[i].code;
        j = i++;
        while (rez[i].code == ca)
          i += 1;
        if (i-j >= 100)
          hist[99] += 1;
        else
          hist[i-j] += 1;
      }

     for (i = 99; i >= 0; i--)
       printf(" %2d: %6d\n",i,hist[i]);
  }
#endif

  if (VERBOSE)
    { if (TooFrequent < INT32_MAX)
        { printf("   Revised kmer count = ");
          Print_Number((int64) kmers,0,stdout);
          printf("\n");
        }
      printf("   Index occupies %.2fGb\n",(1. * kmers) / (0x40000000/sizeof(KmerPos)));
      fflush(stdout);
    }

  if (kmers <= 0)
    { free(rez);
      goto no_mers;
    }

  if (kmers > (int64) (MEM_LIMIT/(4*sizeof(KmerPos))))
    { fprintf(stderr,"Warning: Block size too big, index occupies more than 1/4 of");
      if (MEM_LIMIT == MEM_PHYSICAL)
        fprintf(stderr," physical memory (%.1fGb)\n",(1.*MEM_LIMIT)/0x40000000ll);
      else
        fprintf(stderr," desired memory allocation (%.1fGb)\n",(1.*MEM_LIMIT)/0x40000000ll);
      fflush(stderr);
    }

  *len = kmers;
  return (rez);

no_mers:
  *len = 0;
  return (NULL);
}


/*******************************************************************************************
 *
 *  FILTER MATCH
 *
 ********************************************************************************************/

static int find_tuple(uint64 x, KmerPos *a, int n)
{ int l, r, m;

  // smallest k s.t. a[k].code >= x (or n if does not exist)

  l = 0;
  r = n;
  while (l < r)
    { m = ((l+r) >> 1);
      if (a[m].code < x)
        l = m+1;
      else
        r = m;
    }
  return (l);
}

  //  Determine what *will* be the size of the merged list and histogram of sizes for given cutoffs

static KmerPos  *MG_alist;
static KmerPos  *MG_blist;
static DAZZ_DB  *MG_ablock;
static DAZZ_DB  *MG_bblock;
static SeedPair *MG_hits;
static int       MG_self;

typedef struct
  { int    abeg, aend;
    int    bbeg, bend;
    int64  nhits;
    int    limit;
    int64  hitgram[MAXGRAM];
  } Merge_Arg;

static void *count_thread(void *arg)
{ Merge_Arg  *data   = (Merge_Arg *) arg;
  KmerPos    *asort  = MG_alist;
  KmerPos    *bsort  = MG_blist;
  int64      *gram   = data->hitgram;
  int64       nhits  = 0;
  int         aend   = data->aend;

  int64  ct;
  int    ia, ja;
  uint64 ca, da;

  ia = data->abeg;
  ca = asort[ia].code;
  if (MG_self)
    { uint32 ar;
      int    ka;

      while (1)
        { ja = ka = ia++;
          ct = 0;
          if (IDENTITY)
            while (1)
              { da = asort[ia].code;
                if (da != ca)
                  break;
                ct += (ia-ja);
                ia += 1;
              }
          else
            while (1)
              { da = asort[ia].code;
                if (da != ca)
                  break;
                ar = (asort[ia].read & ~0x1u);
                while (ka < ia && asort[ka].read < ar)
                  ka += 1;
                ct += (ka-ja);
                ia += 1;
              }

          ca = da;
          if (ia > aend)
            { if (ja >= aend)
                break;
              ia = aend;
              ca = asort[ia].code;
              ct -= (ka-ja);
            }

          nhits += ct;
          if (ct < MAXGRAM)
            gram[ct] += 1;
        }
    }
  else
    { int    ib, jb;
      uint64 cb;

      ib = data->bbeg;
      cb = bsort[ib].code;
      while (1)
        { ja = ia++;
          while (1)
            { da = asort[ia].code;
              if (da != ca)
                break;
              ia += 1;
            }

          if (ia > aend)
            { if (ja >= aend)
                break;
              ia  = aend;
              da = asort[ia].code;
            }

          while (cb < ca)
            { ib += 1;
              cb = bsort[ib].code;
            }
          if (cb != ca)
            { ca = da;
              continue;
            }

          jb = ib++;
          while (1)
            { cb = bsort[ib].code;
              if (cb != ca)
                break;
              ib += 1;
            }
          ca = da;

          ct = ((int64) (ia-ja))*(ib-jb);
          nhits += ct;
          if (ct < MAXGRAM)
            gram[ct] += 1;
        }
    }

  data->nhits = nhits;

  return (NULL);
}

  //  Produce the merged list now that the list has been allocated and
  //    the appropriate cutoff determined.

static void *merge_thread(void *arg)
{ Merge_Arg  *data   = (Merge_Arg *) arg;
  KmerPos    *asort  = MG_alist;
  KmerPos    *bsort  = MG_blist;
  DAZZ_READ  *reads  = MG_bblock->reads;
  SeedPair   *hits   = MG_hits;
  int64       nhits  = data->nhits;
  int         aend   = data->aend;
  int         limit  = data->limit;

  int64  ct;
  int    ia, ja;
  uint64 ca, da;
  int    nread = MG_ablock->nreads;

  ia = data->abeg;
  ca = asort[ia].code;
  if (MG_self)
    { uint32 ar, br;
      uint32 ap, bp;
      uint32 as, bs;
      int    a, ka;

      while (1)
        { ja = ka = ia++;
          ct = 0;
          if (IDENTITY)
            while (1)
              { da = asort[ia].code;
                if (da != ca)
                  break;
                ct += (ia-ja);
                ia += 1;
              }
          else
            while (1)
              { da = asort[ia].code;
                if (da != ca)
                  break;
                ar = (asort[ia].read & ~0x1u);
                while (ka < ia && asort[ka].read < ar)
                  ka += 1;
                ct += (ka-ja);
                ia += 1;
              }

          ca = da;
          if (ia > aend)
            { if (ja >= aend)
                break;
              ia = aend;
              ca = asort[ia].code;
              ct -= (ka-ja);
            }

          if (ct >= limit)
            continue;

          if (IDENTITY)
            for (ka = ja+1; ka < ia; ka++)
              { ar = asort[ka].read;
                as = (ar & SIGN_BIT);
                ar >>= 1;
                ap = (asort[ka].rpos & POST_MASK);
                for (a = ja; a < ka; a++)
                  { br = asort[a].read;
                    bs = (br & SIGN_BIT);
                    br >>= 1;
                    bp = asort[a].rpos;
                    if (bs == as)
                      { bp = (bp & POST_MASK);
                        hits[nhits].aread = ar;
                      }
                    else
                      { if ((bp & LONG_BIT) != 0)
                          bp = (reads[br].rlen - (bp & POST_MASK)) + Koff;
                        else
                          bp = (reads[br].rlen - (bp & POST_MASK)) + Kmer;
                        hits[nhits].aread = ar + nread;
                      }
                    hits[nhits].bread = br;
                    hits[nhits].apos  = ap; 
                    hits[nhits].diag  = ap - bp;
                    nhits += 1;
                  }
              }
          else
            for (ka = ja+1; ka < ia; ka++)
              { ar = asort[ka].read;
                as = (ar & SIGN_BIT);
                ar >>= 1;
                ap = (asort[ka].rpos & POST_MASK);
                for (a = ja; a < ka; a++)
                  { br = asort[a].read;
                    bs = (br & SIGN_BIT);
                    br >>= 1;
                    if (br >= ar)
                      break;
                    bp = asort[a].rpos;
                    if (bs == as)
                      { bp = (bp & POST_MASK);
                        hits[nhits].aread = ar;
                      }
                    else
                      { if ((bp & LONG_BIT) != 0)
                          bp = (reads[br].rlen - (bp & POST_MASK)) + Koff;
                        else
                          bp = (reads[br].rlen - (bp & POST_MASK)) + Kmer;
                        hits[nhits].aread = ar + nread;
                      }
                    hits[nhits].bread = br;
                    hits[nhits].apos  = ap; 
                    hits[nhits].diag  = ap - bp;
                    nhits += 1;
                  }
              }
        }
    }
  else
    { int    ib, jb;
      uint64 cb;
      uint32 ar, br;
      uint32 ap, bp;
      uint32 as, bs;
      int    a, b;

      ib = data->bbeg;
      cb = bsort[ib].code;
      while (1)
        { ja = ia++;
          while (1)
            { da = asort[ia].code;
              if (da != ca)
                break;
              ia += 1;
            }

          if (ia > aend)
            { if (ja >= aend)
                break;
              ia = aend;
              da = asort[ia].code;
            }
          
          while (cb < ca)
            { ib += 1;
              cb = bsort[ib].code;
            }
          if (cb != ca)
            { ca = da;
              continue;
            }
          
          jb = ib++;
          while (1) 
            { cb = bsort[ib].code;
              if (cb != ca)
                break;
              ib += 1;
            } 
          ca = da;

          if (((int64) (ia-ja))*(ib-jb) >= limit)
            continue;

          for (a = ja; a < ia; a++)
            { ar = asort[a].read;
              as = (ar & SIGN_BIT);
              ar >>= 1;
              ap = (asort[a].rpos & POST_MASK);
              for (b = jb; b < ib; b++)
                { br = bsort[b].read;
                  bs = (br & SIGN_BIT);
                  br >>= 1;
                  bp = bsort[b].rpos;
                  if (bs == as)
                    { bp = (bp & POST_MASK);
                      hits[nhits].aread = ar;
                    }
                  else
                    { if ((bp & LONG_BIT) != 0)
                        bp = (reads[br].rlen - (bp & POST_MASK)) + Koff;
                      else
                        bp = (reads[br].rlen - (bp & POST_MASK)) + Kmer;
                      hits[nhits].aread = ar + nread;
                    }
                  hits[nhits].bread = br;
                  hits[nhits].apos  = ap; 
                  hits[nhits].diag  = ap - bp;
                  nhits += 1;
                }
            }

        }
    }

  return (NULL);
}

  //  Report threads: given a segment of merged list, find all seeds and from them all alignments.

static DAZZ_DB    *MR_ablock;
static DAZZ_DB    *MR_bblock;
static SeedPair   *MR_hits;
static int         MR_two;
static Align_Spec *MR_spec;
static int         MR_tspace;

typedef struct
  { uint64   max;
    uint64   top;
    uint16  *trace;
  } Trace_Buffer;

#ifdef DO_BRIDGING

static inline int MapToTPAbove(Path *path, int *x, int isA, Trace_Buffer *tbuf)
{ uint16 *trace = tbuf->trace + (uint64) path->trace;
  int a, b, i;

  a = (path->abpos / MR_tspace) * MR_tspace;
  b = path->bbpos;
  for (i = 1; i < path->tlen; i += 2)
    { a += MR_tspace;
      b += trace[i];
      if (a > path->aepos)
        a = path->aepos;
      if (isA)
        { if (a >= *x)
            { *x = a;
              return (b);
            }
        }
      else
        { if (b >= *x)
            { *x = b;
              return (a);
            }
        }
    }
  if (isA)
    { *x = a;
      return (b);
    }
  else
    { *x = b;
      return (a);
    }
}

static inline int MapToTPBelow(Path *path, int *x, int isA, Trace_Buffer *tbuf)
{ uint16 *trace = tbuf->trace + (uint64) path->trace;
  int a, b, i;

  a = ((path->aepos + (MR_tspace-1)) / MR_tspace) * MR_tspace;
  b = path->bepos;
  for (i = path->tlen-1; i >= 0; i -= 2)
    { a -= MR_tspace;
      b -= trace[i];
      if (a < path->abpos)
        a = path->abpos;
      if (isA)
        { if (a <= *x)
            { *x = a;
              return (b);
            }
        }
      else
        { if (b <= *x)
            { *x = b;
              return (a);
            }
        }
    }
  if (isA)
    { *x = a;
      return (b);
    }
  else
    { *x = b;
      return (a);
    }
}

static int Check_Bridge(Path *path)
{ uint16 *trace = (uint16 *) path->trace;
  int     i;

  if (MR_tspace <= TRACE_XOVR)
    { for (i = 0; i < path->tlen; i++)
        if (trace[i] > 250)
          return (1);
    }
  return (0);
}

static void Compute_Bridge_Path(Path *path1, Path *path2, Alignment *align, int comp,
                                int aovl, int bovl, Work_Data *work, Trace_Buffer *tbuf)
{ Path   *apath;
  int     ain, aout;
  int     bin, bout, boff;
  int     err;
  int     i, j, p;
  uint16 *trk;

  apath = align->path;

  if (bovl > aovl)
    { bin  = path2->bbpos;
      bout = path1->bepos;
      ain  = MapToTPBelow(path1,&bin,0,tbuf);
      aout = MapToTPAbove(path2,&bout,0,tbuf);
    }
  else
    { ain  = path2->abpos;
      aout = path1->aepos;
      bin  = MapToTPBelow(path1,&ain,1,tbuf);
      bout = MapToTPAbove(path2,&aout,1,tbuf);
    }

#ifdef TEST_BRIDGE
  printf("\n  Tangle [%5d..%5d] vs [%5d..%5d]  %4d\n",
      path1->abpos,path1->aepos,path2->abpos,path2->aepos,abs(aovl-bovl));
  printf("         [%5d..%5d] vs [%5d..%5d]  %4d vs %4d\n",
      path1->bbpos,path1->bepos,path2->bbpos,path2->bepos,aovl,bovl);
  printf("      (%d,%d) to (%d,%d)\n",ain,bin,aout,bout);
  fflush(stdout);
#endif

  apath->abpos = ain - 2*MR_tspace;
  apath->aepos = aout + 2*MR_tspace;
  apath->bbpos = MapToTPBelow(path1,&(apath->abpos),1,tbuf);
  apath->bepos = MapToTPAbove(path2,&(apath->aepos),1,tbuf);

  if (comp)
    { boff = MR_tspace - apath->aepos % MR_tspace;

      p = align->alen - apath->abpos;
      apath->abpos = align->alen - apath->aepos;
      apath->aepos = p;
      p = align->blen - apath->bbpos;
      apath->bbpos = align->blen - apath->bepos;
      apath->bepos = p;

      boff = boff - apath->abpos % MR_tspace;
      align->aseq  -= boff;
      apath->abpos += boff;
      apath->aepos += boff;
      align->alen  += boff;
    }

#ifdef TEST_BRIDGE
  printf("\n      (%d,%d) to (%d,%d)\n",apath->abpos,apath->bbpos,apath->aepos,apath->bepos);
  fflush(stdout);

  Compute_Alignment(align,work,DIFF_ALIGN,0);
  Print_Reference(stdout,align,work,8,100,10,0,6);
  fflush(stdout);
#endif

  Compute_Alignment(align,work,DIFF_TRACE,MR_tspace);

  trk = (uint16 *) apath->trace;
  if (comp)
    { j = apath->tlen-2;
      i = 0;
      while (i < j)
        { p = trk[i];
          trk[i] = trk[j];
          trk[j] = p;
          p = trk[i+1];
          trk[i+1] = trk[j+1];
          trk[j+1] = p;
          i += 2;
          j -= 2;
        }

      align->aseq  += boff;
      apath->abpos -= boff;
      apath->aepos -= boff;
      align->alen  -= boff;

      p = align->alen - apath->abpos;
      apath->abpos = align->alen - apath->aepos;
      apath->aepos = p;
      p = align->blen - apath->bbpos;
      apath->bbpos = align->blen - apath->bepos;
      apath->bepos = p;
    }

  bin  = apath->bbpos;
  bout = apath->bepos;
  err  = apath->diffs;

  p = 2*(ain / MR_tspace - apath->abpos / MR_tspace);
  for (i = 0; i < p; i += 2)
    { bin += trk[i+1];
      err -= trk[i];
    }

  p = 2*(apath->aepos / MR_tspace - aout / MR_tspace);
  for (i = align->path->tlen, p = i-p; i > p; i -= 2)
    { bout -= trk[i-1];
      err  -= trk[i-2];
    }

#ifdef TEST_BRIDGE
  printf("      (%d,%d) to (%d,%d)\n",ain,bin,aout,bout);
  printf("  Box %d vs %d -> %d %d%%\n",aout-ain,bout-bin,err,
                                       (200*err)/((aout-ain)+(bout-bin)));
  fflush(stdout);
#endif
}

#endif // DO_BRIDGING

static int Entwine(Path *jpath, Path *kpath, Trace_Buffer *tbuf, int *where)
{ int   ac, b2, y2, ae;
  int   i, j, k;
  int   num, den, min;
#ifdef SEE_ENTWINE
  int   strt = 1;
  int   iflare, oflare;
#endif

  uint16 *ktrace = tbuf->trace + (uint64) (kpath->trace);
  uint16 *jtrace = tbuf->trace + (uint64) (jpath->trace);

  min   = 10000;
  num   = 0;
  den   = 0;

#ifdef SEE_ENTWINE
  printf("\n");
#endif

  y2 = jpath->bbpos;
  j  = jpath->abpos/MR_tspace;

  b2 = kpath->bbpos;
  k  = kpath->abpos/MR_tspace;

  if (jpath->abpos == kpath->abpos)
    { min = abs(y2-b2);
      if (min == 0)
        *where = kpath->abpos;
    }

  if (j < k)
    { ac = k*MR_tspace;

      j = 1 + 2*(k-j);
      k = 1;

      for (i = 1; i < j; i += 2)
        y2 += jtrace[i];
    }
  else
    { ac = j*MR_tspace;

      k = 1 + 2*(j-k);
      j = 1;

      for (i = 1; i < k; i += 2)
        b2 += ktrace[i];
    }

  ae = jpath->aepos;
  if (ae > kpath->aepos)
    ae = kpath->aepos;

  while (1)
    { ac += MR_tspace;
      if (ac >= ae)
        break;
      y2 += jtrace[j];
      b2 += ktrace[k];
      j += 2;
      k += 2;

#ifdef SEE_ENTWINE
      printf("   @ %5d : %5d %5d = %4d\n",ac,y2,b2,abs(b2-y2));
#endif

      i = abs(y2-b2);
      if (i <= min)
        { min = i;
          if (i == 0)
            *where = ac;
        }
      num += i;
      den += 1;
#ifdef SEE_ENTWINE
      if (strt)
        { strt   = 0;
          iflare = i;
        }
      oflare = i;
#endif
    }

  if (jpath->aepos == kpath->aepos)
    { i = abs(jpath->bepos-kpath->bepos);
      if (i <= min)
        { min = i;
          if (i == 0)
            *where = kpath->aepos;
        }
    }

#ifdef SEE_ENTWINE
  if (den == 0)
    printf("Nothing\n");
  else
    printf("MINIM = %d AVERAGE = %d  IFLARE = %d  OFLARE = %d\n",min,num/den,iflare,oflare);
#endif

  if (den == 0)
    return (-1);
  else
    return (min);
}


//  Produce the concatentation of path1 and path2 where they are known to meet at
//    the trace point with coordinate ap. Place this result in a big growing buffer,
//    that gets reset when fusion is called with path1 = NULL

static void Fusion(Path *path1, int ap, Path *path2, Trace_Buffer *tbuf)
{ int     k, k1, k2;
  int     len, diff;
  uint16 *trace;

  k1 = 2 * ((ap/MR_tspace) - (path1->abpos/MR_tspace));
  k2 = 2 * ((ap/MR_tspace) - (path2->abpos/MR_tspace));

  len = k1+(path2->tlen-k2);

  if (tbuf->top + len >= tbuf->max)
    { tbuf->max = 1.2*(tbuf->top+len) + 1000;
      tbuf->trace = (uint16 *) Realloc(tbuf->trace,sizeof(uint16)*tbuf->max,"Allocating paths");
      if (tbuf->trace == NULL)
        Clean_Exit(1);
    }

  trace = tbuf->trace + tbuf->top;
  tbuf->top += len;

  diff = 0;
  len  = 0;
  if (k1 > 0)
    { uint16 *t = tbuf->trace + (uint64) (path1->trace);
      for (k = 0; k < k1; k += 2)
        { trace[len++] = t[k];
          trace[len++] = t[k+1];
          diff += t[k];
        }
    }
  if (k2 < path2->tlen)
    { uint16 *t = tbuf->trace + (uint64) (path2->trace);
      for (k = k2; k < path2->tlen; k += 2)
        { trace[len++] = t[k];
          trace[len++] = t[k+1];
          diff += t[k];
        }
    }

  path1->aepos = path2->aepos;
  path1->bepos = path2->bepos;
  path1->diffs = diff;
  path1->trace = (void *) (trace - tbuf->trace);
  path1->tlen  = len;
}

#ifdef DO_BRIDGING

//  Produce the concatentation of path1, path2, and path3 where they are known to meet at
//    the ends of path2 which was produced by Compute-Alignment. Place this result in
//    a big growing buffer.

static void Bridge(Path *path1, Path *path2, Path *path3, Trace_Buffer *tbuf)
{ int     k, k1, k2;
  int     len, diff;
  uint16 *trace;

  k1 = 2 * ((path2->abpos/MR_tspace) - (path1->abpos/MR_tspace));
  if (path2->aepos == path3->aepos)
    k2 = path3->tlen;
  else
    k2 = 2 * ((path2->aepos/MR_tspace) - (path3->abpos/MR_tspace));

  len = k1 + path2->tlen + (path3->tlen-k2);

  if (tbuf->top + len >= tbuf->max)
    { tbuf->max = 1.2*(tbuf->top+len) + 1000;
      tbuf->trace = (uint16 *) Realloc(tbuf->trace,sizeof(uint16)*tbuf->max,"Allocating paths");
      if (tbuf->trace == NULL)
        Clean_Exit(1);
    }

  trace = tbuf->trace + tbuf->top;
  tbuf->top += len;

  diff = 0;
  len  = 0;
  if (k1 > 0)
    { uint16 *t = tbuf->trace + (uint64) (path1->trace);
      for (k = 0; k < k1; k += 2)
        { trace[len++] = t[k];
          trace[len++] = t[k+1];
          diff += t[k];
        }
    }
  if (path2->tlen > 0)
    { uint16 *t = (uint16 *) (path2->trace);
      for (k = 0; k < path2->tlen; k += 2)
        { trace[len++] = t[k];
          trace[len++] = t[k+1];
          diff += t[k];
        }
    }
  if (k2 < path3->tlen)
    { uint16 *t = tbuf->trace + (uint64) (path3->trace);
      for (k = k2; k < path3->tlen; k += 2)
        { trace[len++] = t[k];
          trace[len++] = t[k+1];
          diff += t[k];
        }
    }

  path1->aepos = path3->aepos;
  path1->bepos = path3->bepos;
  path1->diffs = diff;
  path1->trace = (void *) (trace - tbuf->trace);
  path1->tlen  = len;
}

#endif // DO_BRIDGING

static int Handle_Redundancies(Path *amatch, int novls, Path *bmatch,
                               Alignment *align, Work_Data *work, Trace_Buffer *tbuf)
{ Path      *jpath, *kpath, *apath;
  Path      _bpath, *bpath = &_bpath;
  Alignment _blign, *blign = &_blign;

  int   j, k, no;
  int   dist;
  int   awhen = 0, bwhen = 0;
  int   hasB, comp;

#if defined(TEST_CONTAIN) || defined(TEST_BRIDGE)
  for (j = 0; j < novls; j++)
    printf("  %3d: [%5d,%5d] x [%5d,%5d]\n",j,amatch[j].abpos,amatch[j].aepos,
                                              amatch[j].bbpos,amatch[j].bepos);
#endif

  (void) work;

  //  Loop to catch LA's that share a common trace point and fuse them

  hasB = (bmatch != NULL);
  if (hasB)
    { blign->aseq = align->bseq;
      blign->bseq = align->aseq;
      blign->alen = align->blen;
      blign->blen = align->alen;
      blign->path = bpath;
    }
  apath = align->path;
  comp  = COMP(align->flags);

  (void) apath->tlen;   //  Just to shut up stupid compilers

  for (j = 1; j < novls; j++)
    { jpath = amatch+j;
      for (k = j-1; k >= 0; k--)
        { kpath = amatch+k;

          if (kpath->abpos < 0)
            continue;

          if (jpath->abpos < kpath->abpos)

            { if (kpath->abpos <= jpath->aepos && kpath->bbpos <= jpath->bepos)
                { dist = Entwine(jpath,kpath,tbuf,&awhen);
                  if (dist == 0)
                    { if (kpath->aepos > jpath->aepos)
                        { if (hasB)
                            { if (comp)
                                { dist = Entwine(bmatch+k,bmatch+j,tbuf,&bwhen);
                                  if (dist != 0)
                                    continue;
                                  Fusion(jpath,awhen,kpath,tbuf);
                                  Fusion(bmatch+k,bwhen,bmatch+j,tbuf);
                                  bmatch[j] = bmatch[k];
#ifdef TEST_CONTAIN
                                  printf("  Really 1");
#endif
                                }
                              else
                                { dist = Entwine(bmatch+j,bmatch+k,tbuf,&bwhen);
                                  if (dist != 0)
                                    continue;
                                  Fusion(jpath,awhen,kpath,tbuf);
                                  Fusion(bmatch+j,bwhen,bmatch+k,tbuf);
#ifdef TEST_CONTAIN
                                  printf("  Really 2");
#endif
                                }
                            }
                          else
                            { Fusion(jpath,awhen,kpath,tbuf);
#ifdef TEST_CONTAIN
                              printf("  Really 3");
#endif
                            }
                          k = j;
                        }
                      kpath->abpos = -1;
#ifdef TEST_CONTAIN
                      printf("  Fuse! A %d %d\n",j,k);
#endif
                    }
                }
            }

          else // kpath->abpos <= jpath->abpos

            { if (jpath->abpos <= kpath->aepos && jpath->bbpos <= kpath->bepos)
                { dist = Entwine(kpath,jpath,tbuf,&awhen);
                  if (dist == 0)
                    { if (kpath->abpos == jpath->abpos)
                        { if (kpath->aepos > jpath->aepos)
                            { *jpath = *kpath;
                              if (hasB)
                                bmatch[j] = bmatch[k];
                            }
                        }
                      else if (jpath->aepos > kpath->aepos)
                        { if (hasB)
                            { if (comp)
                                { dist = Entwine(bmatch+j,bmatch+k,tbuf,&bwhen);
                                  if (dist != 0)
                                    continue;
                                  Fusion(kpath,awhen,jpath,tbuf);
                                  *jpath = *kpath;
                                  Fusion(bmatch+j,bwhen,bmatch+k,tbuf);
#ifdef TEST_CONTAIN
                                  printf("  Really 4");
#endif
                                }
                              else
                                { dist = Entwine(bmatch+k,bmatch+j,tbuf,&bwhen);
                                  if (dist != 0)
                                    continue;
                                  Fusion(kpath,awhen,jpath,tbuf);
                                  *jpath = *kpath;
                                  Fusion(bmatch+k,bwhen,bmatch+j,tbuf);
                                  bmatch[j] = bmatch[k];
#ifdef TEST_CONTAIN
                                  printf("  Really 5");
#endif
                                }
                            }
                          else
                            { Fusion(kpath,awhen,jpath,tbuf);
                              *jpath = *kpath;
#ifdef TEST_CONTAIN
                              printf("  Really 6");
#endif
                            }
                          k = j;
                        }
                      else
                        { *jpath = *kpath;
                          if (hasB)
                            bmatch[j] = bmatch[k];
                        }
                      kpath->abpos = -1;
#ifdef TEST_CONTAIN
                      printf("  Fuse! B %d %d\n",j,k);
#endif
                    }
                }
            }
        }
    }

#ifdef DO_BRIDGING

  //  Loop to catch LA's that have a narrow parallel overlap and bridge them

  for (j = 1; j < novls; j++)
    { jpath = amatch+j;
      if (jpath->abpos < 0)
        continue;

      for (k = j-1; k >= 0; k--)
        { Path   *path1, *path2;
          Path   *bath1, *bath2;
          int     aovl, bovl;
          Path    jback, kback;

          kpath = amatch+k;
	  if (kpath->abpos < 0)
            continue;

          if (jpath->abpos < kpath->abpos)
            { path1 = jpath;
              path2 = kpath;
            }
          else
            { path1 = kpath;
              path2 = jpath;
            }

          if (path2->abpos >= path1->aepos || path1->aepos >= path2->aepos ||
              path1->bbpos >= path2->bbpos || path2->bbpos >= path1->bepos ||
              path1->bepos >= path2->bepos)
            continue;
          aovl = path1->aepos - path2->abpos;
          bovl = path1->bepos - path2->bbpos;
          if (abs(aovl-bovl) > .2 * (aovl+bovl))
            continue;

          if (hasB)
            { if (comp == (jpath->abpos < kpath->abpos))
                { bath1 = bmatch+k;
                  bath2 = bmatch+j;
                }
              else
                { bath1 = bmatch+j;
                  bath2 = bmatch+k;
                }
              if (bath1->abpos > bath2->abpos)
                { printf("  SYMFAIL %d %d\n",j,k);
                  continue;
                }
            }

          Compute_Bridge_Path(path1,path2,align,0,aovl,bovl,work,tbuf);

          if (Check_Bridge(apath))
            continue;

          jback = *jpath;
          kback = *kpath;

          Bridge(path1,apath,path2,tbuf);
          *jpath = *path1;
          kpath->abpos = -1;

#ifdef TEST_BRIDGE
          { Alignment extra;
            Path      pcopy;

            pcopy  = *jpath;
            extra  = *align;
            pcopy.trace = tbuf->trace + (uint64) jpath->trace;
            extra.path = &pcopy;
            Compute_Trace_PTS(&extra,work,MR_tspace,GREEDIEST);
            Print_Reference(stdout,&extra,work,8,100,10,0,6);
            fflush(stdout);
          }
#endif

          if (hasB)
            { Compute_Bridge_Path(bath1,bath2,blign,comp,bovl,aovl,work,tbuf);

              if (Check_Bridge(bpath))
                { *jpath = jback;
                  *kpath = kback;
                  continue;
                }
             
              Bridge(bath1,bpath,bath2,tbuf);
              bmatch[j] = *bath1;

#ifdef TEST_BRIDGE
              { Alignment extra;
                Path      pcopy;

	        pcopy  = bmatch[j];
	        extra  = *blign;
                pcopy.trace = tbuf->trace + (uint64) bmatch[j].trace;
                extra.path = &pcopy;
                if (comp)
                  { Complement_Seq(extra.aseq,extra.alen);
                    Complement_Seq(extra.bseq,extra.blen);
                  }
                Compute_Trace_PTS(&extra,work,MR_tspace,GREEDIEST);
                Print_Reference(stdout,&extra,work,8,100,10,0,6);
                fflush(stdout);
                if (comp)
                  { Complement_Seq(extra.aseq,extra.alen);
                    Complement_Seq(extra.bseq,extra.blen);
                  }
              }
#endif
            }
        } 
    }

#endif // DO_BRIDGING

  no = 0;
  for (j = 0; j < novls; j++)
    if (amatch[j].abpos >= 0)
      { if (hasB)
          bmatch[no] = bmatch[j];
        amatch[no++] = amatch[j];
      }
  novls = no;

#if defined(TEST_CONTAIN) || defined(TEST_BRIDGE)
  for (j = 0; j < novls; j++)
    printf("  %3d: [%5d,%5d] x [%5d,%5d]\n",j,amatch[j].abpos,amatch[j].aepos,
                                              amatch[j].bbpos,amatch[j].bepos);
#endif

  return (novls);
}

static void Diagonal_Span(Path *path, int *mind, int *maxd)
{ uint16 *points;
  int     i, tlen;
  int     dd, low, hgh;

  points = path->trace;
  tlen   = path->tlen;

  dd = path->abpos - path->bbpos;
  low = hgh = dd;

  dd = path->aepos - path->bepos;
  if (dd < low)
    low = dd;
  else if (dd > hgh)
    hgh = dd;

  dd = (path->abpos/MR_tspace)*MR_tspace - path->bbpos;
  tlen -= 2;
  for (i = 1; i < tlen; i += 2)
    { dd += MR_tspace - points[i];
      if (dd < low)
        low = dd;
      else if (dd > hgh)
        hgh = dd;
    }

  *mind = (low >> Binshift)-1;
  *maxd = (hgh >> Binshift)+1;
}

static void CopyAndComp(char *bcomp,char *bseq, int blen)
{ char *s, *t;

  t = bcomp + (blen-1);
  s = bseq;
  t[1] = 4;
  while (t >= bcomp)
    *t-- = 3-*s++;
  t[0] = 4;
}

typedef struct
  { int64       beg, end;
    int        *score;
    int        *lastp;
    int        *lasta;
    Work_Data  *work;
    FILE       *ofile1;
    FILE       *ofile2;
    int64       nfilt;
    int64       nlas;
#ifdef PROFILE
    int         profyes[MAXHIT+1];
    int         profno[MAXHIT+1];
#endif
  } Report_Arg;

typedef struct
  { uint64 p1;   //  The lower half
    uint64 p2;
  } Double;

static void *report_thread(void *arg)
{ Report_Arg  *data   = (Report_Arg *) arg;
  SeedPair    *hits   = MR_hits;
  Double *hitd = (Double *) MR_hits;
  DAZZ_READ   *bread  = MR_bblock->reads;
  DAZZ_READ   *aread  = MR_ablock->reads;
  char        *aseq   = (char *) (MR_ablock->bases);
  char        *bseq   = (char *) (MR_bblock->bases);
  int         *score  = data->score;
  int         *scorp  = data->score + 1;
  int         *scorm  = data->score - 1;
  int         *lastp  = data->lastp;
  int         *lasta  = data->lasta;
  int          afirst = MR_ablock->tfirst;
  int          bfirst = MR_bblock->tfirst;
  FILE        *ofile1 = data->ofile1;
  FILE        *ofile2 = data->ofile2;
  Work_Data   *work   = data->work;
  int          maxdiag = ( MR_ablock->maxlen >> Binshift);
  int          mindiag = (-MR_bblock->maxlen >> Binshift);
  int          areads  = MR_ablock->nreads;
#ifdef PROFILE
  int         *profyes = data->profyes;
  int         *profno  = data->profno;
  int          maxhit, isyes;
#endif

  Overlap     _ovlb, *ovlb = &_ovlb;
  Overlap     _ovla, *ovla = &_ovla;
  Alignment   _align, *align = &_align;
  Path        *apath = &(ovla->path);
  Path        *bpath;
  char        *bcomp;

  int    AOmax, BOmax;
  int    novla, novlb;
  Path  *amatch, *bmatch;

  Trace_Buffer _tbuf, *tbuf = &_tbuf;
  int          small, tbytes;

  Double *hitc;
  int     minhit;
  uint64  cpair;
  uint64  npair = 0;
  int64   nidx, eidx;

  int64 nfilt = 0;
  int64 nlas  = 0;
  int64 ahits = 0;
  int64 bhits = 0;

  //  In ovl and align roles of A and B are reversed, as the B sequence must be the
  //    complemented sequence !!

  align->path = apath;
  bcomp = New_Read_Buffer(MR_bblock);

  if (MR_tspace <= TRACE_XOVR)
    { small  = 1;
      tbytes = sizeof(uint8);
    }
  else
    { small  = 0;
      tbytes = sizeof(uint16);
    }

  AOmax = BOmax = MATCH_CHUNK;
  amatch = Malloc(sizeof(Path)*AOmax,"Allocating match vector");
  bmatch = Malloc(sizeof(Path)*BOmax,"Allocating match vector");

  tbuf->max   = 2*TRACE_CHUNK;
  tbuf->trace = Malloc(sizeof(short)*tbuf->max,"Allocating trace vector");

  if (amatch == NULL || bmatch == NULL || tbuf->trace == NULL)
    Clean_Exit(1);

  fwrite(&ahits,sizeof(int64),1,ofile1);
  fwrite(&MR_tspace,sizeof(int),1,ofile1);
  if (MR_two)
    { fwrite(&bhits,sizeof(int64),1,ofile2);
      fwrite(&MR_tspace,sizeof(int),1,ofile2);
    }

#ifdef PROFILE
  { int i;
    for (i = 0; i <= MAXHIT; i++)
      profyes[i] = profno[i] = 0;
  }
#endif

  minhit = (Hitmin-1)/Kmer + 1;
  hitc   = hitd + (minhit-1);
  eidx   = data->end - minhit;
  nidx   = data->beg;
  for (cpair = hitd[nidx].p1; nidx <= eidx; cpair = npair)
    if (hitc[nidx].p1 != cpair)
      { nidx += 1;
        while ((npair = hitd[nidx].p1) == cpair)
          nidx += 1;
      }
    else
      { int   ar, br, bc;
        int   alen, blen;
        int   doA, doB;
        int   setaln, amark, amark2;
        int   apos, bpos, diag;
        int64 lidx, sidx;
        int64 f, h2;

        ar = hits[nidx].aread;
        br = hits[nidx].bread;
        if (ar >= areads)
          { bc = 1;
            ar -= areads;
          }
        else
          bc = 0;
        alen = aread[ar].rlen;
        blen = bread[br].rlen;
        if (alen < HGAP_MIN && blen < HGAP_MIN)
          { nidx += 1;
            while ((npair = hitd[nidx].p1) == cpair)
              nidx += 1;
            continue;
          }
#ifdef PROFILE
        maxhit = 0;
        isyes  = 0;
#endif

#ifdef TEST_GATHER
        printf("%5d vs %5d%c : %5d x %5d\n",ar+afirst,br+bfirst,bc?'c':'n',alen,blen);
        fflush(stdout);
#endif
        setaln = 1;
        doA = doB = 0;
        amark2 = 0;
        novla  = novlb = 0;
        tbuf->top = 0;
        for (sidx = nidx; hitd[nidx].p1 == cpair; nidx = h2)
          { amark  = amark2 + PANEL_SIZE;
            amark2 = amark  - PANEL_OVERLAP;

            h2 = lidx = nidx;
            do
              { apos  = hits[nidx].apos;
                npair = hitd[++nidx].p1;
                if (apos <= amark2)
                  h2 = nidx;
              }
            while (npair == cpair && apos <= amark);

            if (nidx-lidx < minhit) continue;

            for (f = lidx; f < nidx; f++)
              { apos = hits[f].apos;
                diag = hits[f].diag >> Binshift;
                if (apos - lastp[diag] >= Kmer)
                  score[diag] += Kmer;
                else
                  score[diag] += apos - lastp[diag];
                lastp[diag] = apos;
              }

#ifdef TEST_GATHER
            printf("  %6lld upto %6d",nidx-lidx,amark);
            fflush(stdout);
#endif

            for (f = lidx; f < nidx; f++)
              { apos = hits[f].apos;
                diag = hits[f].diag;
                bpos = apos - diag;
                diag = diag >> Binshift;
#ifdef PROFILE
                if (score[diag] + scorp[diag] > maxhit)
                  maxhit = score[diag] + scorp[diag];
                if (score[diag] + scorm[diag] > maxhit)
                  maxhit = score[diag] + scorm[diag];
#endif
                if (apos > lasta[diag] &&
                     (score[diag] + scorp[diag] >= Hitmin || score[diag] + scorm[diag] >= Hitmin))
                  { if (setaln)
                      { setaln = 0;
                        align->aseq = aseq + aread[ar].boff;
                        align->bseq = bseq + bread[br].boff;
                        if (bc)
                          { CopyAndComp(bcomp,align->bseq,blen);
                            align->bseq = bcomp;
                          }
                        align->alen = alen;
                        align->blen = blen;
                        align->flags = ovla->flags = ovlb->flags = bc;
                        ovlb->bread = ovla->aread = ar + afirst;
                        ovlb->aread = ovla->bread = br + bfirst;
                        doA = (alen >= HGAP_MIN);
                        doB = (SYMMETRIC && blen >= HGAP_MIN && ! (ar == br && MG_self));
                      }
#ifdef TEST_GATHER
                    else
                      printf("\n                    ");

                    if (scorm[diag] > scorp[diag])
                      printf("  %5d.. x %5d.. %5d (%3d)",
                             bpos,apos,apos-bpos,score[diag]+scorm[diag]);
                    else
                      printf("  %5d.. x %5d.. %5d (%3d)",
                             bpos,apos,apos-bpos,score[diag]+scorp[diag]);
                    fflush(stdout);
#endif
                    nfilt += 1;

#ifdef DO_ALIGNMENT
                    bpath = Local_Alignment(align,work,MR_spec,apos-bpos,apos-bpos,apos+bpos,-1,-1);

                    { int low, hgh, ae;

                      Diagonal_Span(apath,&low,&hgh);
                      if (diag < low)
                        low = diag;
                      else if (diag > hgh)
                        hgh = diag;
                      ae = apath->aepos;
                      for (diag = low; diag <= hgh; diag++)
                        if (ae > lasta[diag])
                          lasta[diag] = ae;
#ifdef TEST_GATHER
                      printf(" %d - %d @ %d",low,hgh,apath->aepos);
                      fflush(stdout);
#endif
                    }

                    if ((apath->aepos-apath->abpos) + (apath->bepos-apath->bbpos) >= MINOVER)
                      { if (doA)
                          { if (novla >= AOmax)
                              { AOmax = 1.2*novla + MATCH_CHUNK;
                                amatch = Realloc(amatch,sizeof(Path)*AOmax,
                                                 "Reallocating match vector");
                                if (amatch == NULL)
                                  Clean_Exit(1);
                              }
                            if (tbuf->top + apath->tlen > tbuf->max)
                              { tbuf->max = 1.2*(tbuf->top+apath->tlen) + TRACE_CHUNK;
                                tbuf->trace = Realloc(tbuf->trace,sizeof(short)*tbuf->max,
                                                      "Reallocating trace vector");
                                if (tbuf->trace == NULL)
                                  Clean_Exit(1);
                              }
                            amatch[novla] = *apath;
                            amatch[novla].trace = (void *) (tbuf->top);
                            memmove(tbuf->trace+tbuf->top,apath->trace,sizeof(short)*apath->tlen);
                            novla += 1;
                            tbuf->top += apath->tlen;
                          }
                        if (doB)
                          { if (novlb >= BOmax)
                              { BOmax = 1.2*novlb + MATCH_CHUNK;
                                bmatch = Realloc(bmatch,sizeof(Path)*BOmax,
                                                        "Reallocating match vector");
                                if (bmatch == NULL)
                                  Clean_Exit(1);
                              }
                            if (tbuf->top + bpath->tlen > tbuf->max)
                              { tbuf->max = 1.2*(tbuf->top+bpath->tlen) + TRACE_CHUNK;
                                tbuf->trace = Realloc(tbuf->trace,sizeof(short)*tbuf->max,
                                                      "Reallocating trace vector");
                                if (tbuf->trace == NULL)
                                  Clean_Exit(1);
                              }
                            bmatch[novlb] = *bpath;
                            bmatch[novlb].trace = (void *) (tbuf->top);
                            memmove(tbuf->trace+tbuf->top,bpath->trace,sizeof(short)*bpath->tlen);
                            novlb += 1;
                            tbuf->top += bpath->tlen;
                          }

#ifdef TEST_GATHER
                        printf("  [%5d,%5d] x [%5d,%5d] = %4d",
                               apath->abpos,apath->aepos,apath->bbpos,apath->bepos,apath->diffs);
                        fflush(stdout);
#endif
#ifdef SHOW_OVERLAP
                        printf("\n\n                    %d(%d) vs %d(%d)\n\n",
                               ovla->aread,ovla->alen,ovla->bread,ovla->blen);
                        Print_ACartoon(stdout,align,ALIGN_INDENT);
#ifdef SHOW_ALIGNMENT
                        Compute_Trace_ALL(align,work);
                        printf("\n                      Diff = %d\n",align->path->diffs);
                        Print_Alignment(stdout,align,work,
                                        ALIGN_INDENT,ALIGN_WIDTH,ALIGN_BORDER,0,5);
#endif
#endif // SHOW_OVERLAP

                      }
#ifdef TEST_GATHER
                    else
                      printf("  No alignment %d",
                              ((apath->aepos-apath->abpos) + (apath->bepos-apath->bbpos))/2);
                    fflush(stdout);
#endif
#endif // DO_ALIGNMENT
                  }
              }

            for (f = lidx; f < nidx; f++)
              { diag = hits[f].diag >> Binshift;
                score[diag] = lastp[diag] = 0;
              }
#ifdef TEST_GATHER
            printf("\n");
            fflush(stdout);
#endif
          }

        for (f = sidx; f < nidx; f++)
          { int d;

            diag = hits[f].diag >> Binshift;
            for (d = diag; d <= maxdiag; d++)
              if (lasta[d] == 0)
                break;
              else
                lasta[d] = 0;
            for (d = diag-1; d >= mindiag; d--)
              if (lasta[d] == 0)
                break;
              else
                lasta[d] = 0;
          }

         
         { int i;

#ifdef TEST_CONTAIN
           if (novla > 1 || novlb > 1)
             printf("\n%5d vs %5d:\n",ar,br);
#endif

           if (novla > 1)
             { if (novlb > 1)
                 novla = novlb = Handle_Redundancies(amatch,novla,bmatch,align,work,tbuf);
               else
                 novla = Handle_Redundancies(amatch,novla,NULL,align,work,tbuf);
             }
           else if (novlb > 1)
             novlb = Handle_Redundancies(bmatch,novlb,NULL,align,work,tbuf);

           for (i = 0; i < novla; i++)
             { ovla->path = amatch[i];
               ovla->path.trace = tbuf->trace + (uint64) (ovla->path.trace);
               if (small)
                 Compress_TraceTo8(ovla,1);
               if (Write_Overlap(ofile1,ovla,tbytes))
                 { fprintf(stderr,"%s: Cannot write to %s too small?\n",SORT_PATH,Prog_Name);
                   Clean_Exit(1);
                 }
             }
           for (i = 0; i < novlb; i++)
             { ovlb->path = bmatch[i];
               ovlb->path.trace = tbuf->trace + (uint64) (ovlb->path.trace);
               if (small)
                 Compress_TraceTo8(ovlb,1);
               if (Write_Overlap(ofile2,ovlb,tbytes))
                 { fprintf(stderr,"%s: Cannot write to %s, too small?\n",SORT_PATH,Prog_Name);
                   Clean_Exit(1);
                 }
             }
          if (doA)
            nlas += novla;
          else
            nlas += novlb;
           ahits += novla;
           bhits += novlb;
#ifdef PROFILE
           isyes = (novla+novlb > 0);
#endif
         }

#ifdef PROFILE
        if (maxhit > MAXHIT)
          maxhit = MAXHIT;
        if (isyes)
          profyes[maxhit] += 1;
        else
          profno[maxhit] += 1;
#endif
      }

  free(tbuf->trace);
  free(bmatch);
  free(amatch);
  free(bcomp-1);

  data->nfilt = nfilt;
  data->nlas  = nlas;

  if (MR_two)
    { rewind(ofile2);
      fwrite(&bhits,sizeof(int64),1,ofile2);
      fclose(ofile2);
    }
  else
    ahits += bhits;

  rewind(ofile1);
  fwrite(&ahits,sizeof(int64),1,ofile1);
  fclose(ofile1);

  return (NULL);
}


/*******************************************************************************************
 *
 *  THE ALGORITHM
 *
 ********************************************************************************************/

static char *NameBuffer(char *aname, char *bname)
{ static char *cat = NULL;
  static int   max = -1;
  int len;

  len = strlen(aname) + strlen(bname) + 100;
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      if ((cat = (char *) realloc(cat,max+1)) == NULL)
        { fprintf(stderr,"%s: Out of memory (Making path name)\n",Prog_Name);
          Clean_Exit(1);
        }
    }
  return (cat);
}

void Match_Filter(char *aname, DAZZ_DB *ablock, char *bname, DAZZ_DB *bblock,
                  void *vasort, int alen, void *vbsort, int blen, Align_Spec *aspec)
{ THREAD     threads[NTHREADS];
  Merge_Arg  parmm[NTHREADS];
  Report_Arg parmr[NTHREADS];
  char      *fname;

  SeedPair *khit, *hhit;
  SeedPair *work1, *work2;
  int64     nhits;
  int64     nfilt, nlas;

  KmerPos  *asort, *bsort;
  int64     atot, btot;

  asort = (KmerPos *) vasort;
  bsort = (KmerPos *) vbsort;

  atot = ablock->totlen;
  btot = bblock->totlen;

  MR_tspace = Trace_Spacing(aspec);

  nfilt = nlas = nhits = 0;

  if (VERBOSE)
    printf("\nComparing %s to %s\n",aname,bname);

  if (alen == 0 || blen == 0)
    goto zerowork;

  { int    i, j, p;
    uint64 c;
    int    limit;

    MG_alist  = asort;
    MG_blist  = bsort;
    MG_ablock = ablock;
    MG_bblock = bblock;
    MG_self   = (ablock == bblock);

    parmm[0].abeg = parmm[0].bbeg = 0;
    for (i = 1; i < NTHREADS; i++)
      { p = (int) ((((int64) alen) * i) / NTHREADS);
        if (p > 0)
          { c = asort[p-1].code;
            while (asort[p].code == c)
              p += 1;
          }
        parmm[i].abeg = parmm[i-1].aend = p;
        parmm[i].bbeg = parmm[i-1].bend = find_tuple(asort[p].code,bsort,blen);
      }
    parmm[NTHREADS-1].aend = alen;
    parmm[NTHREADS-1].bend = blen;

    for (i = 0; i < NTHREADS; i++)
      for (j = 0; j < MAXGRAM; j++)
        parmm[i].hitgram[j] = 0;

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,count_thread,parmm+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

    if (VERBOSE)
      printf("\n");
    if (MEM_LIMIT > 0)
      { int64 histo[MAXGRAM];
        int64 tom, avail;

        for (j = 0; j < MAXGRAM; j++)
          histo[j] = parmm[0].hitgram[j];
        for (i = 1; i < NTHREADS; i++)
          for (j = 0; j < MAXGRAM; j++)
            histo[j] += parmm[i].hitgram[j];

        avail = (int64) (MEM_LIMIT - (sizeof_DB(ablock) + sizeof_DB(bblock))) / sizeof(KmerPos);
        if (asort == bsort || avail > alen + 2*blen)
          avail = (avail - alen) / 2;
        else
          avail = avail - (alen + blen);
        avail *= (.98 * sizeof(KmerPos)) / sizeof(SeedPair);

        tom = 0;
        for (j = 0; j < MAXGRAM; j++)
          { tom += j*histo[j];
            if (tom > avail)
              break;
          }
        limit = j;

        if (limit <= 1)
          { fprintf(stderr,"\nError: Insufficient ");
            if (MEM_LIMIT == MEM_PHYSICAL)
              fprintf(stderr," physical memory (%.1fGb), reduce block size\n",
                             (1.*MEM_LIMIT)/0x40000000ll);
            else
              { fprintf(stderr," memory allocation (%.1fGb),",(1.*MEM_LIMIT)/0x40000000ll);
                fprintf(stderr," reduce block size or increase allocation\n");
              }
            fflush(stderr);
            Clean_Exit(1);
          }
        if (limit < 10)
          { fprintf(stderr,"\nWarning: Sensitivity hampered by low ");
            if (MEM_LIMIT == MEM_PHYSICAL)
              fprintf(stderr," physical memory (%.1fGb), reduce block size\n",
                             (1.*MEM_LIMIT)/0x40000000ll);
            else
              { fprintf(stderr," memory allocation (%.1fGb),",(1.*MEM_LIMIT)/0x40000000ll);
                fprintf(stderr," reduce block size or increase allocation\n");
              }
            fflush(stderr);
          }
        if (VERBOSE)
          { printf("   Capping mutual k-mer matches over %d (effectively -t%d)\n",
                   limit,(int) sqrt(1.*limit));
            fflush(stdout);
          }

        for (i = 0; i < NTHREADS; i++)
          { parmm[i].nhits = 0;
            for (j = 1; j < limit; j++)
              parmm[i].nhits += j * parmm[i].hitgram[j];
            parmm[i].limit = limit;
          }
      }
    else
      for (i = 0; i < NTHREADS; i++)
        parmm[i].limit = INT32_MAX;

    nhits = parmm[0].nhits;
    for (i = 1; i < NTHREADS; i++)
      parmm[i].nhits = nhits += parmm[i].nhits;

    if (VERBOSE)
      { printf("   Hit count = ");
        Print_Number(nhits,0,stdout);
	if (asort == bsort || nhits*sizeof(SeedPair) >= blen*sizeof(KmerPos))
          printf("\n   Highwater of %.2fGb space\n",
                 (1. * (alen*sizeof(KmerPos) + 2*nhits*sizeof(SeedPair)) / 0x40000000ll));
        else
          printf("\n   Highwater of %.2fGb space\n",
                 (1. * ((alen + blen)*sizeof(KmerPos) + nhits*sizeof(SeedPair)) / 0x40000000ll));
        fflush(stdout);
      }

    if (nhits == 0)
      goto zerowork;

    if (asort == bsort)
      hhit = work1 = (SeedPair *) Malloc(sizeof(SeedPair)*(nhits+1),
                                         "Allocating daligner hit vectors");
    else
      { if (nhits*sizeof(SeedPair) >= blen*sizeof(KmerPos))
          bsort = (KmerPos *) Realloc(bsort,sizeof(SeedPair)*(nhits+1),
                                       "Reallocating daligner sort vectors");
        hhit = work1 = (SeedPair *) bsort;
      }
    khit = work2 = (SeedPair *) Malloc(sizeof(SeedPair)*(nhits+1),
                                        "Allocating daligner hit vectors");
    if (hhit == NULL || khit == NULL || bsort == NULL)
      Clean_Exit(1);

    MG_blist = bsort;
    MG_hits  = khit;

    for (i = NTHREADS-1; i > 0; i--)
      parmm[i].nhits = parmm[i-1].nhits;
    parmm[0].nhits = 0;

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,merge_thread,parmm+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

#ifdef TEST_PAIRS
    printf("\nSETUP SORT:\n");
    for (i = 0; i < HOW_MANY && i < nhits; i++)
      printf(" %6d / %6d / %5d / %5d\n",khit[i].aread,khit[i].bread,khit[i].apos,khit[i].diag);
#endif
  }

  { int i, j;
    int pairsort[13];
    int areads = ablock->nreads-1;
    int breads = bblock->nreads-1;
    int maxlen = ablock->maxlen;
    int abits, bbits, pbits;

    abits = 1;
    while (areads > 0)
      { areads >>= 1;
        abits   += 1;
      }

    bbits = 0;
    while (breads > 0)
      { breads >>= 1;
        bbits   += 1;
      }

    pbits = 1;
    while (maxlen > 0)
      { maxlen >>= 1;
        pbits   += 1;
      }

#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
    for (i = 0; i <= (pbits-1)/8; i++)
      pairsort[i] = 8+i;
    j = i;
    for (i = 0; i <= (bbits-1)/8; i++)
      pairsort[j+i] = 4+i;
    j += i;
    for (i = 0; i <= (abits-1)/8; i++)
      pairsort[j+i] = i;
#else
    for (i = 0; i <= (pbits-1)/8; i++)
      pairsort[i] = 11+i;
    j = i;
    for (i = 0; i <= (bbits-1)/8; i++)
      pairsort[j+i] = 7-i;
    j += i;
    for (i = 0; i <= (abits-1)/8; i++)
      pairsort[j+i] = 3-i;
#endif
    pairsort[j+i] = -1;

    khit = (SeedPair *) LSD_Sort(nhits,khit,hhit,16,16,pairsort);

    khit[nhits].aread = 0x7fffffff;
    khit[nhits].bread = 0x7fffffff;
    khit[nhits].apos  = 0x7fffffff;
    khit[nhits].diag  = 0x7fffffff;
  }

#ifdef TEST_CSORT
  { int   i;

    printf("\nCROSS SORT %lld:\n",nhits);
    for (i = 0; i < HOW_MANY && i <= nhits; i++)
      printf(" %6d / %6d / %5d / %5d\n",khit[i].aread,khit[i].bread,khit[i].apos,khit[i].diag);
  }
#endif

  { int  max_diag  = ((ablock->maxlen >> Binshift) - ((-bblock->maxlen) >> Binshift)) + 1;
    int *space;
    int  i;

    MR_ablock = ablock;
    MR_bblock = bblock;
    MR_hits   = khit;
    MR_two    = ! MG_self && SYMMETRIC;
    MR_spec   = aspec;

    { int      p, r;

      parmr[0].beg = 0;
      for (i = 1; i < NTHREADS; i++)
        { p = (nhits * i) / NTHREADS;
          if (p > 0)
            { r = khit[p-1].bread;
              while (khit[p].bread <= r)
                p += 1;
            }
          parmr[i].beg = parmr[i-1].end = p;
        }
      parmr[NTHREADS-1].end = nhits;
    }

    space = (int *) Malloc(NTHREADS*3*max_diag*sizeof(int),"Allocating space for report thread");
    if (space == NULL)
      Clean_Exit(1);

    fname = NameBuffer(aname,bname);

    for (i = 0; i < 3*max_diag*NTHREADS; i++)
      space[i] = 0;
    for (i = 0; i < NTHREADS; i++)
      { if (i == 0)
          parmr[i].score = space - ((-bblock->maxlen) >> Binshift);
        else
          parmr[i].score = parmr[i-1].lasta + max_diag;
        parmr[i].lastp = parmr[i].score + max_diag;
        parmr[i].lasta = parmr[i].lastp + max_diag;
        parmr[i].work  = New_Work_Data();

        sprintf(fname,"%s/%s.%s.N%d.las",SORT_PATH,aname,bname,i+1);
        parmr[i].ofile1 = Fopen(fname,"w");
        if (parmr[i].ofile1 == NULL)
          Clean_Exit(1);

        if (MG_self)
          parmr[i].ofile2 = parmr[i].ofile1;
        else if (SYMMETRIC)
          { sprintf(fname,"%s/%s.%s.N%d.las",SORT_PATH,bname,aname,i+1);
            parmr[i].ofile2 = Fopen(fname,"w");
            if (parmr[i].ofile2 == NULL)
              Clean_Exit(1);
          }
      }

#ifdef NOTHREAD

    for (i = 0; i < NTHREADS; i++)
      report_thread(parmr+i);

#else

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,report_thread,parmr+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

#endif

#ifdef PROFILE
    printf("\n Hit Score distribution:\n");
    for (i = MAXHIT; i >= 0; i--)
      { int   j;
        int64 nyes, nno;

        nyes = 0;
        nno  = 0;
        for (j = 0; j < NTHREADS; j++)
          { nyes += parmr[j].profyes[i];
            nno  += parmr[j].profno[i];
          }
        if (nyes+nno > 0)
          printf(" %4d: %6lld %6lld\n",i,nyes,nno);
      }
#endif

    if (VERBOSE)
      for (i = 0; i < NTHREADS; i++)
        { nfilt += parmr[i].nfilt;
          nlas  += parmr[i].nlas;
        }

    for (i = 0; i < NTHREADS; i++)
      Free_Work_Data(parmr[i].work);
    free(space);
  }

  free(work2);
  free(work1);
  goto epilogue;

zerowork:
  { FILE *ofile;
    int   i;

    fname = NameBuffer(aname,bname);

    nhits  = 0;
    for (i = 0; i < NTHREADS; i++)
      { sprintf(fname,"%s/%s.%s.N%d.las",SORT_PATH,aname,bname,i+1);
        ofile = Fopen(fname,"w");
        fwrite(&nhits,sizeof(int64),1,ofile);
        fwrite(&MR_tspace,sizeof(int),1,ofile);
        fclose(ofile);
        if (! MG_self && SYMMETRIC)
          { sprintf(fname,"%s/%s.%s.N%d.las",SORT_PATH,bname,aname,i+1);
            ofile = Fopen(fname,"w");
            fwrite(&nhits,sizeof(int64),1,ofile);
            fwrite(&MR_tspace,sizeof(int),1,ofile);
            fclose(ofile);
          }
      }
  }

epilogue:

  if (VERBOSE)
    { int width;

      if (nhits <= 0)
        width = 1;
      else
        width = ((int) log10((double) nhits)) + 1;
      width += (width-1)/3;

      printf("\n     ");
      Print_Number(nhits,width,stdout);
      printf(" %d-mers (%e of matrix)\n     ",Kmer,(1.*nhits/atot)/btot);
      Print_Number(nfilt,width,stdout);
      printf(" seed hits (%e of matrix)\n     ",(1.*nfilt/atot)/btot);
      Print_Number(nlas,width,stdout);
      printf(" confirmed hits (%e of matrix)\n",(1.*nlas/atot)/btot);
      fflush(stdout);
    }
}
