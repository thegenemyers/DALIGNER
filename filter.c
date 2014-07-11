/*******************************************************************************************
 *
 *  Fast local alignment filter for long, noisy reads based on "dumbing down" of my RECOMB 2005
 *     filter with Jens Stoye, and a "smarting up" of the k-mer matching by turning it into
 *     a threaded sort and merge paradigm using a super cache coherent radix sort.
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
#include "filter.h"
#include "align.h"

#define THREAD    pthread_t

#define MAX_BIAS  2  //  In -b mode, don't consider tuples with specificity
                     //     <= 4 ^ -(kmer-MAX_BIAS)

#undef  TEST_LSORT
#undef  TEST_KSORT
#undef  TEST_PAIRS
#undef  TEST_CSORT
#undef  TEST_GATHER
#define    HOW_MANY   30000   //  Print first HOW_MANY items for each of the TEST options above

#define TEST_HIT              //  When off, just tuple filtering alone, no check
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

#define CODE(x) ((x) & 0xffffffffll)
#define READ(x) ((x) >> 48)
#define RPOS(x) (((x) >> 32) & 0xffffll)

static  uint64 *Asort = NULL;
static  uint64 *Csort = NULL;

static  int     Alen, Clen;

/*******************************************************************************************
 *
 *  PARAMETER SETUP
 *
 ********************************************************************************************/

static int Kmer;
static int Hitmin;
static int Binshift;
static int Suppress;

static int Kshift;         //  2*Kmer
static int Kmask;          //  4^Kmer-1
static int Kpowr;          //  4^Kmer
static int Binwidth;       //  2^Binshift
static int TooFrequent;    //  (Suppress != 0) ? Suppress : INT32_MAX

int Set_Filter_Params(int kmer, int binshift, int suppress, int hitmin)
{ if (kmer <= 1)
    return (1);

  Kmer     = kmer;
  Binshift = binshift;
  Suppress = suppress;
  Hitmin   = hitmin;

  Kshift = 2*Kmer;
  Kpowr  = (1 << Kshift);
  Kmask  = Kpowr-1;

  Binwidth = (1 << Binshift);

  if (Suppress == 0)
    TooFrequent = INT32_MAX;
  else
    TooFrequent = Suppress;

  return (0);
}


/*******************************************************************************************
 *
 *  LEXICOGRAPHIC SORT
 *
 ********************************************************************************************/

#define BMER      4
#define BSHIFT    8             //  = 2*BMER
#define BPOWR   256             //  = 2^BSHIFT
#define BMASK  0xff             //  = BPOWR-1

static int QSHIFT;              //  = BSHIFT - NSHIFT
static int QPOWR;               //  = 2^(BSHIFT+NSHIFT)
static int QMASK;               //  = QPOWR-1

static int     LEX_shift;
static int     LEX_zsize;
static int     LEX_first;
static int     LEX_last;
static uint64 *LEX_src;
static uint64 *LEX_trg;

typedef struct
  { int  beg;
    int  end;
    int  tptr[BPOWR];
    int  sptr[NTHREADS*BPOWR];
  } Lex_Arg;

static void *lex_thread(void *arg)
{ Lex_Arg    *data  = (Lex_Arg *) arg;
  int        *sptr  = data->sptr;
  int        *tptr  = data->tptr;
  int         shift = LEX_shift;
  int         zsize = LEX_zsize;
  uint64     *src   = LEX_src;
  uint64     *trg   = LEX_trg;
  int         i, n, x;
  uint64      c, b;

  n = data->end;
  if (LEX_first)
    for (i = data->beg; i < n; i++)
      { c = src[i];
        x = tptr[c&BMASK]++;
        trg[x] = c;
        sptr[((c >> QSHIFT) & QMASK) + x/zsize] += 1;
      }
  else if (LEX_last)
    for (i = data->beg; i < n; i++)
      { c = src[i];
        trg[tptr[(c >> shift) & BMASK]++] = c;
      }
  else
    for (i = data->beg; i < n; i++)
      { c = src[i];
        b = (c >> shift);
        x = tptr[b&BMASK]++;
        trg[x] = c;
        sptr[((b >> QSHIFT) & QMASK) + x/zsize] += 1;
      }

  return (NULL);
}

static uint64 *lex_sort(int nbits, int sbit, uint64 *src, uint64 *trg, Lex_Arg *parmx)
{ THREAD  threads[NTHREADS];

  int     len;
  uint64 *xch;
  int     i, j, x, z;

  QSHIFT = BSHIFT - NSHIFT;
  QPOWR  = (BPOWR << NSHIFT);
  QMASK  = (BMASK << NSHIFT);

  len       = parmx[NTHREADS-1].end;
  LEX_zsize = (len-1)/NTHREADS + 1;
  LEX_src   = src;
  LEX_trg   = trg;

  for (LEX_shift = sbit; LEX_shift < nbits; LEX_shift += BSHIFT)
    { LEX_last  = (LEX_shift + BSHIFT >= nbits);
      LEX_first = (LEX_shift == 0);

      if (LEX_shift == sbit)
        { for (i = 0; i < NTHREADS; i++)
            for (z = 0; z < NTHREADS*BPOWR; z++)
              parmx[i].sptr[z] = 0;
        }
      else
        { x = 0;
          for (i = 0; i < NTHREADS; i++)
            { parmx[i].beg = x;
              parmx[i].end = x = LEX_zsize*(i+1);
              for (j = 0; j < BPOWR; j++)
                parmx[i].tptr[j] = 0;
            }
          parmx[NTHREADS-1].end = len;

          for (j = 0; j < BPOWR; j++)
            { x = (j << NSHIFT);
              for (z = 0; z < NTHREADS; z++)
                for (i = 0; i < NTHREADS; i++)
                  { parmx[i].tptr[j] += parmx[z].sptr[x+i];
                    parmx[z].sptr[x+i] = 0;
                  }
            }
        }

      x = 0;
      for (j = 0; j < BPOWR; j++)
        for (i = 0; i < NTHREADS; i++)
          { z = parmx[i].tptr[j];
            parmx[i].tptr[j] = x;
            x += z;
          }

      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,lex_thread,parmx+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);

      xch     = LEX_src;
      LEX_src = LEX_trg;
      LEX_trg = xch;

#ifdef TEST_LSORT
      printf("\nLSORT %d\n",LEX_shift);
      x = (1 << (LEX_shift+BSHIFT))-1;
      for (i = 0; i < len; i++)
        { printf("%5d: %5llx %5llx %5llx %5llx : %4llx",
                 i,LEX_src[i]>>48,(LEX_src[i]>>32)&0xffffll,(LEX_src[i]>>16)&0xffffll,
                 LEX_src[i]&0xffffll,LEX_src[i]&x);
          if (i > 0 && (LEX_src[i] & x) < (LEX_src[i-1] & x))
            printf(" OO");
          printf("\n");
        }
#endif
    }

  return (LEX_src);
}


/*******************************************************************************************
 *
 *  INDEX BUILD
 *
 ********************************************************************************************/

static int *NormShift = NULL;
static int  LogNorm, LogThresh;
static int  LogBase[4];

static HITS_DB    *TA_block;
static uint64     *TA_list;
static HITS_TRACK *TA_dust;

typedef struct
  { int  tnum;
    int *kptr;
    int  fill;
  } Tuple_Arg;


static void *tuple_thread(void *arg)
{ Tuple_Arg  *data  = (Tuple_Arg *) arg;
  int         tnum  = data->tnum;
  int        *kptr  = data->kptr;
  uint64     *list  = TA_list;
  int         n, i, m, c, x, p;
  uint64      h;
  char       *s;

  c  = TA_block->nreads;
  i  = (c * tnum) >> NSHIFT;
  n  = TA_block->reads[i].boff;
  s  = ((char *) (TA_block->bases)) + n;
  n -= Kmer*i;

  if (TA_dust != NULL)
    { HITS_READ *reads = TA_block->reads;
      int       *anno1 = ((int *) (TA_dust->anno)) + 1;
      int       *point = (int *) (TA_dust->data);
      int        a, b, f, q;

      f = anno1[i-1];
      for (m = (c * (tnum+1)) >> NSHIFT; i < m; i++)
        { b = f;
          h = ((uint64) i) << 48;
          f = anno1[i];
          q = reads[i].end-reads[i].beg;
          for (a = b; a <= f; a += 2)
            { if (a == b)
                p = 0;
              else
                p = point[a-1]+1;
              if (a == f)
                q = reads[i].end-reads[i].beg;
              else
                q = point[a];
              if (p+Kmer <= q)
                { c = 0;
                  for (x = 1; x < Kmer; x++)
                    c = (c << 2) | s[p++];
                  while (p < q)
                    { x = s[p];
                      c = ((c << 2) | x) & Kmask;
                      list[n++] = h | (((uint64) (p++)) << 32) | c;
                      kptr[c & BMASK] += 1;
                    }
                }
            }
          s += (q+1);
        }

      m = TA_block->reads[m].boff - Kmer*m;
      kptr[BMASK] += (data->fill = m-n);
      while (n < m)
        list[n++] = 0xffffffffffffffffll;
    }

  else
    for (m = (c * (tnum+1)) >> NSHIFT; i < m; i++)
      { h = ((uint64) i) << 48;
        c = p = 0;
        for (x = 1; x < Kmer; x++)
          c = (c << 2) | s[p++];
        while ((x = s[p]) != 4)
          { c = ((c << 2) | x) & Kmask;
            list[n++] = h | (((uint64) (p++)) << 32) | c;
            kptr[c & BMASK] += 1;
          }
        s += (p+1);
      }

  return (NULL);
}

static void *biased_tuple_thread(void *arg)
{ Tuple_Arg  *data  = (Tuple_Arg *) arg;
  int         tnum  = data->tnum;
  int        *kptr  = data->kptr;
  uint64     *list  = TA_list;
  int         n, i, m;
  int         c, x, a, k, p, d;
  uint64      h;
  char       *s, *t;

  c  = TA_block->nreads;
  i  = (c * tnum) >> NSHIFT;
  n  = TA_block->reads[i].boff;
  s  = ((char *) (TA_block->bases)) + n;
  n -= Kmer*i;

  if (TA_dust != NULL)
    { HITS_READ *reads = TA_block->reads;
      int       *anno1 = ((int *) (TA_dust->anno)) + 1;
      int       *point = (int *) (TA_dust->data);
      int        j, b, f, q;

      f = anno1[i-1];
      for (m = (c * (tnum+1)) >> NSHIFT; i < m; i++)
        { b = f;
          f = anno1[i];
          h = ((uint64) i) << 48;
          t = s+1;
          q = reads[i].end-reads[i].beg;
          for (j = b; j <= f; j += 2)
            { if (j == b)
                p = 0;
              else
                p = point[j-1]+1;
              if (j == f)
                q = reads[i].end-reads[i].beg;
              else
                q = point[j];
              if (p+Kmer <= q)
                { c = a = 0;
                  k = 1;
                  while (p < q)
                    { x = s[p];
                      a += LogBase[x];
                      c  = ((c << 2) | x);
                      while (a < LogNorm && k < Kmer)
                        { if (++p >= q)
                            break;
                          k += 1;
                          x  = s[p];
                          a += LogBase[x];
                          c  = ((c << 2) | x);
                        }
                      while (1)
	                { int u = a-LogBase[(int) t[p-k]];
                          if (u < LogNorm) break;
                          a  = u;
                          k -= 1;
                        }
                      if (a > LogThresh)
                        { d = ((c << NormShift[k]) & Kmask);
                          list[n++] = h | (((uint64) p) << 32) | d;
                          kptr[d & BMASK] += 1;
                        }
                      p += 1;
                      a -= LogBase[(int) s[p-k]];
                    }
                }
            }
          s += (q+1);
	}
    }

  else
    for (m = (c * (tnum+1)) >> NSHIFT; i < m; i++)
      { t = s+1;
        c = p = a = 0;
        k = 1;
        h = ((uint64) i) << 48;
        while ((x = s[p]) != 4)
          { a += LogBase[x];
            c  = ((c << 2) | x);
            while (a < LogNorm && k < Kmer)
              { if ((x = s[++p]) == 4)
                  goto eoread2;
                k += 1;
                a += LogBase[x];
                c  = ((c << 2) | x);
              }
            while (1)
	      { int u = a-LogBase[(int) t[p-k]];
                if (u < LogNorm) break;
                a  = u;
                k -= 1;
              }
            if (a > LogThresh)
              { d = ((c << NormShift[k]) & Kmask);
                list[n++] = h | (((uint64) p) << 32) | d;
                kptr[d & BMASK] += 1;
              }
            p += 1;
            a -= LogBase[(int) s[p-k]];
          }
      eoread2:
        s += (p+1);
      }

  m = TA_block->reads[m].boff - Kmer*m;
  kptr[BMASK] += (data->fill = m-n);
  while (n < m)
    list[n++] = 0xffffffffffffffffll;

  return (NULL);
}

static uint64 *FR_src;
static uint64 *FR_trg;

typedef struct
  { int  beg;
    int  end;
    int  kept;
  } Comp_Arg;

static void *compsize_thread(void *arg)
{ Comp_Arg   *data  = (Comp_Arg *) arg;
  int         end   = data->end;
  uint64     *src   = FR_src;
  int         n, i, c, p;
  uint64      h, g;

  i = data->beg;
  h = CODE(src[i]);
  n = 0;
  while (i < end)
    { p = i++;
      while ((g = CODE(src[i])) == h)
        i += 1;
      if ((c = (i-p)) < TooFrequent)
        n += c;
      h = g;
    }

  data->kept = n;
  return (NULL);
}

static void *compress_thread(void *arg)
{ Comp_Arg   *data  = (Comp_Arg *) arg;
  int         end   = data->end;
  uint64     *src   = FR_src;
  uint64     *trg   = FR_trg;
  int         n, i, p;
  uint64      h, g;

  i = data->beg;
  h = CODE(src[i]);
  n = data->kept;
  while (i < end)
    { p = i++;
      while ((g = CODE(src[i])) == h)
        i += 1;
      if (i-p < TooFrequent)
        { while (p < i)
            trg[n++] = src[p++];
        }
      h = g;
    }

  return (NULL);
}

uint64 *sort_kmers(HITS_DB *block, int *len)
{ THREAD    threads[NTHREADS];
  Tuple_Arg parmt[NTHREADS];
  Comp_Arg  parmf[NTHREADS];
  Lex_Arg   parmx[NTHREADS];

  uint64   *src, *trg, *rez;
  int       kmers;
  int       i, j, x, z;
  uint64    h;

  if (NormShift == NULL && BIASED)
    { NormShift = (int *) Malloc(sizeof(int)*(Kmer+1),"Allocating sort_kmers bias shift");
      if (NormShift == NULL)
        exit (1);
      for (i = 0; i <= Kmer; i++)
        NormShift[i] = Kshift - 2*i;
      LogNorm = 10000 * Kmer;
      LogThresh = 10000 * (Kmer-MAX_BIAS);
      for (i = 0; i < 4; i++)
        LogBase[i] = ceil( -10000. * log(block->freq[i]) / log(4.));
    }

  kmers = block->reads[block->nreads].boff - Kmer * block->nreads;

  if (kmers <= 0)
    goto no_mers;

  if (( (Kshift-1)/BSHIFT + (TooFrequent < INT32_MAX) ) & 0x1)
    { trg = (uint64 *) Malloc(sizeof(uint64)*(kmers+1),"Allocating sort_kmers vectors");
      src = (uint64 *) Malloc(sizeof(uint64)*(kmers+1),"Allocating sort_kmers vectors");
    }
  else
    { src = (uint64 *) Malloc(sizeof(uint64)*(kmers+1),"Allocating sort_kmers vectors");
      trg = (uint64 *) Malloc(sizeof(uint64)*(kmers+1),"Allocating sort_kmers vectors");
    }
  if (src == NULL || trg == NULL)
    exit (1);

  if (VERBOSE)
    { printf("\nKmer count = ");
      Print_Number((int64) kmers,0,stdout);
      printf("\nUsing %.2fGb of space\n",(1. * kmers) / 67108864);
      fflush(stdout);
    }

  TA_block = block;
  TA_list  = src;
  
  { HITS_TRACK *t;

    TA_dust = NULL;
    for (t = block->tracks; t != NULL; t++)
      if (strcmp(t->name,"dust") == 0)
        break;
    TA_dust = t;
  }

  for (i = 0; i < NTHREADS; i++)
    { parmt[i].tnum = i;
      parmt[i].kptr = parmx[i].tptr;
      for (j = 0; j < BPOWR; j++)
        parmt[i].kptr[j] = 0;
    }

  if (BIASED)
    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,biased_tuple_thread,parmt+i);
  else
    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,tuple_thread,parmt+i);

  for (i = 0; i < NTHREADS; i++)
    pthread_join(threads[i],NULL);

  x = z = 0;
  for (i = 0; i < NTHREADS; i++)
    { parmx[i].beg = x;
      z += block->nreads;
      j = z >> NSHIFT;
      parmx[i].end = x = block->reads[j].boff - j*Kmer;
    }

  rez = lex_sort(Kshift,0,src,trg,parmx);
  if (BIASED || TA_dust != NULL)
    for (i = 0; i < NTHREADS; i++)
      kmers -= parmt[i].fill;
  rez[kmers] = Kpowr;

  if (TooFrequent < INT32_MAX && kmers > 0)
    { parmf[0].beg = 0;
      for (i = 1; i < NTHREADS; i++)
        { x = (((int64) i)*kmers) >> NSHIFT;
          h = CODE(rez[x-1]);
          while (CODE(rez[x]) == h)
            x += 1;
          parmf[i-1].end = parmf[i].beg = x;
        }
      parmf[NTHREADS-1].end = kmers;

      if (src == rez)
        { FR_src = src;
          FR_trg = rez = trg;
        }
      else
        { FR_src = trg;
          FR_trg = rez = src;
        }

      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,compsize_thread,parmf+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);

      x = 0;
      for (i = 0; i < NTHREADS; i++)
        { z = parmf[i].kept;
          parmf[i].kept = x;
          x += z;
        }
      kmers = x;

      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,compress_thread,parmf+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);

      rez[kmers] = Kpowr;
    }
    
  if (src != rez)
    free(src);
  else
    free(trg);

#ifdef TEST_KSORT
  { int i;

    printf("\nKMER SORT:\n");
    for (i = 0; i < HOW_MANY && i < kmers; i++)
      { int64 c = rez[i];
        printf(" %5lld / %5lld / %10lld\n",c >> 48,(c >> 32) & 0xffffll,CODE(c));
      }
    fflush(stdout);
  }
#endif

  if (VERBOSE)
    { if (TooFrequent < INT32_MAX)
        { printf("Revised kmer count = ");
          Print_Number((int64) kmers,0,stdout);
          printf("\n");
        }
      printf("Index occupies %.2fGb\n",(1. * kmers) / 134217728);
      fflush(stdout);
    }

  if (kmers <= 0)
    { free(rez);
      goto no_mers;
    }

  *len = kmers;
  return (rez);

no_mers:
  *len = 0;
  return (NULL);
}

void Build_Table(HITS_DB *block)
{ if (Asort == NULL)
    Asort = sort_kmers(block,&Alen);
  else
    Csort = sort_kmers(block,&Clen);
}


/*******************************************************************************************
 *
 *  FILTER MATCH
 *
 ********************************************************************************************/

int count_hits(int b, int e, uint64 *hit)
{ int c, q, p;

  c = Kmer;
  q = (hit[b++] & 0xffffll);
  while (b < e)
    { p = (hit[b++] & 0xffffll);
      if (p < q)
        c += Kmer;
      else if (p-q < Kmer)
        c += p-q;
      else
        c += Kmer;
      q = p;
    }
  return (c);
}

int find_tuple(uint64 x, uint64 *a, int n) // smallest k s.t. a[k] >= x (or n if does not exist)
{ int l, r, m;

  l = 0;
  r = n;
  while (l < r)
    { m = ((l+r) >> 1);
      if (CODE(a[m]) < x)
        l = m+1;
      else
        r = m;
    }
  return (l);
}

static uint64 *MG_alist;
static uint64 *MG_blist;
static uint64 *MG_hits;
static int     MG_self;

typedef struct
  { int    abeg, aend;
    int    bbeg, bend;
    uint64 nhits;
    int   *kptr;
  } Merge_Arg;

static void *count_thread(void *arg)
{ Merge_Arg  *data  = (Merge_Arg *) arg;
  uint64     *asort = MG_alist;
  uint64     *bsort = MG_blist;
  int         nhits = 0;
  int         aend  = data->aend;

  int    ia, ib;
  int    jb, ja;
  uint64 ca, cb;
  uint64 da, db;
  uint64 x;

  ia = data->abeg;
  ca = CODE(asort[ia]);
  ib = data->bbeg;
  cb = CODE(bsort[ib]);
  if (MG_self)
    { while (1)
        { while (cb < ca)
            cb = CODE(bsort[++ib]);
          while (cb > ca)
            ca = CODE(asort[++ia]);
          if (cb == ca)
            { if (ia >= aend) break;

              jb = ib;
              db = cb;
              do
                { x = READ(asort[ia]);
                  while (db == cb && READ(bsort[ib]) < x)
                    db = CODE(bsort[++ib]);
                  nhits += (ib-jb);
                }
              while ((da = CODE(asort[++ia])) == ca);
              while (db == cb)
                db = CODE(bsort[++ib]);

/*
              jb = ib++;
              while ((db = CODE(bsort[ib])) == cb)
                ib += 1;
              do
                { x = READ(asort[ia]);
                  while (jb < ib && READ(bsort[jb]) <= x)
                    jb += 1;
                  nhits += (ib-jb);
                  ia += 1;
                }
              while ((da = CODE(asort[ia])) == ca);
*/

              ca = da;
              cb = db;
            }
        }
    }
  else
    { while (1)
        { while (cb < ca)
            cb = CODE(bsort[++ib]);
          while (cb > ca)
            ca = CODE(asort[++ia]);
          if (cb == ca)
            { if (ia >= aend) break;
              ja = ia++;
              while ((da = CODE(asort[ia])) == ca)
                ia += 1;
              jb = ib++;
              while ((db = CODE(bsort[ib])) == cb)
                ib += 1;
              nhits += (ia-ja)*(ib-jb);
              ca = da;
              cb = db;
            }
        }
    }

  data->nhits = nhits;

  return (NULL);
}

static void *merge_thread(void *arg)
{ Merge_Arg  *data  = (Merge_Arg *) arg;
  int        *kptr  = data->kptr;
  uint64     *asort = MG_alist;
  uint64     *bsort = MG_blist;
  uint64     *hits  = MG_hits;
  int         nhits = data->nhits;
  int         aend  = data->aend;

  int    ia, ib;
  int    jb, ja;
  uint64 ca, cb;
  uint64 da, db;
  uint64 v, x;
  int    a, b, d;

  ia = data->abeg;
  ca = CODE(asort[ia]);
  ib = data->bbeg;
  cb = CODE(bsort[ib]);
  if (MG_self)
    { while (1)
        { while (cb < ca)
            cb = CODE(bsort[++ib]);
          while (cb > ca)
            ca = CODE(asort[++ia]);
          if (cb == ca)
            { if (ia >= aend) break;

              jb = ib;
              db = cb;
              do
                { v = (asort[ia] >> 32);
                  x = (v >> 16);
                  while (db == cb && READ(bsort[ib]) < x)
                    db = CODE(bsort[++ib]);
                  if ((d = ib-jb) > 0)
                    { kptr[v & BMASK] += d;
                      v = (v << 16);
                      for (b = jb; b < ib; b++)
                        { x = bsort[b];
                          hits[nhits++] = (x & 0xffff000000000000ll) | v | ((x >> 32) & 0xffffll);
                        }
                    }
                }
              while ((da = CODE(asort[++ia])) == ca);
              while (db == cb)
                db = CODE(bsort[++ib]);

/*
              jb = ib++;
              while ((db = CODE(bsort[ib])) == cb)
                ib += 1;
              do
                { v = (asort[ia] >> 32);
                  x = (v >> 16);
                  while (jb < ib && READ(bsort[jb]) <= x)
                    jb += 1;
                  if ((d = ib-jb) > 0)
                    { kptr[v & BMASK] += d;
                      v = (v << 16);
                      for (b = jb; b < ib; b++)
                        { x = bsort[b];
                          hits[nhits++] = (x & 0xffff000000000000ll) | v | ((x >> 32) & 0xffffll);
                        }
                    }
                  ia += 1;
                }
              while ((da = CODE(asort[ia])) == ca);
*/

              ca = da;
              cb = db;
            }
        }
    }
  else
    { while (1)
        { while (cb < ca)
            cb = CODE(bsort[++ib]);
          while (cb > ca)
            ca = CODE(asort[++ia]);
          if (cb == ca)
            { if (ia >= aend) break;
              ja = ia++;
              while ((da = CODE(asort[ia])) == ca)
                ia += 1;
              jb = ib++;
              while ((db = CODE(bsort[ib])) == cb)
                ib += 1;
              d = ib-jb;
              for (a = ja; a < ia; a++)
                { v = (asort[a] >> 32);
                  kptr[v & BMASK] += d;
                  v = (v << 16);
                  for (b = jb; b < ib; b++)
                    { x = bsort[b];
                      hits[nhits++] = (x & 0xffff000000000000ll) | v | ((x >> 32) & 0xffffll);
                    }
                }
              ca = da;
              cb = db;
            }
        }
    }

  return (NULL);
}


static HITS_DB *MR_ablock;
static HITS_DB *MR_bblock;
static uint64  *MR_hits;
static int      MR_nhits;
static int      MR_comp;
static int      MR_two;

typedef struct
  { int         beg, end;
    int        *score;
    int        *lastp;
    Work_Data  *work;
    Align_Spec *aspec;
    FILE       *ofile1;
    FILE       *ofile2;
    int         nfilt;
    int         ncheck;
  } Report_Arg;

void *report_thread(void *arg)
{ Report_Arg  *data   = (Report_Arg *) arg;
  uint64      *hits   = MR_hits;
  char        *aseq   = (char *) (MR_ablock->bases);
  char        *bseq   = (char *) (MR_bblock->bases);
  HITS_READ   *aread  = MR_ablock->reads;
  HITS_READ   *bread  = MR_bblock->reads;
  int         *score  = data->score;
  int         *lastp  = data->lastp;
#ifdef TEST_HIT
  Work_Data   *work   = data->work;
#endif
  Align_Spec  *aspec  = data->aspec;
  FILE        *ofile1 = data->ofile1;
  FILE        *ofile2 = data->ofile2;
  int          afirst = MR_ablock->bfirst;
  int          bfirst = MR_bblock->bfirst;

  Overlap     _ovla, *ovla = &_ovla;
  Overlap     _ovlb, *ovlb = &_ovlb;
  Alignment   _align, *align = &_align;
#ifdef TEST_HIT
  Path        *apath;
#endif
  Path        *bpath = &(ovlb->path);
  int          nfilt = 0;
  int64        ohits = 0;
  int          tspace, small, tbytes;

  uint64 *hitc;
  int     minhit, end;

  uint64 p, q;
  int    g, h, f;

  minhit = (Hitmin-1)/Kmer + 1;
  hitc   = hits + (minhit-1);
  end    = data->end - minhit;

  //  In ovl and align roles of A and B are reversed, as the B sequence must be the
  //    complemented sequence !!

  align->flags = ovla->flags = ovlb->flags = MR_comp;
  align->path = bpath;

  tspace = Trace_Spacing(aspec);
  if (tspace <= TRACE_XOVR)
    { small  = 1;
      tbytes = sizeof(uint8);
    }
  else
    { small  = 0;
      tbytes = sizeof(uint16);
    }

  fwrite(&ohits,sizeof(int64),1,ofile1);
  fwrite(&tspace,sizeof(int),1,ofile1);
  if (MR_two)
    { fwrite(&ohits,sizeof(int64),1,ofile2);
      fwrite(&tspace,sizeof(int),1,ofile2);
    }

  h = data->beg;
  for (p = (hits[h] >> 32); h < end; p = q)
    if ((hitc[h] >> 32) != p)
      { h += 1;
        while ((q = (hits[h] >> 32)) == p)
          h += 1;
      }
    else
      { int   apos, bpos, diag;
        int   lasta, lastd;
        int   ar, br;

        g = h;
        do
          { bpos = (hits[h] & 0xffffll);
            apos = (hits[h] >> 16) & 0xffffll;
            diag = (apos - bpos) >> Binshift;
            if (apos - lastp[diag] >= Kmer)
              score[diag] += Kmer;
            else
              score[diag] += apos - lastp[diag];
            lastp[diag] = apos;
            q = (hits[++h] >> 32);
          }
        while (q == p);

#ifdef TEST_GATHER
        printf("%5lld vs %5lld : %3d",(p>>16)&0xffffll,p&0xffffll,h-g);
#endif

        ar = lasta = -1;
        lastd = -(Kmer+1);
        for (f = g; f < h; f++)
          { bpos = (hits[f] & 0xffffll);
            apos = (hits[f] >> 16) & 0xffffll;
            diag = apos - bpos;
            if ((lastd != diag && apos >= lasta) || (lastd == diag && apos > lasta+Kmer))
              { diag >>= Binshift;
                if (score[diag] + score[diag+1] >= Hitmin || score[diag] + score[diag-1] >= Hitmin)
                  { if (ar < 0)
                      { ar = (p & 0xffffll);
                        br = ((p >> 16) & 0xffffll);
                        align->bseq = aseq + aread[ar].boff;
                        align->aseq = bseq + bread[br].boff;
                        align->blen = ovlb->blen = ovla->alen = aread[ar].end-aread[ar].beg;
                        align->alen = ovlb->alen = ovla->blen = bread[br].end-bread[br].beg;
                        ovlb->bread = ovla->aread = ar + afirst;
                        ovlb->aread = ovla->bread = br + bfirst;
                      }
#ifdef TEST_GATHER
                    else
                      printf("\n                    ");
                    if (score[diag-1] > score[diag+1])
                      printf("  %5d.. x %5d.. %5d (%3d)",
                             bpos,apos,apos-bpos,score[diag]+score[diag-1]);
                    else
                      printf("  %5d.. x %5d.. %5d (%3d)",
                             bpos,apos,apos-bpos,score[diag]+score[diag+1]);
#endif
                    nfilt += 1;

#ifdef TEST_HIT

                    apath = Local_Alignment(align,work,aspec,bpos,apos);
                    lasta = bpath->bepos;
                    lastd = lasta - bpath->aepos;
                    if ((apath->aepos - apath->abpos) + (apath->bepos - apath->bbpos) >= MINOVER)
                      { ovla->path = *apath;
                        if (small)
                          { Compress_TraceTo8(ovla);
                            Compress_TraceTo8(ovlb);
                          }
                        Write_Overlap(ofile1,ovla,tbytes);
                        Write_Overlap(ofile2,ovlb,tbytes);
                        ohits += 1;

#ifdef TEST_GATHER
                        printf("  [%5d,%5d] x [%5d,%5d] = %4d",
                               bpath->abpos,bpath->aepos,bpath->bbpos,bpath->bepos,bpath->diffs);
#endif
#ifdef SHOW_OVERLAP
                        printf("\n\n                    %d(%d) vs %d(%d)\n\n",
                               ovlb->aread,ovlb->alen,ovlb->bread,ovlb->blen);
                        Print_ACartoon(stdout,align,ALIGN_INDENT);
#ifdef SHOW_ALIGNMENT
                        if (small)
                          Decompress_TraceTo16(ovlb);
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
                              ((apath->aepos - apath->abpos) + (apath->bepos - apath->bbpos))/2);
#endif

#else // TEST_HIT

                    break;

#endif // TEST_HIT
                  }
              }
          }
#ifdef TEST_GATHER
        printf("\n");
#endif

        for (; g < h; g++)
          { bpos = (hits[g] & 0xffffll);
            apos = (hits[g] >> 16) & 0xffffll;
            diag = (apos - bpos) >> Binshift;
            score[diag] = lastp[diag] = 0;
          }
      }

  data->nfilt  = nfilt;
  data->ncheck = ohits;

  if (MR_two)
    { rewind(ofile2);
      fwrite(&ohits,sizeof(int64),1,ofile2);
      fclose(ofile2);
    }
  else
    ohits *= 2;

  rewind(ofile1);
  fwrite(&ohits,sizeof(int64),1,ofile1);
  fclose(ofile1);

  return (NULL);
}


/*******************************************************************************************
 *
 *  THE ALGORITHM
 *
 ********************************************************************************************/

void Match_Filter(char *aname, HITS_DB *ablock, char *bname, HITS_DB *bblock,
                  int self, int comp, Align_Spec *aspec)
{ THREAD     threads[NTHREADS];
  Merge_Arg  parmm[NTHREADS];
  Lex_Arg    parmx[NTHREADS];
  Report_Arg parmr[NTHREADS];

  uint64   *khit, *hhit;
  uint64   *work1, *work2;
  int64     nhits;
  int64     nfilt, ncheck;

  uint64   *asort, *bsort;
  int       alen, blen;

  int64     atot,   btot;
  int       areads, breads;
  int       amax, bmax;

  areads = ablock->nreads;
  atot   = ablock->totlen;
  amax   = ablock->maxlen;

  breads = bblock->nreads;
  btot   = bblock->totlen;
  bmax   = bblock->maxlen;

  if (comp)
    { asort = Csort;
      alen  = Clen;
    }
  else
    { asort = Asort;
      alen  = Alen;
    }

  if (self)
    { bsort = Asort;
      blen  = Alen;
    }
  else
    bsort = sort_kmers(bblock,&blen);

  nfilt = ncheck = nhits = 0;

  if (alen == 0 || blen == 0)
    goto zerowork;

  { int    i;
    uint64 p, c;

    MG_alist = asort;
    MG_blist = bsort;
    MG_self  = self;

    parmm[0].abeg = parmm[0].bbeg = 0;
    parmm[0].kptr = parmx[0].tptr;
    for (i = 1; i < NTHREADS; i++)
      { p = (alen * ((int64) i)) >> NSHIFT;
        if (p > 0)
          { c = CODE(asort[p-1]);
            while (CODE(asort[p]) == c)
              p += 1;
          }
        parmm[i].abeg = parmm[i-1].aend = p;
        parmm[i].bbeg = parmm[i-1].bend = find_tuple(CODE(asort[p]),bsort,blen);
        parmm[i].kptr = parmx[i].tptr;
      }
    parmm[NTHREADS-1].aend = alen;
    parmm[NTHREADS-1].bend = blen;

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,count_thread,parmm+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

    nhits = parmm[0].nhits;
    for (i = 1; i < NTHREADS; i++)
      parmm[i].nhits = nhits += parmm[i].nhits;

    if (VERBOSE)
      { printf("\nHit count = ");
        Print_Number(nhits,0,stdout);
        if (self || nhits >= blen)
          printf("\nHighwater of %.2fGb space\n",(1. * (Alen + Clen + 2*nhits)) / 134217728);
        else
          printf("\nHighwater of %.2fGb space\n",(1. * (Alen + Clen + blen + nhits)) / 134217728);
        fflush(stdout);
      }

    if (nhits == 0)
      goto zerowork;

    if (self)
      hhit = work1 = (uint64 *) Malloc(sizeof(uint64)*(nhits+1),"Allocating dazzler hit vectors");
    else
      { if (nhits >= blen)
          bsort = (uint64 *) Realloc(bsort,sizeof(uint64)*(nhits+1),
                                       "Reallocating dazzler sort vectors");
        hhit = work1 = bsort;
      }
    khit = work2 = (uint64 *) Malloc(sizeof(uint64)*(nhits+1),"Allocating dazzler hit vectors");
    if (hhit == NULL || khit == NULL || bsort == NULL)
      exit (1);

    MG_blist = bsort;
    MG_hits  = khit;

    for (i = NTHREADS-1; i > 0; i--)
      parmm[i].nhits = parmm[i-1].nhits;
    parmm[0].nhits = 0;

    for (i = 0; i < NTHREADS; i++)
      { parmm[i].kptr = parmx[i].tptr;
        for (p = 0; p < BPOWR; p++)
          parmm[i].kptr[p] = 0;
      }

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,merge_thread,parmm+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

#ifdef TEST_PAIRS
    printf("\nSETUP SORT:\n");
    for (i = 0; i < HOW_MANY && i < nhits; i++)
      { int64 c = khit[i];
        printf(" %5lld / %5lld / %5lld /%5lld\n",
               c >> 48,(c >> 32) & 0xffffll,(c >> 16) & 0xffffll,c & 0xffffll);
      }
#endif
  }

  { int i, x;

    x = 0;
    for (i = 0; i < NTHREADS-1; i++)
      { parmx[i].beg = x;
        parmx[i].end = x = parmm[i+1].nhits;
      }
    parmx[NTHREADS-1].beg = x;
    parmx[NTHREADS-1].end = nhits;

    khit = lex_sort(64,16,khit,hhit,parmx);

    khit[nhits] = 0xffffffffffff0000ll;

#ifdef TEST_CSORT
    printf("\nCROSS SORT:\n");
    for (i = 0; i < HOW_MANY && i <= nhits; i++)
      { uint64 c = khit[i];
        printf(" %5lld / %5lld / %5lld /%5lld\n",
               c >> 48,(c >> 32) & 0xffffll,(c >> 16) & 0xffffll,c & 0xffffll);
      }
#endif
  }

  { int    i, p;
    uint64 d;
    int   *counters;

    MR_ablock = ablock;
    MR_bblock = bblock;
    MR_hits   = khit;
    MR_nhits  = nhits;
    MR_comp   = comp;
    MR_two    = ! self;

    parmr[0].beg  = 0;
    for (i = 1; i < NTHREADS; i++)
      { p = (nhits * ((int64) i)) >> NSHIFT;
        if (p > 0)
          { d = khit[p-1] >> 48;
            while ((khit[p] >> 48) == d)
              p += 1;
          }
        parmr[i].beg  = parmr[i-1].end = p;
      }
    parmr[NTHREADS-1].end = nhits;

    p = (0xffff >> Binshift) + 2;
    counters = (int *) Malloc(NTHREADS*4*p*sizeof(int),"Allocating diagonal buckets");
    if (counters == NULL)
      exit (1);

    for (i = 0; i < 4*p*NTHREADS; i++)
      counters[i] = 0;
    for (i = 0; i < NTHREADS; i++)
      { parmr[i].score = counters + (4*i+1)*p;
        parmr[i].lastp = parmr[i].score + 2*p;
        parmr[i].work  = New_Work_Data();
        parmr[i].aspec = aspec;

        parmr[i].ofile1 =
             Fopen(Catenate(aname,".",bname,Numbered_Suffix((comp?".C":".N"),i,".las")),"w");
        if (parmr[i].ofile1 == NULL)
          exit (1);
        if (self)
          parmr[i].ofile2 = parmr[i].ofile1;
        else
          { parmr[i].ofile2 = 
                Fopen(Catenate(bname,".",aname,Numbered_Suffix((comp?".C":".N"),i,".las")),"w");
            if (parmr[i].ofile2 == NULL)
              exit (1);
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

    if (VERBOSE)
      for (i = 0; i < NTHREADS; i++)
        { nfilt  += parmr[i].nfilt;
          ncheck += parmr[i].ncheck;
        }

    for (i = 0; i < NTHREADS; i++)
      Free_Work_Data(parmr[i].work);
    free(counters);
  }

  free(work2);
  free(work1);

zerowork:

  if (VERBOSE)
    { int width = log10(nhits) + 1;
      width += (width-1)/3;

      printf("\n  ");
      Print_Number(nhits,width,stdout);
      printf(" %d-mers (%e of matrix)\n  ",Kmer,(1.*nhits/atot)/btot);
      Print_Number(nfilt,width,stdout);
      printf(" seed hits (%e of matrix)\n  ",(1.*nfilt/atot)/btot);
      Print_Number(ncheck,width,stdout);
      printf(" confirmed hits (%e of matrix)\n",(1.*ncheck/atot)/btot);
      fflush(stdout);
    }
}
