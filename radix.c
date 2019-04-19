/*******************************************************************************************
 *
 *  Fast threaded lexical sort routine.  Can be compiled to accommodate any element size
 *     (set WORD_SIZE), and makes only n+1 passes to sort n radix bytes.  The radix order
 *     for the bytes of an element may be sorted in any order as listed in the array bytes
 *     (that is -1 terminated).
 *
 *  Author :  Gene Myers
 *  First  :  May 2018
 *
 ********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

#include "radix.h"

typedef unsigned char uint8;
typedef long long     int64;

#undef  TEST_LSORT

static int    NTHREADS;       //  Adjusted downward to nearest power of 2
static int    VERBOSE;        //  Print each byte as it is sorted

void Set_Radix_Params(int nthread, int verbose)
{ NTHREADS = nthread;
  VERBOSE  = verbose;
}

//  Global variables for every "lex_thread"

typedef struct
  { uint8 vector[WORD_SIZE];
  } Element;

static int      LEX_byte;   //  Current byte to sort on
static int      LEX_next;   //  Next byte to sort on (if >= 0)
static int64    LEX_zdiv;   //  Size of thread segments (in bytes)
static Element *LEX_src;    //  Source data goes to ...
static Element *LEX_trg;    //  Target data

//  Thread control record

typedef struct
  { int64  beg;           //  Sort [beg,end) of LEX_src
    int64  end;
    int64  tptr[256];     //  Finger for each 8-bit value
    int64 *sptr;          //  Conceptually [256][NTHREADS].  At end of sorting pass
  } Lex_Arg;              //    sprtr[b][n] = # of occurences of value b in rangd of
                          //    thread n for the *next* pass

//  Threaded sorting pass

static void *lex_thread(void *arg)
{ Lex_Arg    *data  = (Lex_Arg *) arg;
  int64      *sptr  = data->sptr;
  int64      *tptr  = data->tptr;
  Element    *src   = LEX_src;
  Element    *trg   = LEX_trg;
  int64       zdiv  = LEX_zdiv;

  int64       i, n, x;
  uint8      *q;

  n = data->end;
  if (LEX_next < 0)
    for (i = data->beg; i < n; i++)
      { x = tptr[src[i].vector[LEX_byte]]++;
        trg[x] = src[i];
      }
  else
    for (i = data->beg; i < n; i++)
      { q = src[i].vector;
        x = tptr[q[LEX_byte]]++;
        trg[x] = src[i];
        sptr[((x/zdiv) << 8) | q[LEX_next]] += 1;
      }
  return  (NULL);
}

//  Threaded sort initiation pass: count bucket sizes

static void *lexbeg_thread(void *arg)
{ Lex_Arg    *data  = (Lex_Arg *) arg;
  int64      *tptr  = data->tptr;
  Element    *src   = LEX_src;

  int64       i, n;

  n = data->end;
  for (i = data->beg; i < n; i++)
    tptr[src[i].vector[LEX_byte]] += 1;
  return (NULL);
}

//  Radix sort the indicated "bytes" of src, using array trg as the secondary array
//    The arrays contains len elements each of "size" bytes.
//    Return a pointer to the array containing the final result.

void *Radix_Sort(int64 len, void *src, void *trg, int *bytes)
{ pthread_t threads[NTHREADS];
  Lex_Arg   parmx[NTHREADS];   //  Thread control record for sorting

  Element *xch;
  int64    x, y;
  int      i, j, z, b;

  LEX_zdiv = (len-1)/NTHREADS + 1;
  LEX_src   = (Element *) src;
  LEX_trg   = (Element *) trg;

  for (i = 0; i < NTHREADS; i++)
    parmx[i].sptr = (int64 *) alloca(NTHREADS*256*sizeof(int64));

  //  For each requested byte b in order, radix sort

  for (b = 0; bytes[b] >= 0; b++)
    { LEX_byte  = bytes[b];
      LEX_next  = bytes[b+1];

      if (VERBOSE)
        { printf("     Sorting byte %d\n",LEX_byte);
          fflush(stdout);
        }

      //  Setup beg, end, and zero tptr counters

      x = 0;
      for (i = 0; i < NTHREADS; i++)
        { parmx[i].beg = x;
          x = LEX_zdiv*(i+1);
          if (x > len)
            x = len;
          parmx[i].end = x;
          for (j = 0; j < 256; j++)
            parmx[i].tptr[j] = 0;
        }
      parmx[NTHREADS-1].end = len;

      //  If first pass, then explicitly sweep to get tptr counts
      //    otherwise accumulate from sptr counts of last sweep

      if (b == 0)
        { for (i = 0; i < NTHREADS; i++)
            pthread_create(threads+i,NULL,lexbeg_thread,parmx+i);

          for (i = 0; i < NTHREADS; i++)
            pthread_join(threads[i],NULL);
        }
      else
        { int64 *pxt, *pxs;
          for (i = 0; i < NTHREADS; i++)
            { pxt = parmx[i].tptr;
              for (z = 0; z < NTHREADS; z++)
                { pxs = parmx[z].sptr + (i<<8);
                  for (j = 0; j < 256; j++)
                    pxt[j] += pxs[j];
                }
            }
        }

#ifdef TEST_LSORT
      printf(" Next counts\n");
      for (j = 0; j < 256; j++)
        { printf(" %3d:",j);
          for (i = 0; i < NTHREADS; i++)
            printf(" %5lld",parmx[i].tptr[j]);
          printf("\n");
        }
#endif

      //   Zero sptr array counters in preparation of pass

      for (i = 0; i < NTHREADS; i++)
        for (z = (NTHREADS<<8)-1; z >= 0; z--)
          parmx[i].sptr[z] = 0;

      //  Convert tptr from counts to fingers.

      x = 0;
      for (j = 0; j < 256; j++)
        for (i = 0; i < NTHREADS; i++)
          { y = parmx[i].tptr[j];
            parmx[i].tptr[j] = x;
            x += y;
          }

      //  Threaded pass

      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,lex_thread,parmx+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);

      xch     = LEX_src;
      LEX_src = LEX_trg;
      LEX_trg = xch;

#ifdef TEST_LSORT
      printf("\nLSORT %d\n",LEX_byte);
      for (i = 0; i < 1000; i++)
        { printf(" %4d: ",i);
          for (j = 0; j < WORD_SIZE; j++)
            printf("%02x",LEX_src[i].vector[j]);
          printf("\n");
        }
#endif
    }

  return (LEX_src);
}
