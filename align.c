/************************************************************************************\
*                                                                                    *
* Copyright (c) 2014, Dr. Eugene W. Myers (EWM). All rights reserved.                *
*                                                                                    *
* Redistribution and use in source and binary forms, with or without modification,   *
* are permitted provided that the following conditions are met:                      *
*                                                                                    *
*  · Redistributions of source code must retain the above copyright notice, this     *
*    list of conditions and the following disclaimer.                                *
*                                                                                    *
*  · Redistributions in binary form must reproduce the above copyright notice, this  *
*    list of conditions and the following disclaimer in the documentation and/or     *
*    other materials provided with the distribution.                                 *
*                                                                                    *
*  · The name of EWM may not be used to endorse or promote products derived from     *
*    this software without specific prior written permission.                        *
*                                                                                    *
* THIS SOFTWARE IS PROVIDED BY EWM ”AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES,    *
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND       *
* FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL EWM BE LIABLE   *
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS  *
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY      *
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING     *
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN  *
* IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                      *
*                                                                                    *
* For any issues regarding this software and its use, contact EWM at:                *
*                                                                                    *
*   Eugene W. Myers Jr.                                                              *
*   Bautzner Str. 122e                                                               *
*   01099 Dresden                                                                    *
*   GERMANY                                                                          *
*   Email: gene.myers@gmail.com                                                      *
*                                                                                    *
\************************************************************************************/

/*******************************************************************************************
 *
 *  Fast alignment discovery and trace generation along with utilites for displaying alignments
 *     Based on previously unpublished ideas from 2005, subsequently refined in 2013-14.  Basic
 *     idea is to keep a dynamically selected interval of the f.r. waves from my 1986 O(nd) paper.
 *     A recent cool idea is to not record all the details of an alignment while discovering it
 *     but simply record trace points through which the optimal alignment passes every 100bp,
 *     allowing rapid recomputation of the alignment details between trace points.
 *
 *  Author :  Gene Myers
 *  First  :  June 2013
 *  Current:  June 1, 2014
 *
 ********************************************************************************************/

// align1: Derived from the original BOA aligner

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>

#include "DB.h"
#include "align.h"

#undef    DEBUG_PASSES     //  Show forward / backward extension termini for Local_Alignment
#undef    DEBUG_POINTS     //  Show trace points
#undef    DEBUG_WAVE       //  Show waves of Local_Alignment
#undef     SHOW_MATCH_WAVE //  For waves of Local_Alignment also show # of matches
#undef    SHOW_TRAIL       //  Show trace at the end of forward and reverse passes
#undef    SHOW_TPS         //  Show trace points as they are encountered in a wave

#undef  DEBUG_EXTEND       //  Show waves of Extend_Until_Overlap

#undef  DEBUG_ALIGN        //  Show division points of Compute_Trace
#undef  DEBUG_SCRIPT       //  Show trace additions for Compute_Trace
#undef  DEBUG_AWAVE        //  Show F/R waves of Compute_Trace
#define   SMALL_BIT 100

#undef  SHOW_TRACE         //  Show full trace for Print_Alignment

#undef  WAVE_STATS


/****************************************************************************************\
*                                                                                        *
*  Working Storage Abstraction                                                           *
*                                                                                        *
\****************************************************************************************/

typedef struct            //  Hidden from the user, working space for each thread
  { int     vecmax;
    void   *vector;
    int     celmax;
    void   *cells;
    int     pntmax;
    void   *points;
    int     tramax;
    void   *trace;
  } _Work_Data;

Work_Data *New_Work_Data()
{ _Work_Data *work;
  
  work = (_Work_Data *) Malloc(sizeof(_Work_Data),"Allocating work data block");
  if (work == NULL)
    exit (1);
  work->vecmax = 0;
  work->vector = NULL;
  work->pntmax = 0;
  work->points = NULL;
  work->tramax = 0;
  work->trace  = NULL;
  work->celmax = 0;
  work->cells  = NULL;
  return ((Work_Data *) work);
}

static void enlarge_vector(_Work_Data *work, int newmax)
{ work->vecmax = ((int) (newmax*1.2)) + 10000;
  work->vector = Realloc(work->vector,work->vecmax,"Enlarging DP vector");
  if (work->vector == NULL)
    exit (1);
}

static void enlarge_points(_Work_Data *work, int newmax)
{ work->pntmax = ((int) (newmax*1.2)) + 10000;
  work->points = Realloc(work->points,work->pntmax,"Enlarging point vector");
  if (work->points == NULL)
    exit (1);
}

static void enlarge_trace(_Work_Data *work, int newmax)
{ work->tramax = ((int) (newmax*1.2)) + 10000;
  work->trace  = Realloc(work->trace,work->tramax,"Enlarging trace vector");
  if (work->trace == NULL)
    exit (1);
}

void Free_Work_Data(Work_Data *ework)
{ _Work_Data *work = (_Work_Data *) ework;
  if (work->vector != NULL)
    free(work->vector);
  if (work->cells != NULL)
    free(work->cells);
  if (work->trace != NULL)
    free(work->trace);
  if (work->points != NULL)
    free(work->points);
  free(work);
}


/****************************************************************************************\
*                                                                                        *
*  ADAPTIVE PATH FINDING                                                                 *
*                                                                                        *
\****************************************************************************************/

  //  Absolute/Fixed Parameters

#define BVEC  uint64     //  Can be uint32 if PATH_LEN <= 32

#define TRIM_LEN    15   //  Report as the tip, the last wave maximum for which the last
                         //     2*TRIM_LEN edits are prefix-positive at rate ave_corr*f(bias)
                         //     (max value is 20)

#define PATH_LEN    60   //  Follow the last PATH_LEN columns/edges (max value is 63)

  //  Derivative fixed parameters

#define PATH_TOP  0x1000000000000000ll   //  Must be 1 << PATH_LEN
#define PATH_INT  0x0fffffffffffffffll   //  Must be PATH_TOP-1
#define TRIM_MASK 0x7fff                 //  Must be (1 << TRIM_LEN) - 1
#define TRIM_MLAG 200                    //  How far can last trim point be behind best point
#define WAVE_LAG   30                    //  How far can worst point be behind the best point

static double Bias_Factor[10] = { .690, .690, .690, .690, .780,
                                  .850, .900, .933, .966, 1.000 };

  //  Adjustable paramters

typedef struct
  { double ave_corr;
    int    trace_space;
    float  freq[4];
    int    ave_path;
    int16 *score;
    int16 *table;
  } _Align_Spec;
 
/* Fill in bit table: TABLE[x] = 1 iff the alignment modeled by x (1 = match, 0 = mismatch)
     has a non-negative score for every suffix of the alignment under the scoring scheme
     where match = MATCH and mismatch = -1.  MATCH is set so that an alignment with TRIM_PCT
     matches has zero score ( (1-TRIM_PCT) / TRIM_PCT ).                                     */

#define FRACTION 1000  //  Implicit fractional part of scores, i.e. score = x/FRACTION

typedef struct
  { int    mscore;
    int    dscore;
    int16 *table;
    int16 *score;
  } Table_Bits;

static void set_table(int bit, int prefix, int score, int max, Table_Bits *parms)
{ if (bit >= TRIM_LEN)
    { parms->table[prefix] = (int16) (score-max);
      parms->score[prefix] = (int16) score;
    }
  else
    { if (score > max)
        max = score;
      set_table(bit+1,(prefix<<1),score - parms->dscore,max,parms);
      set_table(bit+1,(prefix<<1) | 1,score + parms->mscore,max,parms);
    }
}

/* Create an alignment specification record including path tip tables & values */

Align_Spec *New_Align_Spec(double ave_corr, int trace_space, float *freq)
{ _Align_Spec *spec;
  Table_Bits   parms;
  double       match;
  int          bias;

  spec = (_Align_Spec *) Malloc(sizeof(_Align_Spec),"Allocating alignment specification");
  if (spec == NULL)
    exit (1);

  spec->ave_corr    = ave_corr;
  spec->trace_space = trace_space;
  spec->freq[0]     = freq[0];
  spec->freq[1]     = freq[1];
  spec->freq[2]     = freq[2];
  spec->freq[3]     = freq[3];

  match = freq[0] + freq[3];
  if (match > .5)
    match = 1.-match;
  bias = (int) ((match+.025)*20.-1.);
  if (match < .2)
    { fprintf(stderr,"Warning: Base bias worse than 80/20%% !\n");
      bias = 3;
    }

  spec->ave_path = (int) (PATH_LEN * (1. - Bias_Factor[bias] * (1. - ave_corr)));
  parms.mscore   = (int) (FRACTION * Bias_Factor[bias] * (1. - ave_corr));
  parms.dscore   = FRACTION - parms.mscore;

  parms.score = (int16 *) Malloc(sizeof(int16)*(TRIM_MASK+1)*2,"Allocating trim table");
  if (parms.score == NULL)
    exit (1);
  parms.table = parms.score + (TRIM_MASK+1);

  set_table(0,0,0,0,&parms);

  spec->table = parms.table;
  spec->score = parms.score;

  return ((Align_Spec *) spec);
}

void Free_Align_Spec(Align_Spec *espec)
{ _Align_Spec *spec = (_Align_Spec *) espec;
  free(spec->score);
  free(spec);
}

double Average_Correlation(Align_Spec *espec)
{ return (((_Align_Spec *) espec)->ave_corr); }

int Trace_Spacing(Align_Spec *espec)
{ return (((_Align_Spec *) espec)->trace_space); }

float *Base_Frequencies(Align_Spec *espec)
{ return (((_Align_Spec *) espec)->freq); }


/****************************************************************************************\
*                                                                                        *
*  LOCAL ALIGNMENT FINDER: forward_/reverse_wave and Local_Alignment                     *
*                                                                                        *
\****************************************************************************************/


#ifdef WAVE_STATS

static int64 MAX, TOT, NWV;
static int64 RESTARTS;

void Init_Stats()
{ MAX = TOT = NWV = 0;
  RESTARTS = 0;
}

void Print_Stats()
{ printf("\nMax = %lld  Ave = %.1f  # = %lld\n",MAX,(1.*TOT)/NWV,NWV);
  printf("\nRestarts = %lld\n",RESTARTS);
}

#endif


#ifdef DEBUG_WAVE

static void print_wave(int *V, int *M, int low, int hgh, int besta)
{ int k, bestk;

  (void) M;
  printf("  [%6d,%6d]: ",low,hgh);
  for (k = low; k <= hgh; k++)
    { if (besta == V[k])
        bestk = k;
      // printf(" %3d",(V[k]+k)/2);
      printf(" %3d",besta-V[k]);
    }
  printf(" : %d (%d,%d)\n",besta,(besta+bestk)/2,(besta-bestk)/2);
#ifdef SHOW_MATCH_WAVE
  printf("                   ");
  for (k = low; k <= hgh; k++)
    printf(" %3d",M[k]);
  printf("\n");
#endif
  fflush(stdout);
}

#endif

/* At each furthest reaching point, keep a-coordinate of point (V), bitvector
     recording the last TRIM_LEN columns of the implied alignment (T), and the
     # of matches (1-bits) in the bitvector (M).                               */

typedef struct
  { int ptr;
    int diag;
    int mark;
  } Pebble;

static void forward_wave(_Work_Data *work, _Align_Spec *spec,
                         Alignment *align, Path *bpath,
                         int mind, int maxd, int mida)
{ char *aseq  = align->aseq;
  char *bseq  = align->bseq;
  Path *apath = align->path;

  int     hgh, low, dif, pos;
  int    *V, *M;
  BVEC   *T;

  int    *HA, *HB;
  int    *NA, *NB;
  Pebble *cells;
  int     avail, cmax, boff;

  int     TRACE_SPACE = spec->trace_space;
  int     PATH_AVE    = spec->ave_path;
  int16  *SCORE       = spec->score;
  int16  *TABLE       = spec->table;

  int     besta, besty;
  int     trima, trimy, trimd;
  int     trimha, trimhb;
  int     morea, morey, mored;
  int     moreha, morehb;
  int     more, morem, lasta;
  int     aclip, bclip;

  { int alen = align->alen + 1;
    int blen = align->blen + 1;
    int tlen = alen + blen + 1;

    V  = ((int *) work->vector) + blen;
    M  = V + tlen;
    HA = M + tlen;
    HB = HA + tlen;
    NA = HB + tlen;
    NB = NA + tlen;
    T  = ((BVEC *) (NB + alen)) + blen;

    cells = (Pebble *) (work->cells);
    cmax  = work->celmax;
    avail = 0;

    if (COMP(align->flags))
      boff = align->blen % TRACE_SPACE;
    else
      boff = 0;
  }

  /* Compute 0-wave starting from mid-line */

  hgh   = maxd;
  low   = mind;
  if (aseq == bseq)
    pos = 1;
  else
    pos = -INT32_MAX;
  dif   = 0;

  more  = 1;
  aclip =  INT32_MAX;
  bclip = -INT32_MAX;

  besta  = trima  = morea = lasta = mida;
  besty  = trimy  = morey = (mida-hgh) >> 1;
  trimd  = mored  = 0;
  trimha = moreha = 0;
  trimhb = morehb = 1;
  morem  = -1;

  { int   k;
    char *a;

    a  = aseq + hgh;
    for (k = hgh; k >= low; k--)
      { int     y, c, d;
        int     ha, hb;
        int     na, nb;
        Pebble *pb;

        y = (mida-k) >> 1;

        if (avail >= cmax-1)
          { cmax  = ((int) (avail*1.2)) + 10000;
            cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
            if (cells == NULL)
              exit (1);
          }

        na = ((y+k)/TRACE_SPACE)*TRACE_SPACE;
#ifdef SHOW_TPS
        printf(" A %d: %d,%d,0,%d\n",avail,-1,k,na); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->mark = na;
        ha  = avail++;
        na += TRACE_SPACE;

        nb = ((y+(TRACE_SPACE-boff))/TRACE_SPACE-1)*TRACE_SPACE+boff;
#ifdef SHOW_TPS
        printf(" B %d: %d,%d,0,%d\n",avail,-1,k,nb); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->mark = nb;
        hb  = avail++;
        nb += TRACE_SPACE;

        while (1)
          { c = bseq[y];
            if (c == 4)
              { more = 0;
                if (bclip < k)
                  bclip = k;
                break;
              }
            d = a[y];
            if (c != d)
              { if (d == 4) 
                  { more  = 0;
                    aclip = k;
                  }
                break;
              }
            y += 1;
          }
        c = (y << 1) + k;

        while (y+k >= na)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  exit (1);
              }
#ifdef SHOW_TPS
            printf(" A %d: %d,%d,0,%d\n",avail,ha,k,na); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = ha;
            pb->diag = k;
            pb->mark = na;
            ha  = avail++;
            na += TRACE_SPACE;
          }
        while (y >= nb)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  exit (1);
              }
#ifdef SHOW_TPS
            printf(" B %d: %d,%d,0,%d\n",avail,hb,k,nb); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = hb;
            pb->diag = k;
            pb->mark = nb;
            hb  = avail++;
            nb += TRACE_SPACE;
          }

        if (c > besta)
          { besta  = trima = lasta = c;
            besty  = trimy = y;
            trimha = ha;
            trimhb = hb;
          }

        V[k]  = c;
        T[k]  = PATH_INT;
        M[k]  = PATH_LEN;
        HA[k] = ha;
        HB[k] = hb;
        NA[k] = na;
        NB[k] = nb;

        a -= 1;
      }
  }

  if (more == 0)
    { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
        more = 1;
      if (hgh >= aclip)
        { hgh = aclip-1;
          if (morem <= M[aclip])
            { morem  = M[aclip];
              morea  = V[aclip];
              morey  = (morea - aclip)/2;
              moreha = HA[aclip];
              morehb = HB[aclip];
            }
        }
      if (low <= bclip)
        { low = bclip+1;
          if (morem <= M[bclip])
            { morem  = M[bclip];
              morea  = V[bclip];
              morey  = (morea - bclip)/2;
              moreha = HA[bclip];
              morehb = HB[bclip];
            }
        }
      aclip =  INT32_MAX;
      bclip = -INT32_MAX;
    }

#ifdef DEBUG_WAVE
  printf("\nFORWARD WAVE:\n");
  print_wave(V,M,low,hgh,besta);
#endif

  /* Compute successive waves until no furthest reaching points remain */

  while (more && lasta >= besta - TRIM_MLAG)
    { int     k, n;
      int     ua, ub;
      BVEC    t;
      int     am, ac, ap;
      char   *a;

      hgh += 1;
      if (low > pos)
        { low -= 1;
          NA[low] = NA[low+1];
          NB[low] = NB[low+1];
          V[low]  = -1;
        }
      dif += 1;
      NA[hgh] = NA[hgh-1];
      NB[hgh] = NB[hgh-1];

      am = ac = V[hgh] = V[hgh+1] = V[low-1] = -1;
      a  = aseq + hgh;
      t  = PATH_INT;
      n  = PATH_LEN;
      ua = ub = -1;
      for (k = hgh; k >= low; k--)
        { int     y, m;
          int     ha, hb;
          int     c, d;
          BVEC    b;
          Pebble *pb;

          ap = ac;
          ac = am;
          am = V[d = k-1];

          if (ac < am)
            if (am < ap)
              { c  = ap+1;
                m  = n;
                b  = t;
                ha = ua;
                hb = ub;
              }
            else
              { c  = am+1;
                m  = M[d];
                b  = T[d]; 
                ha = HA[d];
                hb = HB[d];
              }
          else
            if (ac < ap)
              { c  = ap+1;
                m  = n;
                b  = t;
                ha = ua;
                hb = ub;
              }
            else
              { c  = ac+2;
                m  = M[k];
                b  = T[k];
                ha = HA[k];
                hb = HB[k];
              }

          if ((b & PATH_TOP) != 0)
            m -= 1;
          b <<= 1;

          y = (c-k) >> 1;
          while (1)
            { c = bseq[y];
              if (c == 4)
                { more = 0;
                  if (bclip < k)
                    bclip = k;
                  break;
                }
              d = a[y];
              if (c != d)
                { if (d == 4) 
                    { more  = 0;
                      aclip = k;
                    }
                  break;
                }
              y += 1;
              if ((b & PATH_TOP) == 0)
                m += 1;
              b = (b << 1) | 1;
            }
          c = (y << 1) + k;

          while (y+k >= NA[k])
            { if (cells[ha].mark < NA[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        exit (1);
                    }
#ifdef SHOW_TPS
                  printf(" A %d: %d,%d,%d,%d\n",avail,ha,k,dif,NA[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = ha;
                  pb->diag = k;
                  pb->mark = NA[k];
                  ha = avail++;
                }
              NA[k] += TRACE_SPACE;
            }

          while (y >= NB[k])
            { if (cells[hb].mark < NB[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        exit (1);
                    }
#ifdef SHOW_TPS
                  printf(" B %d: %d,%d,%d,%d\n",avail,hb,k,dif,NB[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = hb;
                  pb->diag = k;
                  pb->mark = NB[k];
                  hb = avail++;
                }
              NB[k] += TRACE_SPACE;
            }

          if (c > besta)
            { besta = c;
              besty = y;
              if (m >= PATH_AVE)
                { lasta = c;
                  if (TABLE[b & TRIM_MASK] >= 0)
                    if (TABLE[(b >> TRIM_LEN) & TRIM_MASK] + SCORE[b & TRIM_MASK] >= 0)
                      { trima  = c;
                        trimy  = y;
                        trimd  = dif;
                        trimha = ha;
                        trimhb = hb;
                      }
                }
            }

          t  = T[k];
          n  = M[k];
          ua = HA[k];
          ub = HB[k];
          V[k]  = c;
          T[k]  = b;
          M[k]  = m;
          HA[k] = ha;
          HB[k] = hb;

          a -= 1;
        }

      if (more == 0)
        { if (bseq[besty] != 4 && aseq[besta-besty] != 4)
            more = 1;
          if (hgh >= aclip)
            { hgh = aclip-1;
              if (morem <= M[aclip])
                { morem  = M[aclip];
                  morea  = V[aclip];
                  morey  = (morea - aclip)/2;
                  mored  = dif;
                  moreha = HA[aclip];
                  morehb = HB[aclip];
                }
            }
          if (low <= bclip)
            { low = bclip+1;
              if (morem <= M[bclip])
                { morem  = M[bclip];
                  morea  = V[bclip];
                  morey  = (morea - bclip)/2;
                  mored  = dif;
                  moreha = HA[bclip];
                  morehb = HB[bclip];
                }
            }
          aclip =  INT32_MAX;
          bclip = -INT32_MAX;
        }

      n = besta - WAVE_LAG;
      while (hgh >= low)
        if (V[hgh] < n)
          hgh -= 1;                               
        else
          { while (V[low] < n)
              low += 1;
            break;
          }

#ifdef WAVE_STATS
      k = (hgh-low)+1;
      if (k > MAX)
        MAX = k;
      TOT += k;
      NWV += 1;
#endif

#ifdef DEBUG_WAVE
      print_wave(V,M,low,hgh,besta);
#endif
    }

  { uint16 *atrace = (uint16 *) apath->trace;
    uint16 *btrace = (uint16 *) bpath->trace;
    int     atlen, btlen;
    int     trimx;
    int     a, b, k, h;

    if (morem >= 0)
      { trimx  = morea-morey;
        trimy  = morey;
        trimd  = mored;
        trimha = moreha;
        trimhb = morehb;
      }
    else
      trimx = trima-trimy;

    atlen = btlen = 0;

    a = -1;
    for (h = trimha; h >= 0; h = b)
      { b = cells[h].ptr; 
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = (mida-k)/2;
#ifdef SHOW_TRAIL
    printf("  A path = (%5d,%5d)\n",(mida+k)/2,b); fflush(stdout);
#endif
    for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
      { k = cells[h].diag;
        a = cells[h].mark - k;
        atrace[atlen++] = (uint16) (a-b);
#ifdef SHOW_TRAIL
        printf("     %4d: (%5d,%5d): %3d\n",h,a+k,a,a-b); fflush(stdout);
#endif
        b = a;
      }
    if (b+k != trimx)
      { atrace[atlen++] = (uint16) (trimy-b);
#ifdef SHOW_TRAIL
        printf("           (%5d,%5d): %3d\n",trimx,trimy,trimy-b); fflush(stdout);
#endif
      }
    else if (b != trimy)
      { atrace[atlen-1] = (uint16) (atrace[atlen-1] + (trimy-b));
#ifdef SHOW_TRAIL
        printf("         @ (%5d,%5d): %3d\n",trimx,trimy,trimy-b); fflush(stdout);
#endif
      }

    a = -1;
    for (h = trimhb; h >= 0; h = b)
      { b = cells[h].ptr; 
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = (mida+k)/2;
#ifdef SHOW_TRAIL
    printf("  B path = (%5d,%5d)\n",b,(mida-k)/2); fflush(stdout);
#endif
    for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
      { k = cells[h].diag;
        a = cells[h].mark + k;
        btrace[btlen++] = (uint16) (a-b);
#ifdef SHOW_TRAIL
        printf("     %4d: (%5d,%5d): %3d\n",h,a,a-k,a-b); fflush(stdout);
#endif
        b = a;
      }
    if (b-k != trimy)
      { btrace[btlen++] = (uint16) (trimx-b);
#ifdef SHOW_TRAIL
        printf("           (%5d,%5d): %3d\n",trimx,trimy,trimx-b); fflush(stdout);
#endif
      }
    else if (b != trimx)
      { btrace[btlen-1] = (uint16) (btrace[btlen-1] + (trimx-b));
#ifdef SHOW_TRAIL
        printf("         @ (%5d,%5d): %3d\n",trimx,trimy,trimx-b); fflush(stdout);
#endif
      }

    apath->aepos = (READIDX) trimx;
    apath->bepos = (READIDX) trimy;
    apath->diffs = (READIDX) trimd;
    apath->tlen  = (READIDX) atlen;
    if (COMP(align->flags))
      { bpath->abpos = (READIDX) (align->blen - apath->bepos);
        bpath->bbpos = (READIDX) (align->alen - apath->aepos);
      }
    else
      { bpath->aepos = apath->bepos;
        bpath->bepos = apath->aepos;
      }
    bpath->diffs = (READIDX) trimd;
    bpath->tlen  = (READIDX) btlen;
  }

  work->cells  = (void *) cells;
  work->celmax = cmax;
}

/*** Reverse Wave ***/

static void reverse_wave(_Work_Data *work, _Align_Spec *spec,
                         Alignment *align, Path *bpath, int mind, int maxd, int mida)
{ char *aseq  = align->aseq - 1;
  char *bseq  = align->bseq - 1;
  Path *apath = align->path;

  int     hgh, low, dif, pos;
  int    *V, *M;
  BVEC   *T;

  int    *HA, *HB;
  int    *NA, *NB;
  Pebble *cells;
  int     avail, cmax, boff;

  int     TRACE_SPACE = spec->trace_space;
  int     PATH_AVE    = spec->ave_path;
  int16  *SCORE       = spec->score;
  int16  *TABLE       = spec->table;

  int     besta, besty;
  int     trima, trimy, trimd;
  int     trimha, trimhb;
  int     morea, morey, mored;
  int     moreha, morehb;
  int     more, morem, lasta;
  int     aclip, bclip;

  { int alen = align->alen + 1;
    int blen = align->blen + 1;
    int tlen = alen + blen + 1;

    V  = ((int *) work->vector) + blen;
    M  = V + tlen;
    HA = M + tlen;
    HB = HA + tlen;
    NA = HB + tlen;
    NB = NA + tlen;
    T  = ((BVEC *) (NB + alen)) + blen;

    cells = (Pebble *) (work->cells);
    cmax  = work->celmax;
    avail = 0;

    if (COMP(align->flags))
      boff = align->blen % TRACE_SPACE;
    else
      boff = 0;
  }

  hgh   = maxd;
  low   = mind;
  if (aseq == bseq)
    pos = 1;
  else
    pos = -INT32_MAX;
  dif   = 0;

  more  = 1;
  aclip = -INT32_MAX;
  bclip =  INT32_MAX;

  besta  = trima  = morea = lasta = mida;
  besty  = trimy  = morey = (mida-hgh) >> 1;
  trimd  = mored  = 0;
  trimha = moreha = 0;
  trimhb = morehb = 1;
  morem  = -1;

  { int   k;
    char *a;

    a = aseq + low;
    for (k = low; k <= hgh; k++)
      { int     y, c, d;
        int     ha, hb;
        int     na, nb;
        Pebble *pb;

        y = (mida-k) >> 1;

        if (avail >= cmax-1)
          { cmax  = ((int) (avail*1.2)) + 10000;
            cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
            if (cells == NULL)
              exit (1);
          }

        na = ((y+k-1)/TRACE_SPACE)*TRACE_SPACE;
#ifdef SHOW_TPS
        printf(" A %d: -1,%d,0,%d\n",avail,k,na+TRACE_SPACE); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->mark = y+k;
        ha  = avail++;

        nb = ((y+(TRACE_SPACE-boff)-1)/TRACE_SPACE-1)*TRACE_SPACE+boff;
#ifdef SHOW_TPS
        printf(" B %d: -1,%d,0,%d\n",avail,k,nb+TRACE_SPACE); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->mark = y;
        hb  = avail++;

        while (1)
          { c = bseq[y];
            if (c == 4)
              { more = 0;
                if (bclip > k)
                  bclip = k;
                break;
              }
            d = a[y];
            if (c != d)
              { if (d == 4) 
                  { more  = 0;
                    aclip = k;
                  }
                break;
              }
            y -= 1;
          }
        c = (y << 1) + k;

        while (y+k <= na)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  exit (1);
              }
#ifdef SHOW_TPS
            printf(" A %d: %d,%d,0,%d\n",avail,ha,k,na); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = ha;
            pb->diag = k;
            pb->mark = na;
            ha  = avail++;
            na -= TRACE_SPACE;
          }
        while (y <= nb)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  exit (1);
              }
#ifdef SHOW_TPS
            printf(" B %d: %d,%d,0,%d\n",avail,hb,k,nb); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = hb;
            pb->diag = k;
            pb->mark = nb;
            hb  = avail++;
            nb -= TRACE_SPACE;
          }

        if (c < besta)
          { besta  = trima = lasta = c;
            besty  = trimy = y;
            trimha = ha;
            trimhb = hb;
          }

        V[k]  = c;
        T[k]  = PATH_INT;
        M[k]  = PATH_LEN;
        HA[k] = ha;
        HB[k] = hb;
        NA[k] = na;
        NB[k] = nb;

        a += 1;
      }
  }

  if (more == 0)
    { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
        more = 1;
      if (low <= aclip)
        { low = aclip+1;
          if (morem <= M[aclip])
            { morem  = M[aclip];
              morea  = V[aclip];
              morey  = (morea - aclip)/2;
              moreha = HA[aclip];
              morehb = HB[aclip];
            }
        }
      if (hgh >= bclip)
        { hgh = bclip-1;
          if (morem <= M[bclip])
            { morem  = M[bclip];
              morea  = V[bclip];
              morey  = (morea - bclip)/2;
              moreha = HA[bclip];
              morehb = HB[bclip];
            }
        }
      aclip = -INT32_MAX;
      bclip =  INT32_MAX;
    }

#ifdef DEBUG_WAVE
  printf("\nREVERSE WAVE:\n");
  print_wave(V,M,low,hgh,besta);
#endif

  while (more && lasta <= besta + TRIM_MLAG)
    { int    k, n;
      int    ua, ub;
      BVEC   t;
      int    am, ac, ap;
      char  *a;

      hgh += 1;
      if (low > pos)
        { low -= 1;
          NA[low] = NA[low+1];
          NB[low] = NB[low+1];
          V[low]  = ap = INT32_MAX;
        }
      else
        ap = V[low]; 
      dif += 1;
      NA[hgh] = NA[hgh-1];
      NB[hgh] = NB[hgh-1];

      ac = V[hgh] = V[hgh+1] = V[low-1] = INT32_MAX;
      a  = aseq + low;
      t  = PATH_INT;
      n  = PATH_LEN;
      ua = ub = -1;
      for (k = low; k <= hgh; k++)
        { int     y, m;
          int     ha, hb;
          int     c, d;
          BVEC    b;
          Pebble *pb;

          am = ac;
          ac = ap;
          ap = V[d = k+1];

          if (ac > ap)
            if (ap > am)
              { c = am-1;
                m  = n;
                b  = t;
                ha = ua;
                hb = ub;
              }
            else
              { c  = ap-1;
                m  = M[d];
                b  = T[d];
                ha = HA[d];
                hb = HB[d];
              }
          else
            if (ac > am)
              { c  = am-1;
                m  = n;
                b  = t;
                ha = ua;
                hb = ub;
              }
            else
              { c  = ac-2;
                m  = M[k];
                b  = T[k];
                ha = HA[k];
                hb = HB[k];
              }

          if ((b & PATH_TOP) != 0)
            m -= 1;
          b <<= 1;

          y = (c-k) >> 1;
          while (1)
            { c = bseq[y];
              if (c == 4)
                { more = 0;
                  if (bclip > k)
                    bclip = k;
                  break;
                }
              d = a[y];
              if (c != d)
                { if (d == 4) 
                    { more  = 0;
                      aclip = k;
                    }
                  break;
                }
              y -= 1;
              if ((b & PATH_TOP) == 0)
                m += 1;
              b = (b << 1) | 1;
            }
          c = (y << 1) + k;

          while (y+k <= NA[k])
            { if (cells[ha].mark > NA[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        exit (1);
                    }
#ifdef SHOW_TPS
                  printf(" A %d: %d,%d,%d,%d\n",avail,ha,k,dif,NA[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = ha;
                  pb->diag = k;
                  pb->mark = NA[k];
                  ha = avail++;
                }
              NA[k] -= TRACE_SPACE;
            }
          while (y <= NB[k])
            { if (cells[hb].mark > NB[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        exit (1);
                    }
#ifdef SHOW_TPS
                  printf(" B %d: %d,%d,%d,%d\n",avail,hb,k,dif,NB[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = hb;
                  pb->diag = k;
                  pb->mark = NB[k];
                  hb = avail++;
                }
              NB[k] -= TRACE_SPACE;
            }

          if (c < besta)
            { besta = c;
              besty = y;
              if (m >= PATH_AVE)
                { lasta = c;
                  if (TABLE[b & TRIM_MASK] >= 0)
                    if (TABLE[(b >> TRIM_LEN) & TRIM_MASK] + SCORE[b & TRIM_MASK] >= 0)
                      { trima  = c;
                        trimy  = y;
                        trimd  = dif;
                        trimha = ha;
                        trimhb = hb;
                      }
                }
            }

          t  = T[k];
          n  = M[k];
          ua = HA[k];
          ub = HB[k];
          V[k]  = c;
          T[k]  = b;
          M[k]  = m;
          HA[k] = ha;
          HB[k] = hb;

          a += 1;
        }

      if (more == 0)
        { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
            more = 1;
          if (low <= aclip)
            { low = aclip+1;
              if (morem <= M[aclip])
                { morem  = M[aclip];
                  morea  = V[aclip];
                  morey  = (morea - aclip)/2;
                  mored  = dif;
                  moreha = HA[aclip];
                  morehb = HB[aclip];
                }
            }
          if (hgh >= bclip)
            { hgh = bclip-1;
              if (morem <= M[bclip])
                { morem  = M[bclip];
                  morea  = V[bclip];
                  morey  = (morea - bclip)/2;
                  mored  = dif;
                  moreha = HA[bclip];
                  morehb = HB[bclip];
                }
            }
          aclip = -INT32_MAX;
          bclip =  INT32_MAX;
        }

      n = besta + WAVE_LAG;
      while (hgh >= low)
        if (V[hgh] > n)
          hgh -= 1;                               
        else
          { while (V[low] > n)
              low += 1;
            break;
          }

#ifdef WAVE_STATS
      k = (hgh-low)+1;
      if (k > MAX)
        MAX = k;
      TOT += k;
      NWV += 1;
#endif

#ifdef DEBUG_WAVE
      print_wave(V,M,low,hgh,besta);
#endif
    }

  { uint16 *atrace = (uint16 *) apath->trace;
    uint16 *btrace = (uint16 *) bpath->trace;
    int     atlen, btlen;
    int     trimx;
    int     a, b, k, h;

    if (morem >= 0)
      { trimx  = morea-morey;
        trimy  = morey;
        trimd  = mored;
        trimha = moreha;
        trimhb = morehb;
      }
    else
      trimx = trima-trimy;

    atlen = btlen = 0;

    a = -1;
    for (h = trimha; h >= 0; h = b)
      { b = cells[h].ptr; 
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = cells[h].mark - k;
#ifdef SHOW_TRAIL
    printf("  A path = (%5d,%5d)\n",b+k,b); fflush(stdout);
#endif
    if ((b+k)%TRACE_SPACE != 0)
      { h = cells[h].ptr;
        if (h < 0)
          a = trimy;
        else
          { k = cells[h].diag;
            a = cells[h].mark - k;
          }
#ifdef SHOW_TRAIL
        printf("    +%4d: (%5d,%5d): %3d\n",h,a+k,a,b-a); fflush(stdout);
#endif
        if (apath->tlen == 0)
          atrace[--atlen] = (READIDX) (b-a);
        else
          atrace[0] = (READIDX) (atrace[0] + (b-a));
        b = a;
      }
    if (h >= 0)
      { for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
          { k = cells[h].diag;
            a = cells[h].mark - k;
            atrace[--atlen] = (READIDX) (b-a);
#ifdef SHOW_TRAIL
            printf("     %4d: (%5d,%5d): %3d\n",h,a+k,a,b-a); fflush(stdout);
#endif
            b = a;
          }
        if (b+k != trimx)
          { atrace[--atlen] = (READIDX) (b-trimy);
#ifdef SHOW_TRAIL
            printf("           (%5d,%5d): %3d\n",trimx,trimy,b-trimy); fflush(stdout);
#endif
          }
        else if (b != trimy)
          { atrace[atlen] = (READIDX) (atrace[atlen] + (b-trimy));
#ifdef SHOW_TRAIL
            printf("         @ (%5d,%5d): %3d\n",trimx,trimy,b-trimy); fflush(stdout);
#endif
          }
      }

    a = -1;
    for (h = trimhb; h >= 0; h = b)
      { b = cells[h].ptr; 
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = cells[h].mark + k;
#ifdef SHOW_TRAIL
    printf("  B path = (%5d,%5d)\n",b,b-k); fflush(stdout);
#endif
    if ((b-k)%TRACE_SPACE != boff)
      { h = cells[h].ptr;
        if (h < 0)
          a = trimx;
        else
          { k = cells[h].diag;
            a = cells[h].mark + k;
          }
#ifdef SHOW_TRAIL
        printf("    +%4d: (%5d,%5d): %3d\n",h,a,a-k,b-a); fflush(stdout);
#endif
        if (bpath->tlen == 0)
          btrace[--btlen] = (READIDX) (b-a);
        else
          btrace[0] = (READIDX) (btrace[0] + (b-a));
        b = a;
      }

    if (h >= 0)
      { for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
          { k = cells[h].diag;
            a = cells[h].mark + k;
            btrace[--btlen] = (READIDX) (b-a);
#ifdef SHOW_TRAIL
            printf("     %4d: (%5d,%5d): %3d\n",h,a,a-k,b-a); fflush(stdout);
#endif
            b = a;
          }
        if (b-k != trimy)
          { btrace[--btlen] = (READIDX) (b-trimx);
#ifdef SHOW_TRAIL
            printf("           (%5d,%5d): %3d\n",trimx,trimy,b-trimx); fflush(stdout);
#endif
          }
        else if (b != trimx)
          { btrace[btlen] = (READIDX) (btrace[btlen] + (b-trimx));
#ifdef SHOW_TRAIL
            printf("         @ (%5d,%5d): %3d\n",trimx,trimy,b-trimx); fflush(stdout);
#endif
          }
      }

    apath->abpos = (READIDX) trimx;
    apath->bbpos = (READIDX) trimy;
    apath->diffs = (READIDX) (apath->diffs + trimd);
    apath->tlen  = (READIDX) (apath->tlen  - atlen);
    apath->trace = atrace + atlen;
    if (COMP(align->flags))
      { bpath->aepos = (READIDX) (align->blen - apath->bbpos);
        bpath->bepos = (READIDX) (align->alen - apath->abpos);
      }
    else
      { bpath->abpos = apath->bbpos;
        bpath->bbpos = apath->abpos;
      }
    bpath->diffs = (READIDX) (bpath->diffs + trimd);
    bpath->tlen  = (READIDX) (bpath->tlen  - btlen);
    bpath->trace = btrace + btlen;
  }

  work->cells  = (void *) cells;
  work->celmax = cmax;
}


/* Find the longest local alignment between aseq and bseq through (xcnt,ycnt)
   See associated .h file for the precise definition of the interface.
*/

Path *Local_Alignment(Alignment *align, Work_Data *ework, Align_Spec *espec, int xcnt, int ycnt)
{ _Work_Data  *work = ( _Work_Data *) ework;
  _Align_Spec *spec = (_Align_Spec *) espec;

  Path *apath, *bpath;

  { int alen, blen;
    int maxtp, wsize;

    alen = align->alen;
    blen = align->blen;

    wsize = (6*sizeof(int) + sizeof(BVEC))*(alen+blen+3);
    if (wsize >= work->vecmax)
      enlarge_vector(work,wsize);

    if (alen < blen)
      maxtp = 2*(blen/spec->trace_space+2);
    else
      maxtp = 2*(alen/spec->trace_space+2);
    wsize = 4*maxtp*sizeof(uint16) + sizeof(Path);
    if (wsize > work->pntmax)
      enlarge_points(work,wsize);

    apath = align->path;
    bpath = (Path *) work->points;

    apath->trace = ((uint16 *) (bpath+1)) + maxtp;
    bpath->trace = ((uint16 *) apath->trace) +  2*maxtp;
  }

#ifdef DEBUG_PASSES
  printf("\n");
#endif

  { int l, h, a;

    l = h = xcnt-ycnt;
    a = xcnt+ycnt;

    forward_wave(work,spec,align,bpath,l,h,a);
#ifdef DEBUG_PASSES
    printf("F1 (%d,%d) => (%d,%d) %d\n",xcnt,ycnt,apath->aepos,apath->bepos,apath->diffs);
#endif

    reverse_wave(work,spec,align,bpath,l,h,a);
#ifdef DEBUG_PASSES
    printf("R1 (%d,%d) => (%d,%d) %d\n",xcnt,ycnt,apath->abpos,apath->bbpos,apath->diffs);
#endif
  }

  if (COMP(align->flags))
    { uint16 *trace = (uint16 *) bpath->trace;
      uint16  p;
      int     i, j;

      i = bpath->tlen-1;
      j = 0;
      while (j < i)
        { p = trace[i];
          trace[i] = trace[j];
          trace[j] = p;
          i -= 1;
          j += 1;
        }
    }

#ifdef DEBUG_POINTS
  { uint16 *trace = (uint16 *) apath->trace;
    int     a, h;

    printf("\nA-path (%d,%d)->(%d,%d)",apath->abpos,apath->bbpos,apath->aepos,apath->bepos);
    printf(" %c\n",(COMP(align->flags) ? 'c' : 'n'));
    a = apath->bbpos;
    for (h = 0; h < apath->tlen; h++)
      { int del = trace[h];
        a += del;
        printf("      %d (%d)\n",del,a);
      }
  }

  { uint16 *trace = (uint16 *) bpath->trace;
    int     a, h;

    printf("\nB-path (%d,%d)->(%d,%d)",bpath->abpos,bpath->bbpos,bpath->aepos,bpath->bepos);
    printf(" %c [%d,%d]\n",(COMP(align->flags) ? 'c' : 'n'),align->blen,align->alen);
    a = bpath->bbpos;
    for (h = 0; h < bpath->tlen; h++)
      { int del = trace[h];
        a += del;
        printf("      %d (%d)\n",del,a);
      }
  }
#endif

  return (bpath);
}


/****************************************************************************************\
*                                                                                        *
*  OVERLAP MANIPULATION                                                                  *
*                                                                                        *
\****************************************************************************************/

static int64 PtrSize   = sizeof(void *);
static int64 OvlIOSize = sizeof(Overlap) - sizeof(void *);

int Read_Overlap(FILE *input, Overlap *ovl)
{ if (fread( ((char *) ovl) + PtrSize, OvlIOSize, 1, input) != 1)
    return (1);
  return (0);
}

int Read_Trace(FILE *input, Overlap *ovl, int tbytes)
{ if (tbytes > 0 && ovl->path.tlen > 0)
    { if (fread(ovl->path.trace, tbytes*ovl->path.tlen, 1, input) != 1)
        return (1);
    }
  return (0);
}

void Write_Overlap(FILE *output, Overlap *ovl, int tbytes)
{ fwrite( ((char *) ovl) + PtrSize, OvlIOSize, 1, output);
  if (ovl->path.trace != NULL)
    fwrite(ovl->path.trace,tbytes,ovl->path.tlen,output);
}

void Compress_TraceTo8(Overlap *ovl)
{ uint16 *t16 = (uint16 *) ovl->path.trace;
  uint8  *t8  = (uint8  *) ovl->path.trace;
  int     j;

  for (j = 0; j < ovl->path.tlen; j++)
    t8[j] = (uint8) (t16[j]);
}

void Decompress_TraceTo16(Overlap *ovl)
{ uint16 *t16 = (uint16 *) ovl->path.trace;
  uint8  *t8  = (uint8  *) ovl->path.trace;
  int     j;

  for (j = ovl->path.tlen-1; j >= 0; j--)
    t16[j] = t8[j];
}

void Print_Overlap(FILE *output, Overlap *ovl, int indent)
{ int     i;
  uint16 *trace;

  fprintf(output,"%*s%d vs. ",indent,"",ovl->aread);
  if (COMP(ovl->flags))
    fprintf(output,"c(%d)\n",ovl->bread);
  else
    fprintf(output,"%d\n",ovl->bread);
  fprintf(output,"%*s  [%d,%d] vs [%d,%d] w. %d diffs\n",indent,"",
                 ovl->path.abpos,ovl->path.aepos,ovl->path.bbpos,ovl->path.bepos,ovl->path.diffs);

  trace = (uint16 *) (ovl->path.trace);
  if (trace != NULL)
    { int p = ovl->path.bbpos + trace[0];
      fprintf(output,"%*sTrace: %5d",indent,"",p);
      for (i = 1; i < ovl->path.tlen; i++)
        { if (i%10 == 0)
            fprintf(output,"\n%*s",indent+6,"");
          p += trace[i];
          fprintf(output," %5d",p);
        }
      fprintf(output,"\n");
    }
}

int Check_Trace_Points(Overlap *ovl, int tspace, int verbose, char *fname)
{ int     i, p; 

  if ((ovl->path.aepos-1)/tspace - ovl->path.abpos/tspace != ovl->path.tlen-1)
    { if (verbose) 
        fprintf(stderr,"  %s: Wrong number of trace points\n",fname);
      return (1);
    }         
  p = ovl->path.bbpos;
  if (tspace <= TRACE_XOVR)
    { uint8 *trace8 = (uint8 *) ovl->path.trace;
      for (i = 0; i < ovl->path.tlen; i++)
        p += trace8[i];
    }
  else      
    { uint16 *trace16 = (uint16 *) ovl->path.trace;
      for (i = 0; i < ovl->path.tlen; i++)
        p += trace16[i];
    }
  if (p != ovl->path.bepos)
    { if (verbose)
        fprintf(stderr,"  %s: Trace point sum != aligned interval\n",fname);
      return (1); 
    }         
  return (0);
}


/****************************************************************************************\
*                                                                                        *
*  ALIGNMENT PRINTING                                                                    *
*                                                                                        *
\****************************************************************************************/

static int seqlen(char *seq)
{ int len;

  len = 0;
  while (seq[len] != 4)
    len += 1;
  return (len);
}

/* Complement the sequence in fragment aseq.  The operation does the
   complementation/reversal in place.  Calling it a second time on a
   given fragment restores it to its original state.                */

void Complement_Seq(char *aseq)
{ int   len;
  char *s, *t;
  int   c;

  len = seqlen(aseq);

  s = aseq;
  t = aseq + (len-1);
  while (s < t)
    { c    = 3 - *s;
      *s++ = (char) (3 - *t);
      *t-- = (char) c;
    }
  if (s == t)
    *s = (char) (3 - *s);
}

/* Print an alignment to file between a and b given in trace (unpacked).
   Prefix gives the length of the initial prefix of a that is unaligned.  */

static char ToL[8] = { 'a', 'c', 'g', 't', '.', '[', ']', '-' };
static char ToU[8] = { 'A', 'C', 'G', 'T', '.', '[', ']', '-' };

void Print_Alignment(FILE *file, Alignment *align, Work_Data *ework,
                     int indent, int width, int border, int upper, int coord)
{ _Work_Data *work  = (_Work_Data *) ework;
  int        *trace = align->path->trace;
  int         tlen  = align->path->tlen;

  char *Abuf, *Bbuf, *Dbuf;
  int   i, j, o;
  char *a, *b;
  char  mtag, dtag;
  int   prefa, prefb;
  int   aend, bend;
  int   sa, sb;
  int   match, diff;
  char *N2A;

  if (trace == NULL) return;

#ifdef SHOW_TRACE
  fprintf(file,"\nTrace:\n");
  for (i = 0; i < tlen; i++)
    fprintf(file,"  %3d\n",trace[i]);
#endif

  o = sizeof(char)*3*(width+1);
  if (o > work->vecmax)
    enlarge_vector(work,o);

  if (upper)
    N2A = ToU;
  else
    N2A = ToL;

  Abuf = (char *) work->vector;
  Bbuf = Abuf + (width+1);
  Dbuf = Bbuf + (width+1);

  aend = align->path->aepos;
  bend = align->path->bepos;

  Abuf[width] = Bbuf[width] = Dbuf[width] = '\0';
                                           /* buffer/output next column */
#define COLUMN(x,y)							\
{ int u, v;								\
  if (o >= width)							\
    { fprintf(file,"\n");						\
      fprintf(file,"%*s",indent,"");					\
      if (coord > 0)							\
        { if (sa <= aend)						\
            fprintf(file," %*d",coord,sa);				\
          else								\
            fprintf(file," %*s",coord,"");				\
          fprintf(file," %s\n",Abuf);					\
          fprintf(file,"%*s %*s %s\n",indent,"",coord,"",Dbuf);		\
          fprintf(file,"%*s",indent,"");				\
          if (sb <= bend)						\
            fprintf(file," %*d",coord,sb);				\
          else								\
            fprintf(file," %*s",coord,"");				\
          fprintf(file," %s",Bbuf);					\
        }								\
      else								\
        { fprintf(file," %s\n",Abuf);					\
          fprintf(file,"%*s %s\n",indent,"",Dbuf);			\
          fprintf(file,"%*s %s",indent,"",Bbuf);			\
        }								\
      fprintf(file," %5.1f%%\n",(100.*diff)/(diff+match));		\
      o  = 0;								\
      sa = i;								\
      sb = j;								\
      match = diff = 0;							\
    }									\
  u = (x);								\
  v = (y);								\
  if (u == 4 || v == 4)							\
    Dbuf[o] = ' ';							\
  else if (u == v)							\
    Dbuf[o] = mtag;							\
  else									\
    Dbuf[o] = dtag;							\
  Abuf[o] = N2A[u];							\
  Bbuf[o] = N2A[v];							\
  o += 1;								\
}

  a = align->aseq - 1;
  b = align->bseq - 1;

  o  = 0;
  i = j = 1;

  prefa = align->path->abpos;
  prefb = align->path->bbpos;

  if (prefa > border)
    { i = prefa-(border-1);
      prefa = border;
    }
  if (prefb > border)
    { j = prefb-(border-1);
      prefb = border;
    }

  sa   = i;
  sb   = j;
  mtag = ':';
  dtag = ':';

  while (prefa > prefb)
    { COLUMN(a[i],4)
      i += 1;
      prefa -= 1;
    }
  while (prefb > prefa)
    { COLUMN(4,b[j])
      j += 1;
      prefb -= 1;
    }
  while (prefa > 0)
    { COLUMN(a[i],b[j])
      i += 1;
      j += 1;
      prefa -= 1;
    }

  mtag = '[';
  if (prefb > 0)
    COLUMN(5,5)

  mtag  = '|';
  dtag  = '*';

  match = diff = 0;

  { int p, c;      /* Output columns of alignment til reach trace end */

    for (c = 0; c < tlen; c++)
      if ((p = trace[c]) < 0)
        { p = -p;
          while (i != p)
            { COLUMN(a[i],b[j])
              if (a[i] == b[j])
                match += 1;
              else
                diff += 1;
              i += 1;
              j += 1;
            }
          COLUMN(7,b[j])
          j += 1;
          diff += 1;
        }
      else
        { while (j != p)
            { COLUMN(a[i],b[j])
              if (a[i] == b[j])
                match += 1;
              else
                diff += 1;
              i += 1;
              j += 1;
            }
          COLUMN(a[i],7)
          i += 1;
          diff += 1;
        }
    p = align->path->aepos;
    while (i <= p)
      { COLUMN(a[i],b[j])
        if (a[i] == b[j])
          match += 1;
        else
          diff += 1;
        i += 1;
        j += 1;
      }
  }

  { int c;     /* Output remaining column including unaligned suffix */

    mtag = ']';
    if (a[i] != 4 && b[j] != 4 && border > 0)
      COLUMN(6,6)

    mtag = ':';
    dtag = ':';

    c = 0;
    while (c < border && (a[i] != 4 || b[j] != 4))
      { if (a[i] != 4)
          if (b[j] != 4)
            { COLUMN(a[i],b[j])
              i += 1;
              j += 1;
            }
          else
            { COLUMN(a[i],4)
              i += 1;
            }
        else
          { COLUMN(4,b[j])
            j += 1;
          }
        c += 1;
      }
  }

  /* Print remainder of buffered col.s */

  fprintf(file,"\n");
  fprintf(file,"%*s",indent,"");
  if (coord > 0)
    { if (sa <= aend)
        fprintf(file," %*d",coord,sa);
      else
        fprintf(file," %*s",coord,"");
      fprintf(file," %.*s\n",o,Abuf);
      fprintf(file,"%*s %*s %.*s\n",indent,"",coord,"",o,Dbuf);
      fprintf(file,"%*s",indent,"");
      if (sb <= bend)
        fprintf(file," %*d",coord,sb);
      else
        fprintf(file," %*s",coord,"");
      fprintf(file," %.*s",o,Bbuf);
    }
  else
    { fprintf(file," %.*s\n",o,Abuf);
      fprintf(file,"%*s %.*s\n",indent,"",o,Dbuf);
      fprintf(file,"%*s %.*s",indent,"",o,Bbuf);
    }
  if (diff+match > 0)
    fprintf(file," %5.1f%%\n",(100.*diff)/(diff+match));
  else
    fprintf(file,"\n");

  fflush(file);
}

void Print_Reference(FILE *file, Alignment *align, Work_Data *ework,
                     int indent, int block, int border, int upper, int coord)
{ _Work_Data *work  = (_Work_Data *) ework;
  int        *trace = align->path->trace;
  int         tlen  = align->path->tlen;

  char *Abuf, *Bbuf, *Dbuf;
  int   i, j, o;
  char *a, *b;
  char  mtag, dtag;
  int   prefa, prefb;
  int   aend, bend;
  int   sa, sb, s0;
  int   match, diff;
  char *N2A;

  if (trace == NULL) return;

#ifdef SHOW_TRACE
  fprintf(file,"\nTrace:\n");
  for (i = 0; i < tlen; i++)
    fprintf(file,"  %3d\n",trace[i]);
#endif

  o = sizeof(char)*18*(block+1);
  if (o > work->vecmax)
    enlarge_vector(work,o);

  if (upper)
    N2A = ToU;
  else
    N2A = ToL;

  Abuf = (char *) work->vector;
  Bbuf = Abuf + 6*(block+1);
  Dbuf = Bbuf + 6*(block+1);

  aend = align->path->aepos;
  bend = align->path->bepos;

#define BLOCK(x,y)							\
{ int u, v;								\
  if (i%block == 1 && i != s0 && x != 7 && o > 0)			\
    { fprintf(file,"\n");						\
      fprintf(file,"%*s",indent,"");					\
      if (coord > 0)							\
        { if (sa <= aend)						\
            fprintf(file," %*d",coord,sa);				\
          else								\
            fprintf(file," %*s",coord,"");				\
          fprintf(file," %.*s\n",o,Abuf);				\
          fprintf(file,"%*s %*s %.*s\n",indent,"",coord,"",o,Dbuf);	\
          fprintf(file,"%*s",indent,"");				\
          if (sb <= bend)						\
            fprintf(file," %*d",coord,sb);				\
          else								\
            fprintf(file," %*s",coord,"");				\
          fprintf(file," %.*s",o,Bbuf);					\
        }								\
      else								\
        { fprintf(file," %.*s\n",o,Abuf);				\
          fprintf(file,"%*s %.*s\n",indent,"",o,Dbuf);			\
          fprintf(file,"%*s %.*s",indent,"",o,Bbuf);			\
        }								\
      fprintf(file," %5.1f%%\n",(100.*diff)/(diff+match));		\
      o  = 0;								\
      sa = i;								\
      sb = j;								\
      match = diff = 0;							\
    }									\
  u = (x);								\
  v = (y);								\
  if (u == 4 || v == 4)							\
    Dbuf[o] = ' ';							\
  else if (u == v)							\
    Dbuf[o] = mtag;							\
  else									\
    Dbuf[o] = dtag;							\
  Abuf[o] = N2A[u];							\
  Bbuf[o] = N2A[v];							\
  o += 1;								\
}

  a = align->aseq - 1;
  b = align->bseq - 1;

  o  = 0;
  i = j = 1;

  prefa = align->path->abpos;
  prefb = align->path->bbpos;

  if (prefa > border)
    { i = prefa-(border-1);
      prefa = border;
    }
  if (prefb > border)
    { j = prefb-(border-1);
      prefb = border;
    }

  s0   = i;
  sa   = i;
  sb   = j;
  mtag = ':';
  dtag = ':';

  while (prefa > prefb)
    { BLOCK(a[i],4)
      i += 1;
      prefa -= 1;
    }
  while (prefb > prefa)
    { BLOCK(4,b[j])
      j += 1;
      prefb -= 1;
    }
  while (prefa > 0)
    { BLOCK(a[i],b[j])
      i += 1;
      j += 1;
      prefa -= 1;
    }

  mtag = '[';
  if (prefb > 0)
    BLOCK(5,5)

  mtag  = '|';
  dtag  = '*';

  match = diff = 0;

  { int p, c;      /* Output columns of alignment til reach trace end */

    for (c = 0; c < tlen; c++)
      if ((p = trace[c]) < 0)
        { p = -p;
          while (i != p)
            { BLOCK(a[i],b[j])
              if (a[i] == b[j])
                match += 1;
              else
                diff += 1;
              i += 1;
              j += 1;
            }
          BLOCK(7,b[j])
          j += 1;
          diff += 1;
        }
      else
        { while (j != p)
            { BLOCK(a[i],b[j])
              if (a[i] == b[j])
                match += 1;
              else
                diff += 1;
              i += 1;
              j += 1;
            }
          BLOCK(a[i],7)
          i += 1;
          diff += 1;
        }
    p = align->path->aepos;
    while (i <= p)
      { BLOCK(a[i],b[j])
        if (a[i] == b[j])
		match += 1;
	else
          diff += 1;
        i += 1;
        j += 1;
      }
  }

  { int c;     /* Output remaining column including unaligned suffix */

    mtag = ']';
    if (a[i] != 4 && b[j] != 4 && border > 0)
      BLOCK(6,6)

    mtag = ':';
    dtag = ':';

    c = 0;
    while (c < border && (a[i] != 4 || b[j] != 4))
      { if (a[i] != 4)
          if (b[j] != 4)
            { BLOCK(a[i],b[j])
              i += 1;
              j += 1;
            }
          else
            { BLOCK(a[i],4)
              i += 1;
            }
        else
          { BLOCK(4,b[j])
            j += 1;
          }
        c += 1;
      }
  }

  /* Print remainder of buffered col.s */

  fprintf(file,"\n");
  fprintf(file,"%*s",indent,"");
  if (coord > 0)
    { if (sa <= aend)
        fprintf(file," %*d",coord,sa);
      else
        fprintf(file," %*s",coord,"");
      fprintf(file," %.*s\n",o,Abuf);
      fprintf(file,"%*s %*s %.*s\n",indent,"",coord,"",o,Dbuf);
      fprintf(file,"%*s",indent,"");
      if (sb <= bend)
        fprintf(file," %*d",coord,sb);
      else
        fprintf(file," %*s",coord,"");
      fprintf(file," %.*s",o,Bbuf);
    }
  else
    { fprintf(file," %.*s\n",o,Abuf);
      fprintf(file,"%*s %.*s\n",indent,"",o,Dbuf);
      fprintf(file,"%*s %.*s",indent,"",o,Bbuf);
    }
  if (diff+match > 0)
    fprintf(file," %5.1f%%\n",(100.*diff)/(diff+match));
  else
    fprintf(file,"\n");

  fflush(file);
}

/* Print an ASCII representation of the overlap in align between fragments
   a and b to given file.                                                  */

static void Print_Cartoon(FILE *file, Path *path, int alen, int blen, int comp, int indent)
{
  fprintf(file,"%*s",indent,"");
  if (path->abpos > 0)
    fprintf(file,"   %3d",path->abpos);
  else
    fprintf(file,"      ");
  if (path->aepos < alen)
    fprintf(file,"            %3d",alen-path->aepos);
  fprintf(file,"\n");

  fprintf(file,"%*s",indent,"");
  if (path->abpos > 0)
    fprintf(file,"A =====+--------");
  else
    fprintf(file,"A      ---------");
  if (path->aepos < alen)
    fprintf(file,"+====>");
  else
    fprintf(file,">     ");

  { int asub, bsub;

    asub = path->aepos - path->abpos;
    bsub = path->bepos - path->bbpos;
    fprintf(file,"   dif/(len1+len2) = %d/(%d+%d) = %5.2f%%\n",
                 path->diffs,asub,bsub,(200.*path->diffs)/(asub+bsub));
  }

  { int   sym1e, sym2e;
    int   sym1p, sym2p;

    if (comp > 0)
      { sym1p = '<'; sym2p = '-'; sym1e = '<'; sym2e = '='; }
    else
      { sym1p = '-'; sym2p = '>'; sym1e = '='; sym2e = '>'; }

    fprintf(file,"%*s",indent,"");
    if (path->bbpos > 0)
      fprintf(file,"B %c====+--------",sym1e);
    else
      fprintf(file,"B      %c--------",sym1p);
    if (path->bepos < blen)
      fprintf(file,"+====%c\n",sym2e);
    else
      fprintf(file,"%c\n",sym2p);
  }

  fprintf(file,"%*s",indent,"");
  if (path->bbpos > 0)
    fprintf(file,"   %3d",path->bbpos);
  else
    fprintf(file,"      ");
  if (path->bepos < blen)
    fprintf(file,"            %3d",blen-path->bepos);
  fprintf(file,"\n");

  fflush(file);
} 

void Print_ACartoon(FILE *file, Alignment *align, int indent)
{ int   alen = align->alen;
  int   blen = align->blen;
  Path *path = align->path;
  int   comp = COMP(align->flags);

  Print_Cartoon(file,path,alen,blen,comp,indent);
}

void Print_OCartoon(FILE *file, Overlap *ovl, int indent)
{ int   alen = ovl->alen;
  int   blen = ovl->blen;
  Path *path = &(ovl->path);
  int   comp = COMP(ovl->flags);

  Print_Cartoon(file,path,alen,blen,comp,indent);
}


/****************************************************************************************\
*                                                                                        *
*  O(ND) trace algorithm                                                                 *
*                                                                                        *
\****************************************************************************************/

#ifdef DEBUG_AWAVE

static void print_awave(int *V, int low, int hgh)
{ int k;

  printf("  [%6d,%6d]: ",low,hgh);
  for (k = low; k <= hgh; k++)
    printf(" %3d",V[k]);
  printf("\n");
  fflush(stdout);
}

#endif

#ifdef DEBUG_ALIGN

static int depth = 0;

#endif

typedef struct
  { int  *Stop;          //  Ongoing stack of alignment indels
    char *Aabs, *Babs;   //  Absolute base of A and B sequences

    int  **PVF, **PHF;   //  List of waves for iterative np algorithms
    int   mida,  midb;   //  mid point division for mid-point algorithms

    int   *VF,   *VB;    //  Forward/Reverse waves for nd algorithms
                         //  (defunct: were used for O(nd) algorithms)
  } Trace_Waves;

static int dandc_nd(char *A, int M, char *B, int N, Trace_Waves *wave)
{ int x, y;
  int D;

#ifdef DEBUG_ALIGN
  printf("%*s %ld,%ld: %d vs %d\n",depth,"",A-wave->Aabs,B-wave->Babs,M,N);
#endif

  if (M <= 0)
    { x = (wave->Aabs-A)-1;
      for (y = 1; y <= N; y++)
        { *wave->Stop++ = x;
#ifdef DEBUG_SCRIPT
          printf("%*s *I %ld(%ld)\n",depth,"",y+(B-wave->Babs),(A-wave->Aabs)+1);
#endif
        }
      return (N);
    }
  if (N <= 0)
    { y = (B-wave->Babs)+1;
      for (x = 1; x <= M; x++)
        { *wave->Stop++ = y;
#ifdef DEBUG_SCRIPT
          printf("%*s *D %ld(%ld)\n",depth,"",x+(A-wave->Aabs),(B-wave->Babs)+1);
#endif
        }
      return (M);
    }

  { int  *VF = wave->VF;
    int  *VB = wave->VB;
    int   flow;  //  fhgh == D !
    int   blow, bhgh;
    char *a;

    y = 0;
    if (N < M)
      while (y < N && B[y] == A[y])
        y += 1;
    else
      { while (y < M && B[y] == A[y])
          y += 1;
        if (y >= M && N == M)
          return (0);
      }

    flow   = 0;
    VF[0]  = y;
    VF[-1] = -2;

    x = N-M;
    a = A-x;
    y = N-1;
    if (N > M)
      while (y >= x && B[y] == a[y])
        y -= 1;
    else
      while (y >= 0 && B[y] == a[y])
        y -= 1;

    blow = bhgh = -x;
    VB += x;
    VB[blow]   = y;
    VB[blow-1] = N+1;

    for (D = 1; 1; D += 1)
      { int   k, r;
        int   am, ac, ap;

        //  Forward wave

        flow -= 1;
        am = ac = VF[flow-1] = -2;

        a = A + D;
        x = M - D;
        for (k = D; k >= flow; k--)
          { ap = ac;
            ac = am+1;
            am = VF[k-1];

            if (ac < am)
              if (ap < am)
                y  = am;
              else
                y = ap;
            else
              if (ap < ac)
                y = ac;
              else
                y = ap;

            if (blow <= k && k <= bhgh)
              { r = VB[k];
                if (y > r)
                  { D = (D<<1)-1;
                    if (ap > r)
                      y = ap;
                    else if (ac > r)
                      y = ac;
                    else
                      y = r+1;
                    x = k+y;
                    goto OVERLAP2;
                  }
              }

            if (N < x)
              while (y < N && B[y] == a[y])
                y += 1;
            else
              while (y < x && B[y] == a[y])
                y += 1;
            
            VF[k] = y;
            a -= 1;
            x += 1;
          }

#ifdef DEBUG_AWAVE
        print_awave(VF,flow,D);
#endif

        //  Reverse Wave

        bhgh += 1;
        blow -= 1;
	am = ac = VB[blow-1] = N+1;

        a = A + bhgh;
        x = -bhgh;
        for (k = bhgh; k >= blow; k--)
          { ap = ac+1;
            ac = am;
            am = VB[k-1];

            if (ac > am)
              if (ap > am)
                y  = am;
              else
                y = ap;
            else
              if (ap > ac)
                y = ac;
              else
                y = ap;

            if (flow <= k && k <= D)
              { r = VF[k];
	        if (y <= r)
                  { D = (D << 1);
                    if (ap <= r)
                      y = ap;
                    else if (ac <= r)
                      y = ac;
                    else
                      y = r;
                    x = k+y;
                    goto OVERLAP2;
                  }
              }

            y -= 1;
            if (x > 0)
              while (y >= x && B[y] == a[y])
                y -= 1;
            else
              while (y >= 0 && B[y] == a[y])
                y -= 1;

            VB[k] = y;
            a -= 1;
            x += 1;
          }

#ifdef DEBUG_AWAVE
        print_awave(VB,blow,bhgh);
#endif
      }
  }

OVERLAP2:

#ifdef DEBUG_ALIGN
  printf("%*s (%d,%d) @ %d\n",depth,"",x,y,D);
  fflush(stdout);
#endif
  if (D > 1)
    { 
#ifdef DEBUG_ALIGN
      depth += 2;
#endif
      dandc_nd(A,x,B,y,wave);
      dandc_nd(A+x,M-x,B+y,N-y,wave);
#ifdef DEBUG_ALIGN
      depth -= 2;
#endif
    }
  else if (D == 1)
    { if (M > N)
        { *wave->Stop++ = (B-wave->Babs)+y+1;
#ifdef DEBUG_SCRIPT
          printf("%*s  D %ld(%ld)\n",depth,"",(A-wave->Aabs)+x,(B-wave->Babs)+y+1);
#endif
        }
      else if (M < N)
        { *wave->Stop++ = (wave->Aabs-A)-x-1;
#ifdef DEBUG_SCRIPT
          printf("%*s  I %ld(%ld)\n",depth,"",(B-wave->Babs)+y,(A-wave->Aabs)+x+1);
#endif
        }
#ifdef DEBUG_SCRIPT
      else
        printf("%*s  %ld S %ld\n",depth,"",(wave->Aabs-A)+x,(B-wave->Babs)+y);
#endif
    }

  return (D);
}


static void Compute_Trace_ND_ALL(Alignment *align, Work_Data *ework)
{ _Work_Data *work = (_Work_Data *) ework;
  Trace_Waves wave;

  int   L, D;
  int   asub, bsub;
  Path *path;
  int  *trace;

  path = align->path;
  asub = path->aepos-path->abpos;
  bsub = path->bepos-path->bbpos;

  if (asub < bsub)
    L = bsub;
  else
    L = asub;
  L *= sizeof(int);
  if (L > work->tramax)
    enlarge_trace(work,L);

  trace = wave.Stop = ((int *) work->trace);

  D = 2*(path->diffs + 4)*sizeof(int);
  if (D > work->vecmax)
    enlarge_vector(work,D);
  
  D = (path->diffs+3)/2;
  wave.VF = ((int *) work->vector) + (D+1);
  wave.VB = wave.VF + (2*D+1);

  wave.Aabs = align->aseq;
  wave.Babs = align->bseq;

  path->diffs = (READIDX) dandc_nd(align->aseq+path->abpos,path->aepos-path->abpos,
                                   align->bseq+path->bbpos,path->bepos-path->bbpos,&wave);
  path->trace = trace;
  path->tlen  = wave.Stop - trace;
}


/****************************************************************************************\
*                                                                                        *
*  O(NP) tracing algorithms                                                              *
*                                                                                        *
\****************************************************************************************/

/* Iterative O(np) algorithm for finding the alignment between two substrings (specified
     by a Path record).  The variation includes handling substitutions and guarantees
     to find left-most alignments so that low complexity runs are always aligned in
     the same way.
*/

#ifdef DEBUG_ALIGN

static int ToA[4] = { 'a', 'c', 'g', 't' };

#endif

static int iter_np(char *A, int M, char *B, int N, Trace_Waves *wave)
{ int  **PVF = wave->PVF; 
  int  **PHF = wave->PHF;
  int    D;
  int    del = M-N;

  { int  *F0, *F1, *F2;
    int  *HF;
    int   low, hgh, pos;

#ifdef DEBUG_ALIGN
    printf("\n%*s BASE %ld,%ld: %d vs %d\n",depth,"",A-wave->Aabs,B-wave->Babs,M,N);
    printf("%*s A = ",depth,"");
    for (D = 0; D < M; D++)
      printf("%c",ToA[(int) A[D]]);
    printf("\n");
    printf("%*s B = ",depth,"");
    for (D = 0; D < N; D++)
      printf("%c",ToA[(int) B[D]]);
    printf("\n");
#endif

    if (del >= 0)
      { low = 0;
        hgh = del;
      }
    else
      { low = del;
        hgh = 0;
      }
    if (wave->Aabs == wave->Babs)
      pos = B-A;
    else
      pos = -INT32_MAX;

    F1 = PVF[-2];
    F0 = PVF[-1];

    for (D = low-1; D <= hgh+1; D++)
      F1[D] = F0[D] = -2;
    F0[0] = -1;

    low += 1;
    hgh -= 1;

    for (D = 0; 1; D += 1)
      { int   k, i, j;
        int   am, ac, ap;
        char *a;

        F2 = F1;
        F1 = F0;
        F0 = PVF[D];
        HF = PHF[D];

        if ((D & 0x1) == 0)
          { hgh += 1;
            low -= 1;
            if (low <= pos)
              low += 1;
          }
        F0[hgh+1] = F0[low-1] = -2;

#define FS_MOVE(mdir,pdir)			\
  ac = F1[k]+1;					\
  if (ac < am)					\
    if (ap < am)				\
      { HF[k] = mdir;				\
        j = am;					\
      }						\
    else					\
      { HF[k] = pdir;				\
        j = ap;					\
      }						\
  else						\
    if (ap < ac)				\
      { HF[k] = 0;				\
        j = ac;					\
      }						\
    else					\
      { HF[k] = pdir;				\
        j = ap;					\
      }						\
						\
  if (N < i)					\
    while (j < N && B[j] == a[j])		\
      j += 1;					\
  else						\
    while (j < i && B[j] == a[j])		\
      j += 1;					\
  F0[k] = j;

        j = -2;
        a = A + hgh;
        i = M - hgh;
        for (k = hgh; k > del; k--)
          { ap = j+1;
            am = F2[k-1];
            FS_MOVE(-1,4)
            a -= 1;
            i += 1;
          }

        j = -2;
        a = A + low;
        i = M - low;
        for (k = low; k < del; k++)
          { ap = F2[k+1]+1;
            am = j;
            FS_MOVE(2,1)
            a += 1;
            i -= 1;
          }

        ap = F0[del+1]+1;
        am = j;
        FS_MOVE(2,4)

#ifdef DEBUG_AWAVE
        print_awave(F0,low,hgh);
        print_awave(HF,low,hgh);
#endif

        if (F0[del] >= N)
          break;
      }
  }

  { int   k, h, m, e, c;
    char *a;
    int   ap = (wave->Aabs-A)-1;
    int   bp = (B-wave->Babs)+1;

    PHF[0][0] = 3;

    c = N;
    k = del;
    e = PHF[D][k];
    PHF[D][k] = 3;
    while (e != 3)
      { h = k+e;
        if (e > 1)
          h -= 3;
        else if (e == 0)
          D -= 1;
        else
          D -= 2;
        if (h < k)       // => e = -1 or 2
          { a = A + k;
            if (k < 0)
              m = -k;
            else
              m = 0;
            if (PVF[D][h] <= c)
              c = PVF[D][h]-1;
            while (c >= m && a[c] == B[c])
              c -= 1;
            if (e < 1)  //  => edge is 2, others are 1, and 0
              { if (c <= PVF[D+2][k+1])
                  { e = 4;
                    h = k+1;
                    D = D+2;
                  }
                else if (c == PVF[D+1][k])
                  { e = 0;
                    h = k;
                    D = D+1;
                  }
                else
                  PVF[D][h] = c+1;
              }
            else      //   => edge is 0, others are 1, and 2 (if k != del), 0 (otherwise)
              { if (k == del)
                  m = D;
                else
                  m = D-2;
                if (c <= PVF[m][k+1])
                  { if (k == del)
                      e = 4;
                    else
                      e = 1;
                    h = k+1;
                    D = m;
                  }
                else if (c == PVF[D-1][k])
                  { e = 0;
                    h = k;
                    D = D-1;
                  }
                else
                  PVF[D][h] = c+1;
              }
          }
        m = PHF[D][h];
        PHF[D][h] = e;
        e = m;

        k = h;
      }

    k = D = 0;
    e = PHF[D][k];
    while (e != 3)
      { h = k-e;
        c = PVF[D][k];
        if (e > 1)
          h += 3;
        else if (e == 0)
          D += 1;
        else
          D += 2;
        if (h > k)
          *wave->Stop++ = bp+c;
        else if (h < k)
          *wave->Stop++ = ap-(c+k);
        k = h;
        e = PHF[D][h];
      }

#ifdef DEBUG_SCRIPT
    k = D = 0;
    e = PHF[D][k];
    while (e != 3)
      { h = k-e;
        c = PVF[D][k];
        if (e > 1)
          h += 3;
        else if (e == 0)
          D += 1;
        else
          D += 2;
        if (h > k)
          printf("%*s  D %d(%d)\n",depth,"",(c-k)-(ap-1),c+bp);
        else if (h < k)
          printf("%*s  I %d(%d)\n",depth,"",c+(bp-1),(c+k)-ap);
        else
          printf("%*s  %d S %d\n",depth,"",(c+k)-(ap+1),c+(bp-1));
        k = h;
        e = PHF[D][h];
      }
#endif
  }

  return (D + abs(del));
}

static int middle_np(char *A, int M, char *B, int N, Trace_Waves *wave)
{ int  **PVF = wave->PVF; 
  int  **PHF = wave->PHF;
  int    D;
  int    del = M-N;

  { int  *F0, *F1, *F2;
    int  *HF;
    int   low, hgh, pos;

#ifdef DEBUG_ALIGN
    printf("\n%*s BASE %ld,%ld: %d vs %d\n",depth,"",A-wave->Aabs,B-wave->Babs,M,N);
    printf("%*s A = ",depth,"");
    for (D = 0; D < M; D++)
      printf("%c",ToA[(int) A[D]]);
    printf("\n");
    printf("%*s B = ",depth,"");
    for (D = 0; D < N; D++)
      printf("%c",ToA[(int) B[D]]);
    printf("\n");
#endif

    if (del >= 0)
      { low = 0;
        hgh = del;
      }
    else
      { low = del;
        hgh = 0;
      }
    if (wave->Aabs == wave->Babs)
      pos = B-A;
    else
      pos = -INT32_MAX;

    F1 = PVF[-2];
    F0 = PVF[-1];

    for (D = low-1; D <= hgh+1; D++)
      F1[D] = F0[D] = -2;
    F0[0] = -1;

    low += 1;
    hgh -= 1;

    for (D = 0; 1; D += 1)
      { int   k, i, j;
        int   am, ac, ap;
        char *a;

        F2 = F1;
        F1 = F0;
        F0 = PVF[D];
        HF = PHF[D];

        if ((D & 0x1) == 0)
          { hgh += 1;
            low -= 1;
            if (low <= pos)
              low += 1;
          }
        F0[hgh+1] = F0[low-1] = -2;

        j = -2;
        a = A + hgh;
        i = M - hgh;
        for (k = hgh; k > del; k--)
          { ap = j+1;
            am = F2[k-1];
            FS_MOVE(-1,4)
            a -= 1;
            i += 1;
          }

        j = -2;
        a = A + low;
        i = M - low;
        for (k = low; k < del; k++)
          { ap = F2[k+1]+1;
            am = j;
            FS_MOVE(2,1)
            a += 1;
            i -= 1;
          }

        ap = F0[del+1]+1;
        am = j;
        FS_MOVE(2,4)

#ifdef DEBUG_AWAVE
        print_awave(F0,low,hgh);
        print_awave(HF,low,hgh);
#endif

        if (F0[del] >= N)
          break;
      }
  }

  { int   k, h, m, e, c;
    int   d, f;
    char *a;

    d = D + abs(del);
    c = N;
    k = del;
    for (f = d/2; d > f; d--)
      { e = PHF[D][k];
        h = k+e;
        if (e > 1)
          h -= 3;
        else if (e == 0)
          D -= 1;
        else
          D -= 2;
        if (h < k)       // => e = -1 or 2
          { a = A + k;
            if (k < 0)
              m = -k;
            else
              m = 0;
            if (PVF[D][h] <= c)
              c = PVF[D][h]-1;
            while (c >= m && a[c] == B[c])
              c -= 1;
            if (e < 1)  //  => edge is 2, others are 1, and 0
              { if (c <= PVF[D+2][k+1])
                  { e = 4;
                    h = k+1;
                    D = D+2;
                  }
                else if (c == PVF[D+1][k])
                  { e = 0;
                    h = k;
                    D = D+1;
                  }
                else
                  PVF[D][h] = c+1;
              }
            else      //   => edge is 0, others are 1, and 2 (if k != del), 0 (otherwise)
              { if (k == del)
                  m = D;
                else
                  m = D-2;
                if (c <= PVF[m][k+1])
                  { if (k == del)
                      e = 4;
                    else
                      e = 1;
                    h = k+1;
                    D = m;
                  }
                else if (c == PVF[D-1][k])
                  { e = 0;
                    h = k;
                    D = D-1;
                  }
                else
                  PVF[D][h] = c+1;
              }
          }
        k = h;
      }

    wave->midb = (B-wave->Babs) + PVF[D][k];
    wave->mida = (A-wave->Aabs) + k + PVF[D][k];
  }

  return (1);
}


/****************************************************************************************\
*                                                                                        *
*  COMPUTE_TRACE FLAVORS                                                                 *
*                                                                                        *
\****************************************************************************************/

void Compute_Trace_ALL(Alignment *align, Work_Data *ework)
{ _Work_Data *work = (_Work_Data *) ework;
  Trace_Waves wave;

  Path *path;
  char *aseq, *bseq;
  int   M, N;

  path = align->path;
  aseq = align->aseq;
  bseq = align->bseq;

  M = path->aepos-path->abpos;
  N = path->bepos-path->bbpos;
  
  { int64 s;
    int   d;
    int   dmax;
    int   **PVF, **PHF;

    if (M < N)
      s = N;
    else
      s = M;
    s *= sizeof(int);
    if (s > work->tramax)
      enlarge_trace(work,s);

    dmax = path->diffs - abs(M-N);

    s = (dmax+3)*2*((M+N+3)*sizeof(int) + sizeof(int *));

    if (s > 256000000)
      { Compute_Trace_ND_ALL(align,ework);
        return;
      }

    if (s > work->vecmax)
      enlarge_vector(work,s);

    wave.PVF = PVF = ((int **) (work->vector)) + 2;
    wave.PHF = PHF = PVF + (dmax+3);

    s = M+N+3;
    PVF[-2] = ((int *) (PHF + (dmax+1))) + (N+1);
    for (d = -1; d <= dmax; d++)
      PVF[d] = PVF[d-1] + s;
    PHF[-2] = PVF[dmax] + s;
    for (d = -1; d <= dmax; d++)
      PHF[d] = PHF[d-1] + s;
  }

  wave.Stop = ((int *) work->trace);
  wave.Aabs = aseq;
  wave.Babs = bseq;

  path->diffs = (READIDX) iter_np(aseq+path->abpos,M,bseq+path->bbpos,N,&wave);
  path->trace = work->trace;
  path->tlen  = wave.Stop - ((int *) path->trace);
}

void Compute_Trace_PTS(Alignment *align, Work_Data *ework, int trace_spacing)
{ _Work_Data *work = (_Work_Data *) ework;
  Trace_Waves wave;

  Path   *path;
  char   *aseq, *bseq;
  uint16 *points;
  int     tlen;
  int     ab, bb;
  int     ae, be;
  int     diffs;

  path   = align->path;
  aseq   = align->aseq;
  bseq   = align->bseq;
  tlen   = path->tlen;
  points = (uint16 *) path->trace;

  { int64 s;
    int   d;
    int   M, N;
    int   dmax, nmax;
    int   **PVF, **PHF;

    M = path->aepos-path->abpos;
    N = path->bepos-path->bbpos;
    if (M < N)
      s = N*sizeof(int);
    else
      s = M*sizeof(int);
    if (s > work->tramax)
      enlarge_trace(work,s);

    nmax = 0;
    for (d = 0; d < tlen; d++)
      { if (points[d] > nmax)
          nmax = points[d];
      }
    if (tlen <= 1)
      nmax = N;
    dmax = nmax;

    s = (dmax+3)*2*((trace_spacing+nmax+3)*sizeof(int) + sizeof(int *));

    if (s > work->vecmax)
      enlarge_vector(work,s);

    wave.PVF = PVF = ((int **) (work->vector)) + 2;
    wave.PHF = PHF = PVF + (dmax+3);

    s = trace_spacing+nmax+3;
    PVF[-2] = ((int *) (PHF + (dmax+1))) + (nmax+1);
    for (d = -1; d <= dmax; d++)
      PVF[d] = PVF[d-1] + s;
    PHF[-2] = PVF[dmax] + s;
    for (d = -1; d <= dmax; d++)
      PHF[d] = PHF[d-1] + s;
  }

  wave.Stop = (int *) (work->trace);
  wave.Aabs = aseq;
  wave.Babs = bseq;

  { int i;

    diffs = 0;
    ab = path->abpos;
    ae = (ab/trace_spacing)*trace_spacing;
    bb = path->bbpos;
    tlen -= 1;
    for (i = 0; i < tlen; i++)
      { ae = ae + trace_spacing;
        be = bb + points[i];
        diffs += iter_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave);
        ab = ae;
        bb = be;
      }
    ae = path->aepos;
    be = path->bepos;
    diffs += iter_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave);
  }

  path->trace = work->trace;
  path->tlen  = wave.Stop - ((int *) path->trace);
  path->diffs = (READIDX) diffs;
}

void Compute_Trace_MID(Alignment *align, Work_Data *ework, int trace_spacing)
{ _Work_Data *work = (_Work_Data *) ework;
  Trace_Waves wave;

  Path   *path;
  char   *aseq, *bseq;
  uint16 *points;
  int     tlen;
  int     ab, bb;
  int     ae, be;
  int     diffs;

  path   = align->path;
  aseq   = align->aseq;
  bseq   = align->bseq;
  tlen   = path->tlen;
  points = (uint16 *) path->trace;

  { int64 s;
    int   d;
    int   M, N;
    int   dmax, nmax;
    int   **PVF, **PHF;

    M = path->aepos-path->abpos;
    N = path->bepos-path->bbpos;
    if (M < N)
      s = N*sizeof(int);
    else
      s = M*sizeof(int);
    if (s > work->tramax)
      enlarge_trace(work,s);

    nmax = 0;
    for (d = 0; d < tlen; d++)
      { if (points[d] > nmax)
          nmax = points[d];
      }
    if (tlen <= 1)
      nmax = N;
    dmax = nmax;

    s = (dmax+3)*4*((trace_spacing+nmax+3)*sizeof(int) + sizeof(int *));

    if (s > work->vecmax)
      enlarge_vector(work,s);

    wave.PVF = PVF = ((int **) (work->vector)) + 2;
    wave.PHF = PHF = PVF + (dmax+3);

    s = trace_spacing+nmax+3;
    PVF[-2] = ((int *) (PHF + (dmax+1))) + (nmax+1);
    for (d = -1; d <= dmax; d++)
      PVF[d] = PVF[d-1] + s;
    PHF[-2] = PVF[dmax] + s;
    for (d = -1; d <= dmax; d++)
      PHF[d] = PHF[d-1] + s;
  }

  wave.Stop = ((int *) work->trace);
  wave.Aabs = aseq;
  wave.Babs = bseq;

  { int i;
    int as, bs;
    int af, bf;

    diffs = 0;
    ab = as = af = path->abpos;
    ae = (ab/trace_spacing)*trace_spacing;
    bb = bs = bf = path->bbpos;
    tlen -= 1;
    for (i = 0; i < tlen; i++) 
      { ae = ae + trace_spacing;
        be = bb + points[i];
        if (middle_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave))
          { af = wave.mida;
            bf = wave.midb;
            diffs += iter_np(aseq+as,af-as,bseq+bs,bf-bs,&wave);
            ab = ae;
            bb = be;
            as = af;
            bs = bf;
          }
      }
    ae = path->aepos;
    be = path->bepos;
    if (middle_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave))
      { af = wave.mida;
        bf = wave.midb;
        diffs  += iter_np(aseq+as,af-as,bseq+bs,bf-bs,&wave);
        as = af;
        bs = bf;
      }
    diffs += iter_np(aseq+af,ae-as,bseq+bf,be-bs,&wave);
  }

  path->trace = work->trace;
  path->tlen  = wave.Stop - ((int *) path->trace);
  path->diffs = (READIDX) diffs;
}
