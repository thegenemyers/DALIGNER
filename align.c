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
    EXIT(NULL);
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

static int enlarge_vector(_Work_Data *work, int newmax)
{ void *vec;
  int   max;

  max = ((int) (newmax*1.2)) + 10000;
  vec = Realloc(work->vector,max,"Enlarging DP vector");
  if (vec == NULL)
    EXIT(1);
  work->vecmax = max;
  work->vector = vec;
  return (0);
}

static int enlarge_points(_Work_Data *work, int newmax)
{ void *vec;
  int   max;

  max = ((int) (newmax*1.2)) + 10000;
  vec = Realloc(work->points,max,"Enlarging point vector");
  if (vec == NULL)
    EXIT(1);
  work->pntmax = max;
  work->points = vec;
  return (0);
}

static int enlarge_trace(_Work_Data *work, int newmax)
{ void *vec;
  int   max;

  max = ((int) (newmax*1.2)) + 10000;
  vec = Realloc(work->trace,max,"Enlarging trace vector");
  if (vec == NULL)
    EXIT(1);
  work->tramax = max;
  work->trace  = vec;
  return (0);
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
    EXIT(NULL);

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
    { EPRINTF(EPLACE,"Base bias worse than 80/20%% ! (New_Align_Spec)\n");
      free(spec);
      EXIT(NULL);
    }

  spec->ave_path = (int) (PATH_LEN * (1. - Bias_Factor[bias] * (1. - ave_corr)));
  parms.mscore   = (int) (FRACTION * Bias_Factor[bias] * (1. - ave_corr));
  parms.dscore   = FRACTION - parms.mscore;

  parms.score = (int16 *) Malloc(sizeof(int16)*(TRIM_MASK+1)*2,"Allocating trim table");
  if (parms.score == NULL)
    { free(spec);
      EXIT(NULL);
    }
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
    int diff;
    int mark;
  } Pebble;

static int VectorEl = 6*sizeof(int) + sizeof(BVEC);

static int forward_wave(_Work_Data *work, _Align_Spec *spec, Alignment *align, Path *bpath,
                        int *mind, int maxd, int mida, int minp, int maxp)
{ char *aseq  = align->aseq;
  char *bseq  = align->bseq;
  Path *apath = align->path;

  int     hgh, low, dif;
  int     vlen, vmin, vmax;
  int    *V, *M;
  int    *_V, *_M;
  BVEC   *T;
  BVEC   *_T;

  int    *HA, *HB;
  int    *_HA, *_HB;
  int    *NA, *NB;
  int    *_NA, *_NB;
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

  hgh = maxd;
  low = *mind;
  dif = 0;

  { int span, wing;

    span = (hgh-low)+1;
    vlen = work->vecmax/VectorEl;
    wing = (vlen - span)/2;
    vmin = low - wing;
    vmax = hgh + wing;

    _V  = ((int *) work->vector);
    _M  = _V + vlen;
    _HA = _M + vlen;
    _HB = _HA + vlen;
    _NA = _HB + vlen;
    _NB = _NA + vlen;
    _T  = ((BVEC *) (_NB + vlen));

    V  = _V-vmin;
    M  = _M-vmin;
    HA = _HA-vmin;
    HB = _HB-vmin;
    NA = _NA-vmin;
    NB = _NB-vmin;
    T  = _T-vmin;

    cells = (Pebble *) (work->cells);
    cmax  = work->celmax;
    avail = 0;

    if (COMP(align->flags))
      boff = align->blen % TRACE_SPACE;
    else
      boff = 0;
  }

  /* Compute 0-wave starting from mid-line */

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
              EXIT(1);
            work->celmax = cmax;
            work->cells  = (void *) cells;
          }

        na = ((y+k)/TRACE_SPACE)*TRACE_SPACE;
#ifdef SHOW_TPS
        printf(" A %d: %d,%d,0,%d\n",avail,-1,k,na); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
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
        pb->diff = 0;
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
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" A %d: %d,%d,0,%d\n",avail,ha,k,na); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = ha;
            pb->diag = k;
            pb->diff = 0;
            pb->mark = na;
            ha  = avail++;
            na += TRACE_SPACE;
          }
        while (y >= nb)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" B %d: %d,%d,0,%d\n",avail,hb,k,nb); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = hb;
            pb->diag = k;
            pb->diff = 0;
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

      if (low <= vmin || hgh >= vmax)
        { int   span, wing;
          int64 move;
          int64 vd, md, had, hbd, nad, nbd, td;

          span = (hgh-low)+1;
          if (.8*vlen < span)
            { if (enlarge_vector(work,vlen*VectorEl))
                EXIT(1);

              move = ((void *) _V) - work->vector;
              vlen = work->vecmax/VectorEl;

              _V  = (int *) work->vector;
              _M  = _V + vlen;
              _HA = _M + vlen;
              _HB = _HA + vlen;
              _NA = _HB + vlen;
              _NB = _NA + vlen;
              _T  = ((BVEC *) (_NB + vlen));
            }
          else
            move = 0;

          wing = (vlen - span)/2;

          vd  = ((void *) ( _V+wing)) - (((void *) ( V+low)) - move);
          md  = ((void *) ( _M+wing)) - (((void *) ( M+low)) - move);
          had = ((void *) (_HA+wing)) - (((void *) (HA+low)) - move);
          hbd = ((void *) (_HB+wing)) - (((void *) (HB+low)) - move);
          nad = ((void *) (_NA+wing)) - (((void *) (NA+low)) - move);
          nbd = ((void *) (_NB+wing)) - (((void *) (NB+low)) - move);
          td  = ((void *) ( _T+wing)) - (((void *) ( T+low)) - move);

          if (vd < 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));
          if (md < 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (had < 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (hbd < 0)
            memmove(_HB+wing,  ((void *) (HB+low)) - move, span*sizeof(int));
          if (nad < 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (nbd < 0)
            memmove(_NB+wing,  ((void *) (NB+low)) - move, span*sizeof(int));
          if (td < 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));

          if (td > 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));
          if (nbd > 0)
            memmove(_NB+wing,  ((void *) (NB+low)) - move, span*sizeof(int));
          if (nad > 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (hbd > 0)
            memmove(_HB+wing,  ((void *) (HB+low)) - move, span*sizeof(int));
          if (had > 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (md > 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (vd > 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));

          vmin = low-wing;
          vmax = hgh+wing;

          V  =  _V-vmin;
          M  =  _M-vmin;
          HA = _HA-vmin;
          HB = _HB-vmin;
          NA = _NA-vmin;
          NB = _NB-vmin;
          T  =  _T-vmin;
        }

      if (low > minp)
        { low -= 1;
          NA[low] = NA[low+1];
          NB[low] = NB[low+1];
          V[low]  = -1;
        }
      if (hgh < maxp)
        { hgh += 1;
          NA[hgh] = NA[hgh-1];
          NB[hgh] = NB[hgh-1];
          V[hgh]  = am = -1;
        }
      else
        am = V[hgh];
      dif += 1;

      ac = V[hgh+1] = V[low-1] = -1;
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
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" A %d: %d,%d,%d,%d\n",avail,ha,k,dif,NA[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = ha;
                  pb->diag = k;
                  pb->diff = dif;
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
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" B %d: %d,%d,%d,%d\n",avail,hb,k,dif,NB[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = hb;
                  pb->diag = k;
                  pb->diff = dif;
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
    int     d, e;

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
    e = 0;
#ifdef SHOW_TRAIL
    printf("  A path = (%5d,%5d)\n",(mida+k)/2,b); fflush(stdout);
#endif
    for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
      { k = cells[h].diag;
        a = cells[h].mark - k;
        d = cells[h].diff;
        atrace[atlen++] = (uint16) (d-e);
        atrace[atlen++] = (uint16) (a-b);
#ifdef SHOW_TRAIL
        printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,a-b); fflush(stdout);
#endif
        b = a;
        e = d;
      }
    if (b+k != trimx)
      { atrace[atlen++] = (uint16) (trimd-e);
        atrace[atlen++] = (uint16) (trimy-b);
#ifdef SHOW_TRAIL
        printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimy-b); fflush(stdout);
#endif
      }
    else if (b != trimy)
      { atrace[atlen-1] = (uint16) (atrace[atlen-1] + (trimy-b));
        atrace[atlen-2] = (uint16) (atrace[atlen-2] + (trimd-e));
#ifdef SHOW_TRAIL
        printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimy-b); fflush(stdout);
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
    e = 0;
    low = k;
#ifdef SHOW_TRAIL
    printf("  B path = (%5d,%5d)\n",b,(mida-k)/2); fflush(stdout);
#endif
    for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
      { k = cells[h].diag;
        a = cells[h].mark + k;
        d = cells[h].diff;
        btrace[btlen++] = (uint16) (d-e);
        btrace[btlen++] = (uint16) (a-b);  
#ifdef SHOW_TRAIL
        printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a,a-k,d-e,a-b); fflush(stdout);
#endif
        b = a;
        e = d;
      }
    if (b-k != trimy)
      { btrace[btlen++] = (uint16) (trimd-e);
        btrace[btlen++] = (uint16) (trimx-b);  
#ifdef SHOW_TRAIL
        printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimx-b); fflush(stdout);
#endif
      }
    else if (b != trimx)
      { btrace[btlen-1] = (uint16) (btrace[btlen-1] + (trimx-b));
        btrace[btlen-2] = (uint16) (btrace[btlen-2] + (trimd-e));
#ifdef SHOW_TRAIL
        printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimx-b); fflush(stdout);
#endif
      }

    apath->aepos = trimx;
    apath->bepos = trimy;
    apath->diffs = trimd;
    apath->tlen  = atlen;
    if (COMP(align->flags))
      { bpath->abpos = align->blen - apath->bepos;
        bpath->bbpos = align->alen - apath->aepos;
      }
    else
      { bpath->aepos = apath->bepos;
        bpath->bepos = apath->aepos;
      }
    bpath->diffs = trimd;
    bpath->tlen  = btlen;
  }

  *mind = low;
  return (0);
}

/*** Reverse Wave ***/

static int reverse_wave(_Work_Data *work, _Align_Spec *spec, Alignment *align, Path *bpath,
                        int mind, int maxd, int mida, int minp, int maxp)
{ char *aseq  = align->aseq - 1;
  char *bseq  = align->bseq - 1;
  Path *apath = align->path;

  int     hgh, low, dif;
  int     vlen, vmin, vmax;
  int    *V, *M;
  int    *_V, *_M;
  BVEC   *T;
  BVEC   *_T;

  int    *HA, *HB;
  int    *_HA, *_HB;
  int    *NA, *NB;
  int    *_NA, *_NB;
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

  hgh = maxd;
  low = mind;
  dif = 0;

  { int span, wing;

    span = (hgh-low)+1;
    vlen = work->vecmax/VectorEl;
    wing = (vlen - span)/2;
    vmin = low - wing;
    vmax = hgh + wing;

    _V  = ((int *) work->vector);
    _M  = _V + vlen;
    _HA = _M + vlen;
    _HB = _HA + vlen;
    _NA = _HB + vlen;
    _NB = _NA + vlen;
    _T  = ((BVEC *) (_NB + vlen));

    V  = _V-vmin;
    M  = _M-vmin;
    HA = _HA-vmin;
    HB = _HB-vmin;
    NA = _NA-vmin;
    NB = _NB-vmin;
    T  = _T-vmin;

    cells = (Pebble *) (work->cells);
    cmax  = work->celmax;
    avail = 0;

    if (COMP(align->flags))
      boff = align->blen % TRACE_SPACE;
    else
      boff = 0;
  }

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
              EXIT(1);
            work->celmax = cmax;
            work->cells  = (void *) cells;
          }

        na = ((y+k+TRACE_SPACE-1)/TRACE_SPACE-1)*TRACE_SPACE;
#ifdef SHOW_TPS
        printf(" A %d: -1,%d,0,%d\n",avail,k,na+TRACE_SPACE); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
        pb->mark = y+k;
        ha  = avail++;

        nb = ((y+(TRACE_SPACE-boff)-1)/TRACE_SPACE-1)*TRACE_SPACE+boff;
#ifdef SHOW_TPS
        printf(" B %d: -1,%d,0,%d\n",avail,k,nb+TRACE_SPACE); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
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
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" A %d: %d,%d,0,%d\n",avail,ha,k,na); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = ha;
            pb->diag = k;
            pb->diff = 0;
            pb->mark = na;
            ha  = avail++;
            na -= TRACE_SPACE;
          }
        while (y <= nb)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" B %d: %d,%d,0,%d\n",avail,hb,k,nb); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = hb;
            pb->diag = k;
            pb->diff = 0;
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

      if (low <= vmin || hgh >= vmax)
        { int   span, wing;
          int64 move, vd, md, had, hbd, nad, nbd, td;

          span = (hgh-low)+1;
          if (.8*vlen < span)
            { if (enlarge_vector(work,vlen*VectorEl))
                EXIT(1);

              move = ((void *) _V) - work->vector;
              vlen = work->vecmax/VectorEl;

              _V  = (int *) work->vector;
              _M  = _V + vlen;
              _HA = _M + vlen;
              _HB = _HA + vlen;
              _NA = _HB + vlen;
              _NB = _NA + vlen;
              _T  = ((BVEC *) (_NB + vlen));
            }
          else
            move = 0;

          wing = (vlen - span)/2;

          vd  = ((void *) ( _V+wing)) - (((void *) ( V+low)) - move);
          md  = ((void *) ( _M+wing)) - (((void *) ( M+low)) - move);
          had = ((void *) (_HA+wing)) - (((void *) (HA+low)) - move);
          hbd = ((void *) (_HB+wing)) - (((void *) (HB+low)) - move);
          nad = ((void *) (_NA+wing)) - (((void *) (NA+low)) - move);
          nbd = ((void *) (_NB+wing)) - (((void *) (NB+low)) - move);
          td  = ((void *) ( _T+wing)) - (((void *) ( T+low)) - move);

          if (vd < 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));
          if (md < 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (had < 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (hbd < 0)
            memmove(_HB+wing,  ((void *) (HB+low)) - move, span*sizeof(int));
          if (nad < 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (nbd < 0)
            memmove(_NB+wing,  ((void *) (NB+low)) - move, span*sizeof(int));
          if (td < 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));

          if (td > 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));
          if (nbd > 0)
            memmove(_NB+wing,  ((void *) (NB+low)) - move, span*sizeof(int));
          if (nad > 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (hbd > 0)
            memmove(_HB+wing,  ((void *) (HB+low)) - move, span*sizeof(int));
          if (had > 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (md > 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (vd > 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));

          vmin = low-wing;
          vmax = hgh+wing;

          V  =  _V-vmin;
          M  =  _M-vmin;
          HA = _HA-vmin;
          HB = _HB-vmin;
          NA = _NA-vmin;
          NB = _NB-vmin;
          T  =  _T-vmin;
        }

      if (low > minp)
        { low -= 1;
          NA[low] = NA[low+1];
          NB[low] = NB[low+1];
          V[low]  = ap = INT32_MAX;
        }
      else
        ap = V[low]; 
      if (hgh < maxp)
        { hgh += 1;
          NA[hgh] = NA[hgh-1];
          NB[hgh] = NB[hgh-1];
          V[hgh] = INT32_MAX;
        }
      dif += 1;

      ac = V[hgh+1] = V[low-1] = INT32_MAX;
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
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" A %d: %d,%d,%d,%d\n",avail,ha,k,dif,NA[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = ha;
                  pb->diag = k;
                  pb->diff = dif;
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
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" B %d: %d,%d,%d,%d\n",avail,hb,k,dif,NB[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = hb;
                  pb->diag = k;
                  pb->diff = dif;
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
    int     d, e;

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
    e = 0;
#ifdef SHOW_TRAIL
    printf("  A path = (%5d,%5d)\n",b+k,b); fflush(stdout);
#endif
    if ((b+k)%TRACE_SPACE != 0)
      { h = cells[h].ptr;
        if (h < 0)
          { a = trimy;
            d = trimd;
          }
        else
          { k = cells[h].diag;
            a = cells[h].mark - k;
            d = cells[h].diff;
          }
#ifdef SHOW_TRAIL
        printf("    +%4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,b-a); fflush(stdout);
#endif
        if (apath->tlen == 0)
          { atrace[--atlen] = (uint16) (b-a);
            atrace[--atlen] = (uint16) (d-e);
          }
        else
          { atrace[1] = (uint16) (atrace[1] + (b-a));
            atrace[0] = (uint16) (atrace[0] + (d-e));
          }
        b = a;
        e = d;
      }
    if (h >= 0)
      { for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
          { k = cells[h].diag;
            a = cells[h].mark - k;
            atrace[--atlen] = (uint16) (b-a);
            d = cells[h].diff;
            atrace[--atlen] = (uint16) (d-e);
#ifdef SHOW_TRAIL
            printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,b-a); fflush(stdout);
#endif
            b = a;
            e = d;
          }
        if (b+k != trimx)
          { atrace[--atlen] = (uint16) (b-trimy);
            atrace[--atlen] = (uint16) (trimd-e);
#ifdef SHOW_TRAIL
            printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimy); fflush(stdout);
#endif
          }
        else if (b != trimy)
          { atrace[atlen+1] = (uint16) (atrace[atlen+1] + (b-trimy));
            atrace[atlen]   = (uint16) (atrace[atlen]   + (trimd-e));
#ifdef SHOW_TRAIL
            printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimy); fflush(stdout);
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
    e = 0;
#ifdef SHOW_TRAIL
    printf("  B path = (%5d,%5d)\n",b,b-k); fflush(stdout);
#endif
    if ((b-k)%TRACE_SPACE != boff)
      { h = cells[h].ptr;
        if (h < 0)
          { a = trimx;
            d = trimd;
          } 
        else
          { k = cells[h].diag;
            a = cells[h].mark + k;
            d = cells[h].diff;
          }
#ifdef SHOW_TRAIL
        printf("    +%4d: (%5d,%5d): %3d / %3d\n",h,a,a-k,d-e,b-a); fflush(stdout);
#endif
        if (bpath->tlen == 0)
          { btrace[--btlen] = (uint16) (b-a);
            btrace[--btlen] = (uint16) (b-a);
          }
        else
          { btrace[1] = (uint16) (btrace[1] + (b-a));
            btrace[0] = (uint16) (btrace[0] + (d-e));
          }
        b = a;
        e = d;
      }

    if (h >= 0)
      { for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
          { k = cells[h].diag;
            a = cells[h].mark + k;
            btrace[--btlen] = (uint16) (b-a);
            d = cells[h].diff;
            btrace[--btlen] = (uint16) (d-e);
#ifdef SHOW_TRAIL
            printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a,a-k,d-e,b-a); fflush(stdout);
#endif
            b = a;
            e = d;
          }
        if (b-k != trimy)
          { btrace[--btlen] = (uint16) (b-trimx);
            btrace[--btlen] = (uint16) (trimd-e);
#ifdef SHOW_TRAIL
            printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimx); fflush(stdout);
#endif
          }
        else if (b != trimx)
          { btrace[btlen+1] = (uint16) (btrace[btlen+1] + (b-trimx));
            btrace[btlen]   = (uint16) (btrace[btlen]   + (trimd-e));
#ifdef SHOW_TRAIL
            printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimx); fflush(stdout);
#endif
          }
      }

    apath->abpos = trimx;
    apath->bbpos = trimy;
    apath->diffs = apath->diffs + trimd;
    apath->tlen  = apath->tlen  - atlen;
    apath->trace = atrace + atlen;
    if (COMP(align->flags))
      { bpath->aepos = align->blen - apath->bbpos;
        bpath->bepos = align->alen - apath->abpos;
      }
    else
      { bpath->abpos = apath->bbpos;
        bpath->bbpos = apath->abpos;
      }
    bpath->diffs = bpath->diffs + trimd;
    bpath->tlen  = bpath->tlen  - btlen;
    bpath->trace = btrace + btlen;
  }

  return (0);
}


/* Find the longest local alignment between aseq and bseq through (xcnt,ycnt)
   See associated .h file for the precise definition of the interface.
*/

Path *Local_Alignment(Alignment *align, Work_Data *ework, Align_Spec *espec,
                      int low, int hgh, int anti, int lbord, int hbord)
{ _Work_Data  *work = ( _Work_Data *) ework;
  _Align_Spec *spec = (_Align_Spec *) espec;

  Path *apath, *bpath;
  int   minp, maxp;
  int   selfie;

  { int alen, blen;
    int maxtp, wsize;

    alen = align->alen;
    blen = align->blen;

    if (hgh-low >= 7500)
      wsize = VectorEl*(hgh-low+1);
    else
      wsize = VectorEl*10000;
    if (wsize >= work->vecmax)
      if (enlarge_vector(work,wsize))
        EXIT(NULL);

    if (alen < blen)
      maxtp = 2*(blen/spec->trace_space+2);
    else
      maxtp = 2*(alen/spec->trace_space+2);
    wsize = 4*maxtp*sizeof(uint16) + sizeof(Path);
    if (wsize > work->pntmax)
      if (enlarge_points(work,wsize))
        EXIT(NULL);

    apath = align->path;
    bpath = (Path *) work->points;

    apath->trace = ((uint16 *) (bpath+1)) + maxtp;
    bpath->trace = ((uint16 *) apath->trace) +  2*maxtp;
  }

#ifdef DEBUG_PASSES
  printf("\n");
#endif

  selfie = (align->aseq == align->bseq);
   
  if (lbord < 0)
    { if (selfie && low >= 0)
        minp = 1;
      else
        minp = -INT32_MAX;
    }
  else
    minp = low-lbord;
  if (hbord < 0)
    { if (selfie && hgh <= 0)
        maxp = -1;
      else
        maxp = INT32_MAX;
    }
  else
    maxp = hgh+hbord;

  if (forward_wave(work,spec,align,bpath,&low,hgh,anti,minp,maxp))
    EXIT(NULL);

#ifdef DEBUG_PASSES
  printf("F1 (%d,%d) ~ %d => (%d,%d) %d\n",
         (2*anti+(low+hgh))/4,(anti-(low+hgh))/4,hgh-low,
         apath->aepos,apath->bepos,apath->diffs);
#endif

  if (reverse_wave(work,spec,align,bpath,low,low,anti,minp,maxp))
    EXIT(NULL);

#ifdef DEBUG_PASSES
  printf("R1 (%d,%d) => (%d,%d) %d\n",
         (anti+low)/2,(anti-low)/2,apath->abpos,apath->bbpos,apath->diffs);
#endif

  if (COMP(align->flags))
    { uint16 *trace = (uint16 *) bpath->trace;
      uint16  p;
      int     i, j;

      i = bpath->tlen-2;
      j = 0;
      while (j < i)
        { p = trace[i];
          trace[i] = trace[j];
          trace[j] = p;
          p = trace[i+1];
          trace[i+1] = trace[j+1];
          trace[j+1] = p;
          i -= 2;
          j += 2;
        }
    }

#ifdef DEBUG_POINTS
  { uint16 *trace = (uint16 *) apath->trace;
    int     a, h;

    printf("\nA-path (%d,%d)->(%d,%d)",apath->abpos,apath->bbpos,apath->aepos,apath->bepos);
    printf(" %c\n",(COMP(align->flags) ? 'c' : 'n'));
    a = apath->bbpos;
    for (h = 1; h < apath->tlen; h += 2)
      { int dif = trace[h-1];
        int del = trace[h];
        a += del;
        printf("      %d / %d (%d)\n",dif,del,a);
      }
  }

  { uint16 *trace = (uint16 *) bpath->trace;
    int     a, h;

    printf("\nB-path (%d,%d)->(%d,%d)",bpath->abpos,bpath->bbpos,bpath->aepos,bpath->bepos);
    printf(" %c [%d,%d]\n",(COMP(align->flags) ? 'c' : 'n'),align->blen,align->alen);
    a = bpath->bbpos;
    for (h = 1; h < bpath->tlen; h += 2)
      { int dif = trace[h-1];
        int del = trace[h];
        a += del;
        printf("      %d / %d (%d)\n",dif,del,a);
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

void Print_Overlap(FILE *output, Overlap *ovl, int tbytes, int indent)
{ int     i;

  fprintf(output,"%*s%d vs. ",indent,"",ovl->aread);
  if (COMP(ovl->flags))
    fprintf(output,"c(%d)\n",ovl->bread);
  else
    fprintf(output,"%d\n",ovl->bread);
  fprintf(output,"%*s  [%d,%d] vs [%d,%d] w. %d diffs\n",indent,"",
                 ovl->path.abpos,ovl->path.aepos,ovl->path.bbpos,ovl->path.bepos,ovl->path.diffs);

  if (tbytes == 1)
    { uint8 *trace = (uint8 *) (ovl->path.trace);
      if (trace != NULL)
        { int p = ovl->path.bbpos + trace[1];
          fprintf(output,"%*sTrace: %3d/%5d",indent,"",trace[0],p);
          for (i = 3; i < ovl->path.tlen; i += 2)
            { if (i%10 == 0)
                fprintf(output,"\n%*s",indent+6,"");
              p += trace[i];
              fprintf(output," %3d/%5d",trace[i-1],p);
            }
          fprintf(output,"\n");
        }
    }
  else
    { uint16 *trace = (uint16 *) (ovl->path.trace);
      if (trace != NULL)
        { int p = ovl->path.bbpos + trace[1];
          fprintf(output,"%*sTrace: %3d/%5d",indent,"",trace[0],p);
          for (i = 3; i < ovl->path.tlen; i += 2)
            { if (i%10 == 0)
                fprintf(output,"\n%*s",indent+6,"");
              p += trace[i];
              fprintf(output," %3d/%5d",trace[i-1],p);
            }
          fprintf(output,"\n");
        }
    }
}

int Check_Trace_Points(Overlap *ovl, int tspace, int verbose, char *fname)
{ int     i, p; 

  if (((ovl->path.aepos-1)/tspace - ovl->path.abpos/tspace)*2 != ovl->path.tlen-2)
    { if (verbose) 
        EPRINTF(EPLACE,"  %s: Wrong number of trace points\n",fname);
      return (1);
    }         
  p = ovl->path.bbpos;
  if (tspace <= TRACE_XOVR)
    { uint8 *trace8 = (uint8 *) ovl->path.trace;
      for (i = 1; i < ovl->path.tlen; i += 2)
        p += trace8[i];
    }
  else      
    { uint16 *trace16 = (uint16 *) ovl->path.trace;
      for (i = 1; i < ovl->path.tlen; i += 2)
        p += trace16[i];
    }
  if (p != ovl->path.bepos)
    { if (verbose)
        EPRINTF(EPLACE,"  %s: Trace point sum != aligned interval\n",fname);
      return (1); 
    }         
  return (0);
}


void Flip_Alignment(Alignment *align, int full)
{ char *aseq  = align->aseq;
  char *bseq  = align->bseq;
  int   alen  = align->alen;
  int   blen  = align->blen;
  Path *path  = align->path;
  int   comp  = COMP(align->flags);

  int  *trace = (int *) path->trace;
  int   tlen  = path->tlen;

  int   i, j, p;

  if (comp)
    { p = path->abpos;
      path->abpos = blen - path->bepos;
      path->bepos = alen - p;
      p = path->aepos;
      path->aepos = blen - path->bbpos;
      path->bbpos = alen - p;

      if (full)
        { alen += 2;
          blen += 2;

          for (i = 0; i < tlen; i++)
            if ((p = trace[i]) < 0)
              trace[i] = alen + p;
            else
              trace[i] = p - blen;

          i = tlen-1;
          j = 0;
          while (j < i)
            { p = trace[i];
              trace[i] = trace[j];
              trace[j] = p;
              i -= 1;
              j += 1;
            }

          alen -= 2;
          blen -= 2;
        }
    }
  else
    { p = path->abpos;
      path->abpos = path->bbpos;
      path->bbpos = p;
      p = path->aepos;
      path->aepos = path->bepos;
      path->bepos = p;

      if (full)
        for (i = 0; i < tlen; i++)
          trace[i] = - (trace[i]);
    }

  align->aseq  = bseq;
  align->bseq  = aseq;
  align->alen  = blen;
  align->blen  = alen;
}


/****************************************************************************************\
*                                                                                        *
*  ALIGNMENT PRINTING                                                                    *
*                                                                                        *
\****************************************************************************************/

/* Complement the sequence in fragment aseq.  The operation does the
   complementation/reversal in place.  Calling it a second time on a
   given fragment restores it to its original state.                */

void Complement_Seq(char *aseq, int len)
{ char *s, *t;
  int   c;

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

int Print_Alignment(FILE *file, Alignment *align, Work_Data *ework,
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

  if (trace == NULL) return (0);

#ifdef SHOW_TRACE
  fprintf(file,"\nTrace:\n");
  for (i = 0; i < tlen; i++)
    fprintf(file,"  %3d\n",trace[i]);
#endif

  o = sizeof(char)*3*(width+1);
  if (o > work->vecmax)
    if (enlarge_vector(work,o))
      EXIT(1);

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
  return (0);
}

int Print_Reference(FILE *file, Alignment *align, Work_Data *ework,
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
  int   vmax;

  if (trace == NULL) return (0);

#ifdef SHOW_TRACE
  fprintf(file,"\nTrace:\n");
  for (i = 0; i < tlen; i++)
    fprintf(file,"  %3d\n",trace[i]);
#endif

  vmax = work->vecmax/3;
  o = sizeof(char)*6*(block+1);
  if (o > vmax)
    { if (enlarge_vector(work,3*o))
        EXIT(1);
      vmax = work->vecmax/3;
    }

  Abuf = (char *) work->vector;
  Bbuf = Abuf + vmax;
  Dbuf = Bbuf + vmax;

  if (upper)
    N2A = ToU;
  else
    N2A = ToL;

  aend = align->path->aepos;
  bend = align->path->bepos;

#define BLOCK(x,y)							\
{ int u, v;								\
  if (i%block == 1 && i != s0 && x < 4 && o > 0)			\
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
  if (o >= vmax)							\
    { if (enlarge_vector(work,3*o))					\
        EXIT(1);							\
      vmax = work->vecmax/3;						\
      memmove(work->vector+2*vmax,Dbuf,o);				\
      memmove(work->vector+vmax,Bbuf,o);				\
      memmove(work->vector,Abuf,o);					\
      Abuf = (char *) work->vector;					\
      Bbuf = Abuf + vmax;						\
      Dbuf = Bbuf + vmax;						\
    }									\
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
  return (0);
}

/* Print an ASCII representation of the overlap in align between fragments
   a and b to given file.                                                  */

static inline void repchar(FILE *file, int symbol, int rep)
{ while (rep-- > 0)
    fputc(symbol,file);
}

void Alignment_Cartoon(FILE *file, Alignment *align, int indent, int coord)
{ int   alen = align->alen;
  int   blen = align->blen;
  Path *path = align->path;
  int   comp = COMP(align->flags);
  int   w;

  fprintf(file,"%*s",indent,"");
  if (path->abpos > 0)
    fprintf(file,"    %*d ",coord,path->abpos);
  else
    fprintf(file,"%*s",coord+5,"");
  if (path->aepos < alen)
    fprintf(file,"%*s%d",coord+8,"",alen-path->aepos);
  fprintf(file,"\n");

  fprintf(file,"%*s",indent,"");
  if (path->abpos > 0)
    { fprintf(file,"A ");
      w = Number_Digits((int64) path->abpos);
      repchar(file,' ',coord-w);
      repchar(file,'=',w+3);
      fputc('+',file);
      repchar(file,'-',coord+5);
    }
  else
    { fprintf(file,"A %*s",coord+4,"");
      repchar(file,'-',coord+5);
    }

  if (path->aepos < alen)
    { fputc('+',file);
      w = Number_Digits((int64) (alen-path->aepos));
      repchar(file,'=',w+2);
      fputc('>',file);
      repchar(file,' ',w);
    }
  else
    { fputc('>',file);
      repchar(file,' ',coord+3);
    }

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
      { fprintf(file,"B ");
        w = Number_Digits((int64) path->bbpos);
        repchar(file,' ',coord-w);
        fputc(sym1e,file);
        repchar(file,'=',w+2);
        fputc('+',file);
        repchar(file,'-',coord+5);
      }
    else
      { fprintf(file,"B ");
        repchar(file,' ',coord+3);
        fputc(sym1p,file);
        repchar(file,'-',coord+5);
      }
    if (path->bepos < blen)
      { fprintf(file,"+");
        w = Number_Digits((int64) (blen-path->bepos));
        repchar(file,'=',w+2);
        fprintf(file,"%c\n",sym2e);
      }
    else
      fprintf(file,"%c\n",sym2p);
  }

  fprintf(file,"%*s",indent,"");
  if (path->bbpos > 0)
    fprintf(file,"    %*d ",coord,path->bbpos);
  else
    fprintf(file,"%*s",coord+5,"");
  if (path->bepos < blen)
    fprintf(file,"%*s%d",coord+8,"",blen-path->bepos);
  fprintf(file,"\n");

  fflush(file);
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


static int Compute_Trace_ND_ALL(Alignment *align, Work_Data *ework)
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
    if (enlarge_trace(work,L))
      EXIT(1);

  trace = wave.Stop = ((int *) work->trace);

  D = 2*(path->diffs + 4)*sizeof(int);
  if (D > work->vecmax)
    if (enlarge_vector(work,D))
      EXIT(1);
  
  D = (path->diffs+3)/2;
  wave.VF = ((int *) work->vector) + (D+1);
  wave.VB = wave.VF + (2*D+1);

  wave.Aabs = align->aseq;
  wave.Babs = align->bseq;

  path->diffs = dandc_nd(align->aseq+path->abpos,path->aepos-path->abpos,
                         align->bseq+path->bbpos,path->bepos-path->bbpos,&wave);
  path->trace = trace;
  path->tlen  = wave.Stop - trace;
  return (0);
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
    int   low, hgh;
    int   posl, posh;

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

    posl = -INT32_MAX;
    posh =  INT32_MAX;
    if (wave->Aabs == wave->Babs)
      { if (B == A)
          { EPRINTF(EPLACE,"Error: self comparison starts on diagonal 0 (Compute_Trace)\n");
            EXIT(-1);
          }
        else if (B < A)
          posl = (B-A)+1;
        else
          posh = (B-A)-1;
      }

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
          { if (low > posl)
              low -= 1;
            if (hgh < posh)
              hgh += 1;
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
    int   low, hgh;
    int   posl, posh;

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

    posl = -INT32_MAX;
    posh =  INT32_MAX;
    if (wave->Aabs == wave->Babs)
      { if (B == A)
          { EPRINTF(EPLACE,"Error: self comparison starts on diagonal 0 (Compute_Trace)\n");
            EXIT(1);
          }
        else if (B < A)
          posl = (B-A)+1;
        else
          posh = (B-A)-1;
      }

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
          { if (low > posl)
              low -= 1;
            if (hgh < posh)
              hgh += 1;
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

  return (0);
}


/****************************************************************************************\
*                                                                                        *
*  COMPUTE_TRACE FLAVORS                                                                 *
*                                                                                        *
\****************************************************************************************/

int Compute_Trace_ALL(Alignment *align, Work_Data *ework)
{ _Work_Data *work = (_Work_Data *) ework;
  Trace_Waves wave;

  Path *path;
  char *aseq, *bseq;
  int   M, N, D;

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
      if (enlarge_trace(work,s))
        EXIT(1);

    dmax = path->diffs - abs(M-N);

    s = (dmax+3)*2*((M+N+3)*sizeof(int) + sizeof(int *));

    if (s > 256000000)
      return (Compute_Trace_ND_ALL(align,ework));

    if (s > work->vecmax)
      if (enlarge_vector(work,s))
        EXIT(1);

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

  D = iter_np(aseq+path->abpos,M,bseq+path->bbpos,N,&wave);
  if (D < 0)
    EXIT(1);
  path->diffs = D;
  path->trace = work->trace;
  path->tlen  = wave.Stop - ((int *) path->trace);

  return (0);
}

int Compute_Trace_PTS(Alignment *align, Work_Data *ework, int trace_spacing)
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
      if (enlarge_trace(work,s))
        EXIT(1);

    nmax = 0;
    dmax = 0;
    for (d = 1; d < tlen; d += 2)
      { if (points[d-1] > dmax)
          dmax = points[d-1];
        if (points[d] > nmax)
          nmax = points[d];
      }
    if (tlen <= 1)
      nmax = N;

    s = (dmax+3)*2*((trace_spacing+nmax+3)*sizeof(int) + sizeof(int *));

    if (s > work->vecmax)
      if (enlarge_vector(work,s))
        EXIT(1);

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

  { int i, d;

    diffs = 0;
    ab = path->abpos;
    ae = (ab/trace_spacing)*trace_spacing;
    bb = path->bbpos;
    tlen -= 2;
    for (i = 1; i < tlen; i += 2)
      { ae = ae + trace_spacing;
        be = bb + points[i];
        d  = iter_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave);
        if (d < 0)
          EXIT(1);
        diffs += d;
        ab = ae;
        bb = be;
      }
    ae = path->aepos;
    be = path->bepos;
    d  = iter_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave);
    if (d < 0)
      EXIT(1);
    diffs += d;
  }

  path->trace = work->trace;
  path->tlen  = wave.Stop - ((int *) path->trace);
  path->diffs = diffs;

  return (0);
}

int Compute_Trace_MID(Alignment *align, Work_Data *ework, int trace_spacing)
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
      if (enlarge_trace(work,s))
        EXIT(1);

    nmax = 0;
    dmax = 0;
    for (d = 1; d < tlen; d += 2)
      { if (points[d-1] > dmax)
          dmax = points[d-1];
        if (points[d] > nmax)
          nmax = points[d];
      }
    if (tlen <= 1)
      nmax = N;

    s = (dmax+3)*4*((trace_spacing+nmax+3)*sizeof(int) + sizeof(int *));

    if (s > work->vecmax)
      if (enlarge_vector(work,s))
        EXIT(1);

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

  { int i, d;
    int as, bs;
    int af, bf;

    diffs = 0;
    ab = as = af = path->abpos;
    ae = (ab/trace_spacing)*trace_spacing;
    bb = bs = bf = path->bbpos;
    tlen -= 2;
    for (i = 1; i < tlen; i += 2) 
      { ae = ae + trace_spacing;
        be = bb + points[i];
        if (middle_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave))
          EXIT(1);
        af = wave.mida;
        bf = wave.midb;
        d  = iter_np(aseq+as,af-as,bseq+bs,bf-bs,&wave);
        if (d < 0)
          EXIT(1);
        diffs += d;
        ab = ae;
        bb = be;
        as = af;
        bs = bf;
      }

    ae = path->aepos;
    be = path->bepos;

    if (middle_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave))
      EXIT(1);
    af = wave.mida;
    bf = wave.midb;
    d  = iter_np(aseq+as,af-as,bseq+bs,bf-bs,&wave);
    if (d < 0)
      EXIT(1);
    diffs += d;
    as = af;
    bs = bf;
    
    d += iter_np(aseq+af,ae-as,bseq+bf,be-bs,&wave);
    if (d < 0)
      EXIT(1);
    diffs += d;
  }

  path->trace = work->trace;
  path->tlen  = wave.Stop - ((int *) path->trace);
  path->diffs = diffs;

  return (0);
}
