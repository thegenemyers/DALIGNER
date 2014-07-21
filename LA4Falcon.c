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
 *  Utility for displaying the overlaps in a .las file in a variety of ways including
 *    a minimal listing of intervals, a cartoon, and a full out alignment.
 *
 *  Author:    Gene Myers
 *  Creation:  July 2013
 *  Last Mod:  Jan 2015
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"
#include "align.h"

static char *Usage[] =
    { "[-mfsocarUFM] [-i<int(4)>] [-w<int(100)>] [-b<int(10)>] ",
      "    <src1:db|dam> [ <src2:db|dam> ] <align:las> [ <reads:FILE> | <reads:range> ... ]"
    };

#define LAST_READ_SYMBOL  '$'

static int ORDER(const void *l, const void *r)
{ int x = *((int32 *) l);
  int y = *((int32 *) r);
  return (x-y);
}

int main(int argc, char *argv[])
{ HITS_DB   _db1, *db1 = &_db1;
  HITS_DB   _db2, *db2 = &_db2;
  Overlap   _ovl, *ovl = &_ovl;
  Alignment _aln, *aln = &_aln;

  FILE   *input;
  int64   novl;
  int     tspace, tbytes, small;
  int     reps, *pts;
  int     input_pts;

  int     ALIGN, CARTOON, REFERENCE, FLIP;
  int     INDENT, WIDTH, BORDER, UPPERCASE;
  int     ISTWO;
  int     MAP;
  int     FALCON, OVERLAP, M4OVL;
  int     SEED_MIN, SKIP;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("LA4Falcon")

    INDENT    = 4;
    WIDTH     = 100;
    BORDER    = 10;

    FALCON    = 0;
    M4OVL     = 0;
    SEED_MIN  = 8000;
    SKIP      = 0;

    ALIGN     = 0;
    REFERENCE = 0;
    CARTOON   = 0;
    FLIP      = 0;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("smfocarUFM")
            break;
          case 'i':
            ARG_NON_NEGATIVE(INDENT,"Indent")
            break;
          case 'w':
            ARG_POSITIVE(WIDTH,"Alignment width")
            break;
          case 'b':
            ARG_NON_NEGATIVE(BORDER,"Alignment border")
            break;
          case 'H':
            ARG_POSITIVE(SEED_MIN,"seed threshold (in bp)")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    UPPERCASE = flags['U'];
    ALIGN     = flags['a'];
    REFERENCE = flags['r'];
    CARTOON   = flags['c'];
    FLIP      = flags['F'];
    MAP       = flags['M'];
    OVERLAP   = flags['o'];
    M4OVL     = flags['m'];
    FALCON    = flags['f'];
    SKIP      = flags['s'];

    if (argc <= 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        exit (1);
      }
  }

  //  Open trimmed DB or DB pair

  { int   status;
    char *pwd, *root;
    FILE *input;

    ISTWO  = 0;
    status = Open_DB(argv[1],db1);
    if (status < 0)
      exit (1);
    if (db1->part > 0)
      { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[1]);
        exit (1);
      }

    if (argc > 3)
      { pwd   = PathTo(argv[3]);
        root  = Root(argv[3],".las");
        if ((input = fopen(Catenate(pwd,"/",root,".las"),"r")) != NULL)
          { ISTWO = 1;
            fclose(input);
            status = Open_DB(argv[2],db2);
            if (status < 0)
              exit (1);
            if (db2->part > 0)
              { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[2]);
                exit (1);
              }
            Trim_DB(db2);
          }
        else
          db2 = db1;
        free(root);
        free(pwd);
      }
    else
      db2 = db1;
    Trim_DB(db1);
  }

  //  Process read index arguments into a sorted list of read ranges

  input_pts = 0;
  if (argc == ISTWO+4)
    { if (argv[ISTWO+3][0] != LAST_READ_SYMBOL || argv[ISTWO+3][1] != '\0')
        { char *eptr, *fptr;
          int   b, e;

          b = strtol(argv[ISTWO+3],&eptr,10);
          if (eptr > argv[ISTWO+3] && b > 0)
            { if (*eptr == '-')
                { if (eptr[1] != LAST_READ_SYMBOL || eptr[2] != '\0')
                    { e = strtol(eptr+1,&fptr,10);
                      input_pts = (fptr <= eptr+1 || *fptr != '\0' || e <= 0);
                    }
                }
              else
                input_pts = (*eptr != '\0');
            }
          else
            input_pts = 1;
        }
    }

  if (input_pts)
    { int v, x;
      FILE *input;

      input = Fopen(argv[ISTWO+3],"r");
      if (input == NULL)
        exit (1);

      reps = 0;
      while ((v = fscanf(input," %d",&x)) != EOF)
        if (v == 0)
          { fprintf(stderr,"%s: %d'th item of input file %s is not an integer\n",
                           Prog_Name,reps+1,argv[2]);
            exit (1);
          }
        else
          reps += 1;

      reps *= 2;
      pts   = (int *) Malloc(sizeof(int)*reps,"Allocating read parameters");
      if (pts == NULL)
        exit (1);

      rewind(input);
      for (v = 0; v < reps; v += 2)
        { fscanf(input," %d",&x);
          pts[v] = pts[v+1] = x;
        }

      fclose(input);
    }

  else
    { pts  = (int *) Malloc(sizeof(int)*2*argc,"Allocating read parameters");
      if (pts == NULL)
        exit (1);

      reps = 0;
      if (argc > 3+ISTWO)
        { int   c, b, e;
          char *eptr, *fptr;

          for (c = 3+ISTWO; c < argc; c++)
            { if (argv[c][0] == LAST_READ_SYMBOL)
                { b = db1->nreads;
                  eptr = argv[c]+1;
                }
              else
                b = strtol(argv[c],&eptr,10);
              if (eptr > argv[c])
                { if (b <= 0)
                    { fprintf(stderr,"%s: %d is not a valid index\n",Prog_Name,b);
                      exit (1);
                    }
                  if (*eptr == '\0')
                    { pts[reps++] = b;
                      pts[reps++] = b;
                      continue;
                    }
                  else if (*eptr == '-')
                    { if (eptr[1] == LAST_READ_SYMBOL)
                        { e = INT32_MAX;
                          fptr = eptr+2;
                        }
                      else
                        e = strtol(eptr+1,&fptr,10);
                      if (fptr > eptr+1 && *fptr == 0 && e > 0)
                        { pts[reps++] = b;
                          pts[reps++] = e;
                          if (b > e)
                            { fprintf(stderr,"%s: Empty range '%s'\n",Prog_Name,argv[c]);
                              exit (1);
                            }
                          continue;
                        }
                    }
                }
              fprintf(stderr,"%s: argument '%s' is not an integer range\n",Prog_Name,argv[c]);
              exit (1);
            }

          qsort(pts,reps/2,sizeof(int64),ORDER);

          b = 0;
          for (c = 0; c < reps; c += 2)
            if (b > 0 && pts[b-1] >= pts[c]-1)
              { if (pts[c+1] > pts[b-1])
                  pts[b-1] = pts[c+1];
              }
            else
              { pts[b++] = pts[c];
                pts[b++] = pts[c+1];
              }
          pts[b++] = INT32_MAX;
          reps = b;
        }
      else
        { pts[reps++] = 1;
          pts[reps++] = INT32_MAX;
        }
    }

  //  Initiate file reading and read (novl, tspace) header

  { char  *over, *pwd, *root;

    pwd   = PathTo(argv[2+ISTWO]);
    root  = Root(argv[2+ISTWO],".las");
    over  = Catenate(pwd,"/",root,".las");
    input = Fopen(over,"r");
    if (input == NULL)
      exit (1);

    if (fread(&novl,sizeof(int64),1,input) != 1)
      SYSTEM_ERROR
    if (fread(&tspace,sizeof(int),1,input) != 1)
      SYSTEM_ERROR

    if (tspace <= TRACE_XOVR)
      { small  = 1;
        tbytes = sizeof(uint8);
      }
    else
      { small  = 0;
        tbytes = sizeof(uint16);
      }

    if (!(FALCON || M4OVL))
      {
        printf("\n%s: ",root);
        Print_Number(novl,0,stdout);
        printf(" records\n");
      }

    free(pwd);
    free(root);
  }

  //  Read the file and display selected records

  { int        j;
    uint16    *trace;
    Work_Data *work;
    int        tmax;
    int        in, npt, idx, ar;
    int64      tps;
    int64      p_aread = -1;
    char       buffer[131072];
    int        skip_rest = 0;

    char      *abuffer, *bbuffer;
    int        ar_wide, br_wide;
    int        ai_wide, bi_wide;
    int        mn_wide, mx_wide;
    int        tp_wide;
    int        blast, match, seen, lhalf, rhalf;

    aln->path = &(ovl->path);
    if (ALIGN || REFERENCE || FALCON)
      { work = New_Work_Data();
        abuffer = New_Read_Buffer(db1);
        bbuffer = New_Read_Buffer(db2);
      }
    else
      { abuffer = NULL;
        bbuffer = NULL;
        work = NULL;
      }

    tmax  = 1000;
    trace = (uint16 *) Malloc(sizeof(uint16)*tmax,"Allocating trace vector");
    if (trace == NULL)
      exit (1);

    in  = 0;
    npt = pts[0];
    idx = 1;

    ar_wide = Number_Digits((int64) db1->nreads);
    br_wide = Number_Digits((int64) db2->nreads);
    ai_wide = Number_Digits((int64) db1->maxlen);
    bi_wide = Number_Digits((int64) db2->maxlen);
    if (db1->maxlen < db2->maxlen)
      { mn_wide = ai_wide;
        mx_wide = bi_wide;
        tp_wide = Number_Digits((int64) db1->maxlen/tspace+2);
      }
    else
      { mn_wide = bi_wide;
        mx_wide = ai_wide;
        tp_wide = Number_Digits((int64) db2->maxlen/tspace+2);
      }
    ar_wide += (ar_wide-1)/3;
    br_wide += (br_wide-1)/3;
    ai_wide += (ai_wide-1)/3;
    bi_wide += (bi_wide-1)/3;
    mn_wide += (mn_wide-1)/3;
    tp_wide += (tp_wide-1)/3;

    if (FLIP)
      { int x;
        x = ar_wide; ar_wide = br_wide; br_wide = x;
        x = ai_wide; ai_wide = bi_wide; bi_wide = x;
      }

    //  For each record do

    blast = -1;
    match = 0;
    seen  = 0;
    lhalf = rhalf = 0;
    for (j = 0; j < novl; j++)

       //  Read it in

      { Read_Overlap(input,ovl);
        if (ovl->path.tlen > tmax)
          { tmax = ((int) 1.2*ovl->path.tlen) + 100;
            trace = (uint16 *) Realloc(trace,sizeof(uint16)*tmax,"Allocating trace vector");
            if (trace == NULL)
              exit (1);
          }
        ovl->path.trace = (void *) trace;
        Read_Trace(input,ovl,tbytes);

        //  Determine if it should be displayed

        ar = ovl->aread+1;
        if (in)
          { while (ar > npt)
              { npt = pts[idx++];
                if (ar < npt)
                  { in = 0;
                    break;
                  }
                npt = pts[idx++];
              }
          }
        else
          { while (ar >= npt)
              { npt = pts[idx++];
                if (ar <= npt)
                  { in = 1;
                    break;
                  }
                npt = pts[idx++];
              }
          }
        if (!in)
          continue;

        //  Display it

        aln->alen  = db1->reads[ovl->aread].rlen;
        aln->blen  = db2->reads[ovl->bread].rlen;
        aln->flags = ovl->flags;
        tps        = ((ovl->path.aepos-1)/tspace - ovl->path.abpos/tspace);

        if (OVERLAP && !FALCON)
          { if (ovl->path.abpos != 0 && ovl->path.bbpos != 0)
              continue;
            if (ovl->path.aepos != aln->alen && ovl->path.bepos != aln->blen)
              continue;
          }

        if (MAP)
          { while (ovl->bread != blast)
              { if (!match && seen && !(lhalf && rhalf))
                  { printf("Missing ");
                    Print_Number((int64) blast+1,br_wide+1,stdout);
                    printf(" %d ->%lld\n",db2->reads[blast].rlen,db2->reads[blast].coff);
                  }
                match = 0;
                seen  = 0;
                lhalf = rhalf = 0;
                blast += 1;
              }
            seen = 1;
            if (ovl->path.abpos == 0)
              rhalf = 1;
            if (ovl->path.aepos == aln->alen)
              lhalf = 1;
            if (ovl->path.bbpos != 0 || ovl->path.bepos != aln->blen)
              continue;
            match = 1;
          }

        // printf(" %7d %7d\n",ovl->path.abpos,ovl->path.aepos);
        // continue;

        if (!(FALCON || M4OVL) ) {
            if (ALIGN || CARTOON || REFERENCE)
              printf("\n");
            if (FLIP)
              { Flip_Alignment(aln,0);
                Print_Number((int64) ovl->bread+1,ar_wide+1,stdout);
                printf("  ");
                Print_Number((int64) ovl->aread+1,br_wide+1,stdout);
              }
            else
              { Print_Number((int64) ovl->aread+1,ar_wide+1,stdout);
                printf("  ");
                Print_Number((int64) ovl->bread+1,br_wide+1,stdout);
              }
            if (COMP(ovl->flags))
              printf(" c");
            else
              printf(" n");
            printf("   [");
            Print_Number((int64) ovl->path.abpos,ai_wide,stdout);
            printf("..");
            Print_Number((int64) ovl->path.aepos,ai_wide,stdout);
            printf("] x [");
            Print_Number((int64) ovl->path.bbpos,bi_wide,stdout);
            printf("..");
            Print_Number((int64) ovl->path.bepos,bi_wide,stdout);
            printf("]");
        }

        //  Display it
        if (M4OVL)
          {
            int64 bbpos, bepos;
            double acc;



            if (COMP(ovl->flags)) {
                bbpos = (int64) aln->blen - (int64) ovl->path.bepos;
                bepos = (int64) aln->blen - (int64) ovl->path.bbpos;
            } else {
                bbpos = (int64) ovl->path.bbpos;
                bepos = (int64) ovl->path.bepos;
            }
            acc = 100-(200. * ovl->path.diffs)/( ovl->path.aepos - ovl->path.abpos + ovl->path.bepos - ovl->path.bbpos
);
            printf("%09lld %09lld %lld %5.2f ", (int64) ovl->aread, (int64) ovl->bread,  (int64) bbpos - (int64) bepos
, acc);
            printf("0 %lld %lld %lld ", (int64) ovl->path.abpos, (int64) ovl->path.aepos, (int64) aln->alen);
            printf("%d %lld %lld %lld ", COMP(ovl->flags), bbpos, bepos, (int64) aln->blen);
            if ( ((int64) aln->blen < (int64) aln->alen) && ((int64) ovl->path.bbpos < 1) && ((int64) aln->blen - (int64) ovl->path.bepos < 1) )
              {
                printf("contains\n");
              }
            else if ( ((int64) aln->alen < (int64) aln->blen) && ((int64) ovl->path.abpos < 1) && ((int64) aln->alen -
 (int64) ovl->path.aepos < 1) )
              {
                printf("contained\n");
              }
            else
              {
                printf("overlap\n");
              }

          }
        if (FALCON)
          {
            if (p_aread == -1) {
                Load_Read(db1, ovl->aread, abuffer, 2);
                printf("%08d %s\n", ovl->aread, abuffer);
                p_aread = ovl->aread;
                skip_rest = 0;
            }
            if (p_aread != ovl -> aread ) {
                printf("+ +\n");
                Load_Read(db1, ovl->aread,abuffer, 2);
                printf("%08d %s\n", ovl->aread, abuffer);
                p_aread = ovl->aread;
                skip_rest = 0;
            }

            if (skip_rest == 0) {
                Load_Read(db2, ovl->bread, bbuffer, 0);
                p_aread = ovl->aread;
                if (COMP(aln->flags))
                  Complement_Seq(bbuffer, aln->blen);
                Upper_Read(bbuffer);
                strncpy( buffer, bbuffer + ovl->path.bbpos, (int64) ovl->path.bepos - (int64) ovl->path.bbpos );
                buffer[ (int64) ovl->path.bepos - (int64) ovl->path.bbpos - 1] = '\0';
                printf("%08d %s", ovl->bread, buffer);
#undef TEST_ALN_OUT
#ifdef TEST_ALN_OUT
                {
                    tps = ((ovl->path.aepos-1)/tspace - ovl->path.abpos/tspace);
                    if (small)
                        Decompress_TraceTo16(ovl);
                    Load_Read(db1, ovl->aread, abuffer, 0);
                    Load_Read(db2, ovl->bread, bbuffer, 0);
                    if (COMP(aln->flags))
                        Complement_Seq(bbuffer, aln->blen);
                    Compute_Trace_PTS(aln,work,tspace);
                    int tlen  = aln->path->tlen;
                    int *trace = aln->path->trace;
                    int u;
                    printf(" ");
                    for (u = 0; u < tlen; u++)
                        printf("%d,", (int16) trace[u]);
                }
#endif
                printf("\n");
                if (SKIP == 1) {  //if SKIP = 0, then skip_rest is always 0
                    if ( ((int64) aln->alen < (int64) aln->blen) && ((int64) ovl->path.abpos < 1) && ((int64) aln->alen - (int64) ovl->path.aepos < 1) ) {
                        printf("* *\n");
                        skip_rest = 1;
                    }
                }
            }
         }

        if (ALIGN || CARTOON || REFERENCE)
          { if (ALIGN || REFERENCE)
              { char *aseq, *bseq;
                int   amin,  amax;
                int   bmin,  bmax;

                if (FLIP)
                  Flip_Alignment(aln,0);
                if (small)
                  Decompress_TraceTo16(ovl);

                amin = ovl->path.abpos - BORDER;
                if (amin < 0) amin = 0;
                amax = ovl->path.aepos + BORDER;
                if (amax > aln->alen) amax = aln->alen;
                if (COMP(aln->flags))
                  { bmin = (aln->blen-ovl->path.bepos) - BORDER;
                    if (bmin < 0) bmin = 0;
                    bmax = (aln->blen-ovl->path.bbpos) + BORDER;
                    if (bmax > aln->blen) bmax = aln->blen;
                  }
                else
                  { bmin = ovl->path.bbpos - BORDER;
                    if (bmin < 0) bmin = 0;
                    bmax = ovl->path.bepos + BORDER;
                    if (bmax > aln->blen) bmax = aln->blen;
                  }

                aseq = Load_Subread(db1,ovl->aread,amin,amax,abuffer,0);
                bseq = Load_Subread(db2,ovl->bread,bmin,bmax,bbuffer,0);

                aln->aseq = aseq - amin;
                if (COMP(aln->flags))
                  { Complement_Seq(bseq,bmax-bmin);
                    aln->bseq = bseq - (aln->blen - bmax);
                  }
                else
                  aln->bseq = bseq - bmin;

                Compute_Trace_PTS(aln,work,tspace);

                if (FLIP)
                  { if (COMP(aln->flags))
                      { Complement_Seq(aseq,amax-amin);
                        Complement_Seq(bseq,bmax-bmin);
                        aln->aseq = aseq - (aln->alen - amax);
                        aln->bseq = bseq - bmin;
                      }
                    Flip_Alignment(aln,1);
                  }
              }
            if (CARTOON)
              { printf("  (");
                Print_Number(tps,tp_wide,stdout);
                printf(" trace pts)\n\n");
                Alignment_Cartoon(stdout,aln,INDENT,mx_wide);
              }
            else
              { printf(" :   = ");
                Print_Number((int64) ovl->path.diffs,mn_wide,stdout);
                printf(" diffs  (");
                Print_Number(tps,tp_wide,stdout);
                printf(" trace pts)\n");
              }
            if (REFERENCE)
              Print_Reference(stdout,aln,work,INDENT,WIDTH,BORDER,UPPERCASE,mx_wide);
            if (ALIGN)
              Print_Alignment(stdout,aln,work,INDENT,WIDTH,BORDER,UPPERCASE,mx_wide);
          }
        else if (!(FALCON || M4OVL) )
          { printf(" :   < ");
            Print_Number((int64) ovl->path.diffs,mn_wide,stdout);
            printf(" diffs  (");
            Print_Number(tps,tp_wide,stdout);
            printf(" trace pts)\n");
          }
      }

    if (FALCON)
      { printf("+ +\n");
        printf("- -\n");
      }


    free(trace);
    if (ALIGN || FALCON)
      { free(bbuffer-1);
        free(abuffer-1);
        Free_Work_Data(work);
      }
  }

  Close_DB(db1);
  if (ISTWO)
    Close_DB(db2);

  exit (0);
}
