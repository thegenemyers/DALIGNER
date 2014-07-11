/*******************************************************************************************
 *
 *  Utility for displaying the overlaps in a .las file in a variety of ways including
 *    a minimal listing of intervals, a cartoon, and a full out alignment.  One can
 *    also list only those overlaps whose A read is in a set of supplied read ranges.
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
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

int ORDER(const void *l, const void *r)
{ int x = *((int32 *) l);
  int y = *((int32 *) r);
  return (x-y);
}

static char *Usage[] =
    { "[-coU] [-(a|r):<db>] [-i<int(4)>] [-w<int(100)>] [-b<int(10)>] ",
      "       <align:las> [ <reads:range> ... ]"
    };

int main(int argc, char *argv[])
{ HITS_DB   _db,  *db  = &_db; 
  Overlap   _ovl, *ovl = &_ovl;
  Alignment _aln, *aln = &_aln;

  FILE   *input;
  int64   novl;
  int     tspace, tbytes, small;
  int     reps, *pts;

  int     ALIGN, CARTOON, OVERLAP, REFERENCE;
  int     INDENT, WIDTH, BORDER, UPPERCASE;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("LAshow")

    ALIGN     = 0;
    REFERENCE = 0;
    INDENT    = 4;
    WIDTH     = 100;
    BORDER    = 10;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("coU")
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
          case 'r':
            REFERENCE = 1;
          case 'a':
            ALIGN = 1;
            if (argv[i][2] != ':')
              { fprintf(stderr,"%s: Unrecognizable option %s\n",Prog_Name,argv[i]);
                exit (1);
              }
            if (Open_DB(argv[i]+3,db))
              exit (1);
            Trim_DB(db);
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    CARTOON   = flags['c'];
    OVERLAP   = flags['o'];
    UPPERCASE = flags['U'];

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        exit (1);
      }
  }

  //  Process read index arguments into a sorted list of read ranges

  pts  = (int *) Malloc(sizeof(int)*2*argc,"Allocating read parameters");
  if (pts == NULL)
    exit (1);

  reps = 0;
  if (argc > 2)
    { int   c, b, e;
      char *eptr, *fptr;

      for (c = 2; c < argc; c++)
        { if (argv[c][0] == '#')
            { fprintf(stderr,"%s: # is not allowed as range start, '%s'\n",
                      Prog_Name,argv[c]);
              exit (1);
            }
          else
            { b = strtol(argv[c],&eptr,10);
              if (b < 1)
                { fprintf(stderr,"%s: Non-positive index?, '%d'\n",Prog_Name,b);
                  exit (1);
                }
            }
          if (eptr > argv[c])
            { if (*eptr == '\0')
                { pts[reps++] = b;
                  pts[reps++] = b;
                  continue;
                }
              else if (*eptr == '-')
                { if (eptr[1] == '#')
                    { e = INT32_MAX;
                      fptr = eptr+2;
                    }
                  else
                    e = strtol(eptr+1,&fptr,10);
                  if (fptr > eptr+1 && *fptr == 0 && eptr[1] != '-')
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

  //  Initiate file reading and read (novl, tspace) header
  
  { char  *over, *pwd, *root;

    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".las");
    over  = Catenate(pwd,"/",root,".las");
    input = Fopen(over,"r");
    if (input == NULL)
      exit (1);

    fread(&novl,sizeof(int64),1,input);
    fread(&tspace,sizeof(int),1,input);

    if (tspace <= TRACE_XOVR)
      { small  = 1;
        tbytes = sizeof(uint8);
      }
    else
      { small  = 0;
        tbytes = sizeof(uint16);
      }

    printf("\n%s: ",root);
    Print_Number(novl,0,stdout);
    printf(" records\n");

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

    if (ALIGN)
      { work = New_Work_Data();
  
        aln->path = &(ovl->path);
        aln->aseq = New_Read_Buffer(db);
        aln->bseq = New_Read_Buffer(db);
      }
    else
      work = NULL;

    tmax  = 1000;
    trace = (uint16 *) Malloc(sizeof(uint16)*tmax,"Allocating trace vector");
    if (trace == NULL)
      exit (1);

    in  = 0;
    npt = pts[0];
    idx = 1;

    //  For each record do

    for (j = 0; j < novl; j++)

       //  Read it in

      { Read_Overlap(input,ovl);
        if (ovl->path.tlen > tmax)
          { tmax = 1.2*ovl->path.tlen + 100;
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

        if (OVERLAP)
          { if (ovl->path.abpos != 0 && ovl->path.bbpos != 0)
              continue;
            if (ovl->path.aepos != ovl->alen && ovl->path.bepos != ovl->blen)
              continue;
          }

        //  Display it

        if (CARTOON || ALIGN)
          printf("\n");
        Print_Number((int64) ovl->aread+1,10,stdout);
        printf("  ");
        Print_Number((int64) ovl->bread+1,9,stdout);
        if (COMP(ovl->flags))
          printf(" c");
        else
          printf(" n");
        printf("   [");
        Print_Number((int64) ovl->path.abpos,6,stdout);
        printf("..");
        Print_Number((int64) ovl->path.aepos,6,stdout);
        printf("] x [");
        Print_Number((int64) ovl->path.bbpos,6,stdout);
        printf("..");
        Print_Number((int64) ovl->path.bepos,6,stdout);
        printf("]");

/*
{ int u;
  if (small)
    for (u = 0; u < ovl->path.tlen-1; u++)
    // for (u = 1; u < ovl->path.tlen-2; u += 2)
      printf("  %3d\n",((uint8 *) ovl->path.trace)[u]);
  else
    for (u = 0; u < ovl->path.tlen-1; u++)
    // for (u = 1; u < ovl->path.tlen-2; u += 2)
      printf("  %3d\n",((uint16 *) ovl->path.trace)[u]);
}
*/

        // tps = ((ovl->path.aepos-1)/tspace - ovl->path.abpos/tspace) + 1;
        tps = ((ovl->path.aepos-1)/tspace - ovl->path.abpos/tspace);
        if (ALIGN)
          { if (small)
              Decompress_TraceTo16(ovl);
            aln->alen  = ovl->alen;
            aln->blen  = ovl->blen;
            aln->flags = ovl->flags;
            Load_Read(db,ovl->aread,aln->aseq,0);
            Load_Read(db,ovl->bread,aln->bseq,0);
            if (COMP(aln->flags))
              Complement_Seq(aln->bseq);
            Compute_Trace_PTS(aln,work,tspace);
            if (CARTOON)
              { printf("  (");
                Print_Number(tps,3,stdout);
                printf(" trace pts)\n\n");
                Print_ACartoon(stdout,aln,INDENT);
              }
            else
              { printf(" :   = ");
                Print_Number((int64) ovl->path.diffs,6,stdout);
                printf(" diffs  (");
                Print_Number(tps,3,stdout);
                printf(" trace pts)\n");
              }
            if (REFERENCE)
              Print_Reference(stdout,aln,work,INDENT,WIDTH,BORDER,UPPERCASE,5);
            else
              Print_Alignment(stdout,aln,work,INDENT,WIDTH,BORDER,UPPERCASE,5);
          }
        else if (CARTOON)
          { printf("  (");
            Print_Number(tps,3,stdout);
            printf(" trace pts)\n\n");
            Print_OCartoon(stdout,ovl,INDENT);
          }
        else
          { printf(" :   < ");
            Print_Number((int64) ovl->path.diffs,6,stdout);
            printf(" diffs  (");
            Print_Number(tps,3,stdout);
            printf(" trace pts)\n");
          }
      }

    free(trace);
    if (ALIGN)
      { free(aln->bseq-1);
        free(aln->aseq-1);
        Free_Work_Data(work);
      }
  }

  exit (0);
}
