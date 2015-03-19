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

/*********************************************************************************************\
 *
 *  Produce a script to compute overlaps for all block pairs of a DB, and then sort and merge
 *    them into as many .las files as their are blocks.
 *
 *  Author:  Gene Myers
 *  Date  :  June 1, 2014
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"
#include "filter.h"

#undef  LSF  //  define if want a directly executable LSF script

static char *Usage[] =
  { "[-vbAI] [-k<int(14)>] [-w<int(6)>] [-h<int(35)>] [-t<int>] [-M<int>]",
    "        [-e<double(.70)] [-l<int(1000)>] [-s<int(100)] [-H<int>]",
    "        [-m<track>]+ [-dal<int(4)>] [-deg<int(25)>]",
    "        <path:db|dam> [<block:int>[-<range:int>]"
  };

static int power(int base, int exp)
{ int i, pow;

  pow = 1;
  for (i = 0; i < exp; i++)
    pow *= base;
  return (pow);
}

#define LSF_ALIGN "bsub -q medium -n 4 -o ALIGN.out -e ALIGN.err -R span[hosts=1] -J align#%d"
#define LSF_MERGE "bsub -q short -n 12 -o MERGE.out -e MERGE.err -R span[hosts=1] -J merge#%d"

int main(int argc, char *argv[])
{ int   nblocks;
  int   useblock;
  int   fblock, lblock;
#ifdef LSF
  int   jobid;
#endif

  char *pwd, *root;

  int    MUNIT, DUNIT;
  int    VON, BON, AON, ION;
  int    WINT, TINT, HGAP, HINT, KINT, SINT, LINT, MINT;
  double EREL;
  int    MMAX, MTOP;
  char **MASK;

  { int    i, j, k;         //  Process options
    int    flags[128];
    char  *eptr;

    ARG_INIT("HPCdaligner")

    DUNIT = 4;
    MUNIT = 25;
    KINT  = 14;
    WINT  = 6;
    HINT  = 35;
    TINT  = 0;
    HGAP  = 0;
    EREL  = 0.;
    LINT  = 1000;
    SINT  = 100;
    MINT  = -1;

    MTOP = 0;
    MMAX = 10;
    MASK = (char **) Malloc(MMAX*sizeof(char *),"Allocating mask track array");
    if (MASK == NULL)
      exit (1);

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vbAI");
            break;
          case 'k':
            ARG_POSITIVE(KINT,"K-mer length")
            break;
          case 'w':
            ARG_POSITIVE(WINT,"Log of bin width")
            break;
          case 'h':
            ARG_POSITIVE(HINT,"Hit threshold (in bp.s)")
            break;
          case 't':
            ARG_POSITIVE(TINT,"Tuple suppression frequency")
            break;
          case 'H':
            ARG_POSITIVE(HGAP,"HGAP threshold (in bp.s)")
            break;
          case 'e':
            ARG_REAL(EREL)
            if (EREL < .7 || EREL >= 1.)
              { fprintf(stderr,"%s: Average correlation must be in [.7,1.) (%g)\n",Prog_Name,EREL);
                exit (1);
              }
            break;
          case 'l':
            ARG_POSITIVE(LINT,"Minimum ovlerap length")
            break;
          case 's':
            ARG_POSITIVE(SINT,"Trace spacing")
            break;
          case 'M':
            ARG_NON_NEGATIVE(MINT,"Memory allocation (in Gb)")
            break;
          case 'm':
            if (MTOP >= MMAX)
              { MMAX = 1.2*MTOP + 10;
                MASK = (char **) Realloc(MASK,MMAX*sizeof(char *),"Reallocating mask track array");
                if (MASK == NULL)
                  exit (1);
              }
            MASK[MTOP++] = argv[i]+2;
            break;
          case 'd':
            if (argv[i][2] == 'e' && argv[i][3] == 'g')
              { MUNIT = strtol(argv[i]+4,&eptr,10);
                if (*eptr != '\0' || argv[i][4] == '\0')
                  { fprintf(stderr,"%s: -mrg argument is not an integer\n",Prog_Name);
                    exit (1);
                  }
                if (MUNIT <= 0)
                  { fprintf(stderr,"%s: Files per merge must be positive (%d)\n",
                                   Prog_Name,MUNIT);
                    exit (1);
                  }
                if (MUNIT < 3)
                  { fprintf(stderr,"%s: Files per merge must be at least 3 (%d)\n",
                                   Prog_Name,MUNIT);
                    exit (1);
                  }
              }
            else if (argv[i][2] == 'a' && argv[i][3] == 'l')
              { DUNIT = strtol(argv[i]+4,&eptr,10);
                if (*eptr != '\0' || argv[i][4] == '\0')
                  { fprintf(stderr,"%s: -dal argument is not an integer\n",Prog_Name);
                    exit (1);
                  }
                if (DUNIT <= 0)
                  { fprintf(stderr,"%s: Blocks per daligner call must be positive (%d)\n",
                                   Prog_Name,DUNIT);
                    exit (1);
                  }
              }
            else
              { fprintf(stderr,"%s: -%.3s is an illegal option\n",Prog_Name,argv[i]+1);
                exit (1);
              }
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VON = flags['v'];
    BON = flags['b'];
    AON = flags['A'];
    ION = flags['I'];

    if (argc < 2 || argc > 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[3]);
        exit (1);
      }
  }

  //  Make sure DB exists and is partitioned, get number of blocks in partition

  pwd = PathTo(argv[1]);
  if (strcmp(argv[1]+(strlen(argv[1])-4),".dam") == 0)
    root = Root(argv[1],".dam");
  else
    root = Root(argv[1],".db");

  { int i, nfiles;
    FILE *dbvis;

    dbvis = fopen(Catenate(pwd,"/",root,".dam"),"r");
    if (dbvis == NULL)
      { dbvis = Fopen(Catenate(pwd,"/",root,".db"),"r");
        if (dbvis == NULL)
          exit (1);
      }

    if (fscanf(dbvis,"files = %d\n",&nfiles) != 1)
      SYSTEM_ERROR
    for (i = 0; i < nfiles; i++)
      { char buffer[30001];

        if (fgets(buffer,30000,dbvis) == NULL)
          SYSTEM_ERROR
      }

    useblock = 1;
    if (fscanf(dbvis,"blocks = %d\n",&nblocks) != 1)
      { useblock = 0;
        nblocks  = 1;
      }
  }

  //  Set range fblock-lblock checking that DB.<fblock-1>.las exists & DB.<fblock>.las does not

  { char *eptr, *fptr;
    FILE *file;

    if (argc == 3)
      { fblock = strtol(argv[2],&eptr,10);
        if (*eptr != '\0' && *eptr != '-')
          { fprintf(stderr,"%s: final argument '%s' does not start with an integer\n",
                           Prog_Name,argv[2]);
            exit (1);
          }
        if (*eptr == '-')
          { lblock = strtol(eptr+1,&fptr,10);
            if (*fptr != '\0')
              { fprintf(stderr,"%s: second part of range '%s' is not an integer\n",
                               Prog_Name,eptr+1);
                exit (1);
              }
          }
        else
          lblock = fblock;
        if (fblock < 1 || lblock > nblocks || fblock > lblock)
          { fprintf(stderr,"%s: range %d-%d is empty or out of bounds\n",Prog_Name,fblock,lblock);
            exit (1);
          }
      }
    else
      { fblock = 1;
        lblock = nblocks;
      }

    if (fblock > 1)
      { file = fopen(Catenate(root,Numbered_Suffix(".",fblock-1,".las"),"",""),"r");
        if (file == NULL)
          { fprintf(stderr,"%s: File %s.%d.las should already be present!\n",
                           Prog_Name,root,fblock-1);
            exit (1);
          }
        else
          fclose(file);
      }
    file = fopen(Catenate(root,Numbered_Suffix(".",fblock,".las"),"",""),"r");
    if (file != NULL)
      { fprintf(stderr,"%s: File %s.%d.las should not yet exist!\n",
                       Prog_Name,root,fblock);
        exit (1);
      }
  }

  { int level, njobs;
    int i, j, k;
    int usepath;

    //  Produce all necessary daligner jobs ...

    usepath = (strcmp(pwd,".") != 0);

    njobs = 0;
    for (i = fblock; i <= lblock; i++)
      njobs += (i-1)/DUNIT+1;

    printf("# Daligner jobs (%d)\n",njobs);

#ifdef LSF
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; i++)
      { int bits;
        int low, hgh;

        bits = (i-1)/DUNIT+1;
        low  = 1;
        for (j = 1; j <= bits; j++)
          {
#ifdef LSF
            printf(LSF_ALIGN,jobid++);
            printf(" \"");
#endif
            printf("daligner");
            if (VON)
              printf(" -v");
            if (BON)
              printf(" -b");
            if (AON)
              printf(" -A");
            if (ION)
              printf(" -I");
            if (KINT != 14)
              printf(" -k%d",KINT);
            if (WINT != 6)
              printf(" -w%d",WINT);
            if (HINT != 35)
              printf(" -h%d",HINT);
            if (TINT > 0)
              printf(" -t%d",TINT);
            if (HGAP > 0)
              printf(" -H%d",HGAP);
            if (EREL > .1)
              printf(" -e%g",EREL);
            if (LINT != 1000)
              printf(" -l%d",LINT);
            if (SINT != 100)
              printf(" -s%d",SINT);
            if (MINT >= 0)
              printf(" -M%d",MINT);
            for (k = 0; k < MTOP; k++)
              printf(" -m%s",MASK[k]);
            if (useblock)
              if (usepath)
                printf(" %s/%s.%d",pwd,root,i);
              else
                printf(" %s.%d",root,i);
            else
              if (usepath)
                printf(" %s/%s",pwd,root);
              else
                printf(" %s",root);
            hgh = (i*j)/bits + 1;
            for (k = low; k < hgh; k++)
              if (useblock)
                if (usepath)
                  printf(" %s/%s.%d",pwd,root,k);
                else
                  printf(" %s.%d",root,k);
              else
                if (usepath)
                  printf(" %s/%s",pwd,root);
                else
                  printf(" %s",root);
#ifdef LSF
            printf("\"");
#endif
            printf("\n");
            low = hgh;
          }
      }

    //  ... and then all the initial sort & merge jobs for each block pair

    printf("# Initial sort jobs (%d)\n", lblock*lblock - (fblock-1)*(fblock-1) );

#ifdef LSF
    jobid = 1;
#endif
    for (i = 1; i <= lblock; i++)
      for (j = (i < fblock ? fblock : 1); j <= lblock; j++)
        {
#ifdef LSF
          printf(LSF_MERGE,jobid++);
          printf(" \"");
#endif
          printf("LAsort");
          if (VON)
            printf(" -v");
          for (k = 0; k < NTHREADS; k++)
            if (useblock)
              { printf(" %s.%d.%s.%d.C%d",root,i,root,j,k);
                printf(" %s.%d.%s.%d.N%d",root,i,root,j,k);
              }
            else
              { printf(" %s.%s.C%d",root,root,k);
                printf(" %s.%s.N%d",root,root,k);
              }
          printf(" && LAmerge");
          if (VON)
            printf(" -v");
          if (lblock == 1)
            printf(" %s.%d",root,i);
          else if (i < fblock)
            printf(" L1.%d.%d",i,(j-fblock)+1);
          else
            printf(" L1.%d.%d",i,j);
          for (k = 0; k < NTHREADS; k++)
            if (useblock)
              { printf(" %s.%d.%s.%d.C%d.S",root,i,root,j,k);
                printf(" %s.%d.%s.%d.N%d.S",root,i,root,j,k);
              }
            else
              { printf(" %s.%s.C%d.S",root,root,k);
                printf(" %s.%s.N%d.S",root,root,k);
              }
          printf(" && rm");
          for (k = 0; k < NTHREADS; k++)
            if (useblock)
              { printf(" %s.%d.%s.%d.C%d.S.las",root,i,root,j,k);
                printf(" %s.%d.%s.%d.N%d.S.las",root,i,root,j,k);
              }
            else
              { printf(" %s.%s.C%d.S.las",root,root,k);
                printf(" %s.%s.N%d.S.las",root,root,k);
              }
#ifdef LSF
          printf("\"");
#endif
          printf("\n");
        }

    //  Higher level merges (if lblock > 1)

    if (lblock > 1)
      { int pow, mway;

        //  Determine most balance mway for merging in ceil(log_mrg lblock) levels

        pow = 1;
        for (level = 0; pow < lblock; level++)
          pow *= MUNIT;

        for (mway = MUNIT; mway >= 3; mway--)
          if (power(mway,level) < lblock)
            break;
        mway += 1;

        //  Issue the commands for each merge level

        { int  p, cnt, dnt;

          cnt = lblock;
          dnt = (lblock-fblock)+1;
          for (i = 1; i <= level; i++)
            { int bits, dits;
              int low, hgh;

              bits = (cnt-1)/mway+1;
              dits = (dnt-1)/mway+1;

              //  Incremental update merges

#ifdef LSF
              jobid = 1;
#endif
              if (dnt >= 1)
                { int last;

                  last = (dnt == 1 || i == level);
                  printf("# Level %d jobs (%d)\n",i,bits*((lblock-fblock)+1) + dits*(fblock-1));
                  for (j = 1; j < fblock; j++)
                    {
#ifdef LSF
                      printf(LSF_MERGE,jobid++);
                      printf(" \"");
#endif
                      if (last)
                        printf("mv %s.%d.las L%d.%d.0.las && ",root,j,i,j);
                      low = 1;
                      for (p = 1; p <= dits; p++)
                        { hgh = (dnt*p)/dits;
#ifdef LSF
                          if (p > 1)
                            { printf(LSF_MERGE,jobid++);
                              printf(" \"");
                            }
#endif
                          printf("LAmerge");
                          if (VON)
                            printf(" -v");
                          if (last)
                            printf(" %s.%d L%d.%d.0",root,j,i,j);
                          else
                            printf(" L%d.%d.%d",i+1,j,p);
                          for (k = low; k <= hgh; k++)
                            printf(" L%d.%d.%d",i,j,k);
                          printf(" && rm");
                          if (last)
                            printf(" L%d.%d.0.las",i,j);
                          for (k = low; k <= hgh; k++)
                            printf(" L%d.%d.%d.las",i,j,k);
#ifdef LSF
                          printf("\"");
#endif
                          printf("\n");
                          low = hgh+1;
                        }
                    }
                  if (dnt > 1)
                    dnt = dits;
                  else
                    dnt = 0;
                }
              else
                printf("# Level %d jobs (%d)\n",i,bits*((lblock-fblock)+1));

              //  New block merges

              for (j = fblock; j <= lblock; j++) 
                { low = 1;
                  for (p = 1; p <= bits; p++)
                    { hgh = (cnt*p)/bits;
#ifdef LSF
                      printf(LSF_MERGE,jobid++);
                      printf(" \"");
#endif
                      printf("LAmerge");
                      if (VON)
                        printf(" -v");
                      if (i == level)
                        printf(" %s.%d",root,j);
                      else
                        printf(" L%d.%d.%d",i+1,j,p);
                      for (k = low; k <= hgh; k++)
                        printf(" L%d.%d.%d",i,j,k);
                      printf(" && rm");
                      for (k = low; k <= hgh; k++)
                        printf(" L%d.%d.%d.las",i,j,k);
#ifdef LSF
                      printf("\"");
#endif
                      printf("\n");
                      low = hgh+1;
                    }
                }

              cnt  = bits;
            }
        }
    }
  }

  free(root);
  free(pwd);

  exit (0);
}
