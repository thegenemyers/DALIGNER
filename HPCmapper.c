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
 *  Date  :  December 31, 2014
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
  { "[-vb] [-k<int(20)>] [-w<int(6)>] [-h<int(50)>] [-t<int>] [-M<int>]",
    "      [-e<double(.85)] [-l<int(1000)>] [-s<int(100)] [-H<int>]",
    "      [-m<track>]+ [-dal<int(4)>] [-deg<int(25)>]",
    "      <ref:db|dam> <reads:db|dam> [<first:int>[-<last:int>]]"
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
{ int   nblocks1, nblocks2;
  int   useblock1, useblock2;
  int   fblock, lblock;
#ifdef LSF
  int   jobid;
#endif

  char *pwd1, *root1;
  char *pwd2, *root2;

  int    MUNIT, DUNIT;
  int    VON, BON, CON;
  int    WINT, TINT, HGAP, HINT, KINT, SINT, LINT, MINT;
  double EREL;
  int    MMAX, MTOP;
  char **MASK;

  { int    i, j, k;         //  Process options
    int    flags[128];
    char  *eptr;

    ARG_INIT("HPCmapper")

    DUNIT = 4;
    MUNIT = 25;
    KINT  = 20;
    WINT  = 6;
    HINT  = 50;
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
            ARG_FLAGS("vbc");
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
    CON = flags['c'];

    if (argc < 3 || argc > 4)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[3]);
        exit (1);
      }
  }

  //  Make sure DAM and DB exist and the DB is partitioned, get number of blocks in partition

  pwd1   = PathTo(argv[1]);
  if (strcmp(argv[1]+(strlen(argv[1])-4),".dam") == 0)
    root1 = Root(argv[1],".dam");
  else
    root1 = Root(argv[1],".db");

  { int i, nfiles;
    FILE *dbvis;

    dbvis = fopen(Catenate(pwd1,"/",root1,".dam"),"r");
    if (dbvis == NULL)
      { dbvis = Fopen(Catenate(pwd1,"/",root1,".db"),"r");
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

    useblock1 = 1;
    if (fscanf(dbvis,"blocks = %d\n",&nblocks1) != 1)
      { useblock1 = 0;
        nblocks1  = 1;
      }

    fclose(dbvis);
  }

  pwd2   = PathTo(argv[2]);
  if (strcmp(argv[2]+(strlen(argv[2])-4),".dam") == 0)
    root2 = Root(argv[2],".dam");
  else
    root2 = Root(argv[2],".db");

  if (strcmp(root2,root1) == 0 && strcmp(pwd1,pwd2) == 0)
    { fprintf(stderr,"%s: Comparing the same data base %s/%s against itself, use HPCdaligner\n",
                     Prog_Name,pwd1,root1);
      exit (1);
    }

  { int i, nfiles;
    FILE *dbvis;

    dbvis = fopen(Catenate(pwd2,"/",root2,".dam"),"r");
    if (dbvis == NULL)
      { dbvis = Fopen(Catenate(pwd2,"/",root2,".db"),"r");
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

    useblock2 = 1;
    if (fscanf(dbvis,"blocks = %d\n",&nblocks2) != 1)
      { useblock2 = 0;
        nblocks2  = 1;
      }

    fclose(dbvis);
  }

  //  Set range fblock-lblock checking that DB.<fblock-1>.las exists & DB.<fblock>.las does not

  { char *eptr, *fptr;
    FILE *file;

    if (argc == 4)
      { fblock = strtol(argv[3],&eptr,10);
        if (*eptr != '\0' && *eptr != '-')
          { fprintf(stderr,"%s: final argument '%s' does not start with an integer\n",
                           Prog_Name,argv[3]);
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
        if (fblock < 1 || lblock > nblocks2 || fblock > lblock)
          { fprintf(stderr,"%s: range %d-%d is empty or out of bounds\n",Prog_Name,fblock,lblock);
            exit (1);
          }
      }
    else
      { fblock = 1;
        lblock = nblocks2;
      }

    if (fblock > 1)
      { file = fopen(Catenate(root1,".",root2,Numbered_Suffix(".",fblock-1,".las")),"r");
        if (file == NULL)
          { fprintf(stderr,"%s: File %s.%s.%d.las should already be present!\n",
                           Prog_Name,root1,root2,fblock-1);
            exit (1);
          }
        else
          fclose(file);
      }
    if (useblock2)
      { file = fopen(Catenate(root1,".",root2,Numbered_Suffix(".",fblock,".las")),"r");
        if (file != NULL)
          { fprintf(stderr,"%s: File %s.%s.%d.las should not yet exist!\n",
                           Prog_Name,root1,root2,fblock);
            exit (1);
          }
      }
    else
      { file = fopen(Catenate(root1,".",root2,".las"),"r");
        if (file != NULL)
          { fprintf(stderr,"%s: File %s.%s.las should not yet exist!\n",
                           Prog_Name,root1,root2);
            exit (1);
          }
      }
  }

  { int level, njobs;
    int i, j, k;
    int usepath1, usepath2;

    //  Produce all necessary daligner jobs ...

    usepath1 = (strcmp(pwd1,".") != 0);
    usepath2 = (strcmp(pwd2,".") != 0);

    njobs = nblocks1 * ( (lblock-fblock)/DUNIT + 1);

    printf("# Daligner jobs (%d)\n",njobs);

#ifdef LSF
    jobid = 1;
#endif
    for (i = 1; i <= nblocks1; i++)
      { int bits;
        int low, hgh;

        bits = (lblock-fblock)/DUNIT+1;
        low  = fblock;
        for (j = 1; j <= bits; j++)
          {
#ifdef LSF
            printf(LSF_ALIGN,jobid++);
            printf(" \"");
#endif
            printf("daligner -A");
            if (VON)
              printf(" -v");
            if (BON)
              printf(" -b");
            printf(" -k%d",KINT);
            if (WINT != 6)
              printf(" -w%d",WINT);
            printf(" -h%d",HINT);
            if (TINT > 0)
              printf(" -t%d",TINT);
            if (HGAP > 0)
              printf(" -H%d",HGAP);
            if (EREL > .1)
              printf(" -e%g",EREL);
            else
              printf(" -e.85");
            if (LINT != 1000)
              printf(" -l%d",LINT);
            if (SINT != 100)
              printf(" -s%d",SINT);
            if (MINT >= 0)
              printf(" -M%d",MINT);
            for (k = 0; k < MTOP; k++)
              printf(" -m%s",MASK[k]);
            if (useblock1)
              if (usepath1)
                printf(" %s/%s.%d",pwd1,root1,i);
              else
                printf(" %s.%d",root1,i);
            else
              if (usepath1)
                printf(" %s/%s",pwd1,root1);
              else
                printf(" %s",root1);
	    hgh = fblock + (((lblock-fblock)+1)*j)/bits;
            for (k = low; k < hgh; k++)
              if (useblock2)
                if (usepath2)
                  printf(" %s/%s.%d",pwd2,root2,k);
                else
                  printf(" %s.%d",root2,k);
              else
                if (usepath2)
                  printf(" %s/%s",pwd2,root2);
                else
                  printf(" %s",root2);
#ifdef LSF
            printf("\"");
#endif
            printf("\n");
            low = hgh;
          }
      }

    //  ... and then all the initial sort & merge jobs for each block pair

    printf("# Initial sort jobs (%d)\n", nblocks1*((lblock-fblock)+1));

#ifdef LSF
    jobid = 1;
#endif
    for (i = 1; i <= nblocks1; i++)
      for (j = fblock; j <= lblock; j++)
        {
#ifdef LSF
          printf(LSF_MERGE,jobid++);
          printf(" \"");
#endif
          printf("LAsort");
          if (VON)
            printf(" -v");
          if (CON)
            printf(" -c");
          for (k = 0; k < NTHREADS; k++)
            { if (useblock1)
                printf(" %s.%d",root1,i);
              else
                printf(" %s",root1);
              if (useblock2)
                printf(".%s.%d.C%d",root2,j,k);
              else
                printf(".%s.C%d",root2,k);
              if (useblock1)
                printf(" %s.%d",root1,i);
              else
                printf(" %s",root1);
              if (useblock2)
                printf(".%s.%d.N%d",root2,j,k);
              else
                printf(".%s.N%d",root2,k);
            }
          printf(" && LAmerge");
          if (VON)
            printf(" -v");
          if (CON)
            printf(" -c");
          if (nblocks1 == 1)
            if (useblock2)
              printf(" %s.%s.%d",root1,root2,j);
            else
              printf(" %s.%s",root1,root2);
          else
            printf(" L1.%d.%d",i,j);
          for (k = 0; k < NTHREADS; k++)
            { if (useblock1)
                printf(" %s.%d",root1,i);
              else
                printf(" %s",root1);
              if (useblock2)
                printf(".%s.%d.C%d.S",root2,j,k);
              else
                printf(".%s.C%d.S",root2,k);
              if (useblock1)
                printf(" %s.%d",root1,i);
              else
                printf(" %s",root1);
              if (useblock2)
                printf(".%s.%d.N%d.S",root2,j,k);
              else
                printf(".%s.N%d.S",root2,k);
            }
          printf(" && rm");
          for (k = 0; k < NTHREADS; k++)
            { if (useblock1)
                printf(" %s.%d",root1,i);
              else
                printf(" %s",root1);
              if (useblock2)
                printf(".%s.%d.C%d.S.las",root2,j,k);
              else
                printf(".%s.C%d.S.las",root2,k);
              if (useblock1)
                printf(" %s.%d",root1,i);
              else
                printf(" %s",root1);
              if (useblock2)
                printf(".%s.%d.N%d.S.las",root2,j,k);
              else
                printf(".%s.N%d.S.las",root2,k);
              if (useblock1)
                printf(" %s.%d",root1,i);
              else
                printf(" %s",root1);
              if (useblock2)
                printf(".%s.%d.C%d.las",root2,j,k);
              else
                printf(".%s.C%d.las",root2,k);
              if (useblock1)
                printf(" %s.%d",root1,i);
              else
                printf(" %s",root1);
              if (useblock2)
                printf(".%s.%d.N%d.las",root2,j,k);
              else
                printf(".%s.N%d.las",root2,k);
            }
#ifdef LSF
          printf("\"");
#endif
          printf("\n");
        }

    //  Higher level merges (if lblock > 1)

    if (nblocks1 > 1)
      { int pow, mway;

        //  Determine most balance mway for merging in ceil(log_mrg nblock1) levels

        pow = 1;
        for (level = 0; pow < nblocks1; level++)
          pow *= MUNIT;

        for (mway = MUNIT; mway >= 3; mway--)
          if (power(mway,level) < nblocks1)
            break;
        mway += 1;

        //  Issue the commands for each merge level

        { int  p, cnt;

          cnt = nblocks1;
          for (i = 1; i <= level; i++)
            { int bits;
              int low, hgh;

              bits = (cnt-1)/mway+1;
              printf("# Level %d jobs (%d)\n",i,bits*((lblock-fblock)+1));

              //  Block merges

#ifdef LSF
              jobid = 1;
#endif
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
                      if (CON)
                        printf(" -c");
                      if (i == level)
                        if (useblock2)
                          printf(" %s.%s.%d",root1,root2,j);
                        else
                          printf(" %s.%s",root1,root2);
                      else
                        printf(" L%d.%d.%d",i+1,j,p);
                      for (k = low; k <= hgh; k++)
                        printf(" L%d.%d.%d",i,k,j);
                      printf(" && rm");
                      for (k = low; k <= hgh; k++)
                        printf(" L%d.%d.%d.las",i,k,j);
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

  free(root2);
  free(pwd2);
  free(root1);
  free(pwd1);

  exit (0);
}
