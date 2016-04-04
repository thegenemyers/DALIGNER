/*********************************************************************************************\
 *
 *  Produce a script to compute overlaps for all block pairs between two DBs, and then sort
 *    and merge *    them into as many .las files as their are blocks of the 1st DB.
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
    "      [-e<double(.85)] [-l<int(1000)>] [-s<int(100)]",
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

#define LSF_ALIGN "bsub -q medium -n 4 -o MAPALL.out -e MAPALL.err -R span[hosts=1] -J align#%d"
#define LSF_SORT  "bsub -q short -n 12 -o SORT.ALL.out -e SORT.ALL.err -R span[hosts=1] -J sort#%d"
#define LSF_MERGE \
            "bsub -q short -n 12 -o MERGE.ALL%d.out -e MERGE.ALL%d.err -R span[hosts=1] -J merge#%d"

int main(int argc, char *argv[])
{ int   nblocks1, nblocks2;
  int   useblock1, useblock2;
  int   usepath1, usepath2;
  int   fblock, lblock;
#ifdef LSF
  int   jobid;
#endif

  char *pwd1, *root1;
  char *pwd2, *root2;

  int    MUNIT, DUNIT;
  int    VON, BON, CON;
  int    WINT, TINT, HINT, KINT, SINT, LINT, MINT;
  double EREL;
  int    MMAX, MTOP;
  char **MASK;

  { int    i, j, k;         //  Process options
    int    flags[128];
    char  *eptr;

    ARG_INIT("HPC.mapper")

    DUNIT = 4;
    MUNIT = 25;
    KINT  = 20;
    WINT  = 6;
    HINT  = 50;
    TINT  = 0;
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
    if (fscanf(dbvis,"blocks = %d\n",&nblocks1) != 1 || nblocks1 == 1)
      { useblock1 = 0;
        nblocks1  = 1;
      }

    usepath1 = (strcmp(pwd1,".") != 0);

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
    if (fscanf(dbvis,"blocks = %d\n",&nblocks2) != 1 || nblocks2 == 1)
      { useblock2 = 0;
        nblocks2  = 1;
      }

    usepath2 = (strcmp(pwd2,".") != 0);

    fclose(dbvis);
  }

  //  Set range fblock-lblock checking that DB.<fblock-1>.las exists & DB.<fblock>.las does not

  { char *eptr, *fptr, *src2;
    FILE *file;

    if (argc == 4)
      { fblock = strtol(argv[3],&eptr,10);
        if (*eptr != '\0' && *eptr != '-')
          { fprintf(stderr,"%s: final argument '%s' does not start with an integer\n",
                           Prog_Name,argv[3]);
            exit (1);
          }
        useblock2 = 1;
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

    if (usepath2)
      src2 = Strdup(Catenate(pwd2,"/",root2,""),"Allocating small string!");
    else
      src2 = Strdup(root2,"Allocating small string!");
    if (src2 == NULL)
      exit (1);

    if (fblock > 1)
      { file = fopen(Catenate(src2,".",root1,Numbered_Suffix(".",fblock-1,".las")),"r");
        if (file == NULL)
          { fprintf(stderr,"%s: File %s.%s.%d.las should already be present!\n",
                           Prog_Name,src2,root1,fblock-1);
            exit (1);
          }
        else
          fclose(file);
      }
    if (useblock2)
      { file = fopen(Catenate(src2,".",root1,Numbered_Suffix(".",fblock,".las")),"r");
        if (file != NULL)
          { fprintf(stderr,"%s: File %s.%s.%d.las should not yet exist!\n",
                           Prog_Name,src2,root1,fblock);
            exit (1);
          }
      }
    else
      { file = fopen(Catenate(src2,".",root1,".las"),"r");
        if (file != NULL)
          { fprintf(stderr,"%s: File %s.%s.las should not yet exist!\n",
                           Prog_Name,src2,root1);
            exit (1);
          }
      }

    free(src2);
  }

  { int level, njobs;
    int i, j, k, t;
    char orient[2] = { 'C', 'N' };

    //  Produce all necessary daligner jobs ...

    njobs = nblocks1 * ( (lblock-fblock)/DUNIT + 1);

    printf("# Daligner jobs (%d)\n",njobs);

#ifdef LSF
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; i++)
      { int bits;
        int low, hgh;

        bits = (nblocks1-1)/DUNIT+1;
        low  = 1;
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

            printf(" ");
            if (usepath2)
              printf("%s/",pwd2);
            printf("%s",root2);
            if (useblock2)
              printf(".%d",i);

	    hgh = 1 + (nblocks1*j)/bits;
            for (k = low; k < hgh; k++)
              { printf(" ");
                if (usepath1)
                  printf("%s/",pwd1);
                printf("%s",root1);
                if (useblock1)
                  printf(".%d",k);
              }
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
    for (j = fblock; j <= lblock; j++)
      for (i = 1; i <= nblocks1; i++)
        {
#ifdef LSF
          printf(LSF_SORT,jobid++);
          printf(" \"");
#endif
          printf("LAsort ");
          if (VON)
            printf("-v ");
          if (CON)
            printf("-c ");
          for (k = 0; k < NTHREADS; k++)
            for (t = 0; t < 2; t++)
              { printf("%s",root2);
                if (useblock2)
                  printf(".%d",j);
                printf(".%s",root1);
                if (useblock1)
                  printf(".%d",i);
                printf(".%c%d ",orient[t],k);
              }

          printf("&& LAmerge ");
          if (VON)
            printf("-v ");
          if (CON)
            printf("-c ");
          if (nblocks1 == 1)
            { if (usepath2)
                printf("%s/",pwd2);
              printf("%s.%s",root2,root1);
              if (useblock2)
                printf(".%d",j);
            }
          else
            printf("L1.%d.%d",j,i);

          for (k = 0; k < NTHREADS; k++)
            for (t = 0; t < 2; t++)
              { printf(" %s",root2);
                if (useblock2)
                  printf(".%d",j);
                printf(".%s",root1);
                if (useblock1)
                  printf(".%d",i);
                printf(".%c%d.S",orient[t],k);
              }
#ifdef LSF
          printf("\"");
#endif
          printf("\n");
        }

    //  Check .las files (optional)

    printf("# Check all level 1 .las files (optional but recommended)\n");

    for (j = fblock; j <= lblock; j++)
      for (i = 1; i <= nblocks1; i++)
        { printf("LAcheck -vS ");
          if (usepath2)
            printf("%s/%s ",pwd2,root2);
          else
            printf("%s ",root2);
          if (usepath1)
            printf("%s/%s ",pwd1,root1);
          else
            printf("%s ",root1);
          if (nblocks1 == 1)
            { if (usepath2)
                printf("%s/",pwd2);
              printf("%s.%s",root2,root1);
              if (useblock2)
                printf(".%d",j);
            }
          else
            printf("L1.%d.%d",j,i);
          printf("\n");
        }

    //  Clean up (optional)

    printf("# Remove initial .las files (optional)\n");

    for (j = fblock; j <= lblock; j++)
      for (i = 1; i <= nblocks1; i++)
        for (t = 0; t < 4; t++)
          { printf("rm");
            for (k = 0; k < NTHREADS; k++)
              { printf(" %s",root2);
                if (useblock2)
                  printf(".%d",j);
                printf(".%s",root1);
                if (useblock1)
                  printf(".%d",i);
                printf(".%c%d",orient[t%2],k);
                if (t >= 2)
                  printf(".S");
                printf(".las");
              }
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
                      printf(LSF_MERGE,i,i,jobid++);
                      printf(" \"");
#endif
                      printf("LAmerge ");
                      if (VON)
                        printf("-v ");
                      if (CON)
                        printf("-c ");
                      if (i == level)
                        { if (usepath2)
                            printf("%s/",pwd2);
                          printf("%s.%s",root2,root1);
                          if (useblock2)
                            printf(".%d",j);
                        }
                      else
                        printf("L%d.%d.%d",j,i+1,p);
                      for (k = low; k <= hgh; k++)
                        printf(" L%d.%d.%d",i,j,k);

#ifdef LSF
                      printf("\"");
#endif
		      printf("\n");
                      low = hgh+1;
                    }
                }

              //  Check new .las (optional)

              printf("# Check all level %d .las files (optional but recommended)\n",i+1);

              for (j = fblock; j <= lblock; j++) 
                { low = 1;
                  for (p = 1; p <= bits; p++)
                    { hgh = (cnt*p)/bits;
                      printf("LAcheck -vS ");
                      if (usepath2)
                        printf("%s/%s ",pwd2,root2);
                      else
                        printf("%s ",root2);
                      if (usepath1)
                        printf("%s/%s ",pwd1,root1);
                      else
                        printf("%s ",root1);
                      if (i == level)
                        { if (usepath2)
                            printf("%s/",pwd2);
                          printf("%s.%s",root2,root1);
                          if (useblock2)
                            printf(".%d",j);
                        }
                      else
                        printf("L%d.%d.%d",j,i+1,p);
                      printf("\n");
                      low = hgh+1;
                    }
                }

              //  Cleanup (optional)

              printf("# Remove level %d .las files (optional)\n",i);

              for (j = fblock; j <= lblock; j++) 
                { low = 1;
                  for (p = 1; p <= bits; p++)
                    { hgh = (cnt*p)/bits;
                      printf("rm");
                      for (k = low; k <= hgh; k++)
                        printf(" L%d.%d.%d.las",i,j,k);
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
