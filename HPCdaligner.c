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

static char *Usage[] =
  { "[-vbd] [-k<int(14)>] [-w<int(6)>] [-h<int(35)>] [-t<int>] [-H<int>]",
    "       [-e<double(.70)] [-m<double(.55)>] [-l<int(1000)>]",
    "       [-s<int(100)>] [-dal<int(4)>] [-mrg<int(25)>]",
    "       <path:db> [<block:int>[-<range:int>]"
  };

int power(int base, int exp)
{ int i, pow;

  pow = 1;
  for (i = 0; i < exp; i++)
    pow *= base;
  return (pow);
}

int main(int argc, char *argv[])
{ int   nblocks;
  int   fblock, lblock;

  char *pwd, *root;

  int    MUNIT, DUNIT;
  int    VON, BON, DON;
  int    WINT, TINT, HGAP, HINT, KINT, SINT, LINT;
  double EREL, MREL;

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
    EREL  = .70;
    MREL  = .55;
    SINT  = 100;
    LINT  = 1000;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
          optflags:
            ARG_FLAGS("vbd");
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
          case 's':
            ARG_POSITIVE(SINT,"Trace spacing")
            break;
          case 'l':
            ARG_POSITIVE(LINT,"Minimum ovlerap length")
            break;
          case 'm':
            if (argv[i][2] == 'r' && argv[i][3] == 'g')
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
                break;
              }
            ARG_REAL(MREL)
            if (MREL < .55 || MREL >= 1.)
              { fprintf(stderr,"%s: Minimum correlation must be in [.55,1.) (%g)\n",Prog_Name,MREL);
                exit (1);
              }
            break;
          case 'd':
            if (argv[i][2] == 'a' && argv[i][3] == 'l')
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
              goto optflags;
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VON = flags['v'];
    BON = flags['b'];
    DON = flags['d'];

    if (argc < 2 || argc > 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[3]);
        exit (1);
      }
  }

  //  Make sure DB exists and is partitioned, get number of blocks in partition

  pwd  = PathTo(argv[1]);
  root = Root(argv[1],".db");

  { int i, nfiles;
    FILE *dbvis;

    dbvis = Fopen(Catenate(pwd,"/",root,".db"),"r");
    if (dbvis == NULL)
      exit (1);

    fscanf(dbvis,"files = %d\n",&nfiles);
    for (i = 0; i < nfiles; i++)
      { char buffer[30001];
        fgets(buffer,30000,dbvis);
      }

    if (fscanf(dbvis,"blocks = %d\n",&nblocks) != 1)
      { fprintf(stderr,"%s: Database %s has not been split\n",Prog_Name,root);
        exit (0);
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

    for (i = fblock; i <= lblock; i++)
      { int bits;
        int low, hgh;

        bits = (i-1)/DUNIT+1;
        low  = 1;
        for (j = 1; j <= bits; j++)
          { printf("daligner");
            if (VON)
              printf(" -v");
            if (BON)
              printf(" -b");
            if (DON)
              printf(" -d");
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
            if (EREL != .7)
              printf(" -e%g",EREL);
            if (MREL != .55)
              printf(" -m%g",MREL);
            if (LINT != 1000)
              printf(" -l%d",LINT);
            if (SINT != 100)
              printf(" -s%d",SINT);
            if (usepath)
              printf(" %s/%s.%d",pwd,root,i);
            else
              printf(" %s.%d",root,i);
            hgh = (i*j)/bits + 1;
            for (k = low; k < hgh; k++)
              if (usepath)
                printf(" %s/%s.%d",pwd,root,k);
              else
                printf(" %s.%d",root,k);
            printf("\n");
            low = hgh;
          }
      }

    //  ... and then all the initial sort & merge jobs for each block pair

    printf("# Initial sort jobs (%d)\n", lblock*lblock - (fblock-1)*(fblock-1) );

    for (i = 1; i <= lblock; i++)
      for (j = (i < fblock ? fblock : 1); j <= lblock; j++)
        { printf("LAsort");
          if (VON)
            printf(" -v");
          for (k = 0; k < NTHREADS; k++)
            { printf(" %s.%d.%s.%d.C%d",root,i,root,j,k);
              printf(" %s.%d.%s.%d.N%d",root,i,root,j,k);
            }
          printf(" && LAmerge");
          if (VON)
            printf(" -v");
          if (lblock == 1)
            printf(" %s.%d",root,i);
          else if (i < fblock)
            printf(" %s.L1.%d.%d",root,i,(j-fblock)+1);
          else
            printf(" %s.L1.%d.%d",root,i,j);
          for (k = 0; k < NTHREADS; k++)
            { printf(" %s.%d.%s.%d.C%d.S",root,i,root,j,k);
              printf(" %s.%d.%s.%d.N%d.S",root,i,root,j,k);
            }
          printf(" && rm");
          for (k = 0; k < NTHREADS; k++)
            { printf(" %s.%d.%s.%d.C%d.S.las",root,i,root,j,k);
              printf(" %s.%d.%s.%d.N%d.S.las",root,i,root,j,k);
            }
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

              if (dnt >= 1)
                { int last;

                  last = (dnt == 1 || i == level);
                  printf("# Level %d jobs (%d)\n",i,bits*((lblock-fblock)+1) + dits*(fblock-1));
                  for (j = 1; j < fblock; j++)
                    { if (last)
                        printf("mv %s.%d.las %s.L%d.%d.0.las && ",root,j,root,i,j);
                      low = 1;
                      for (p = 1; p <= dits; p++)
                        { hgh = (dnt*p)/dits;
                          printf("LAmerge");
                          if (VON)
                            printf(" -v");
                          if (last)
                            printf(" %s.%d %s.L%d.%d.0",root,j,root,i,j);
                          else
                            printf(" %s.L%d.%d.%d",root,i+1,j,p);
                          for (k = low; k <= hgh; k++)
                            printf(" %s.L%d.%d.%d",root,i,j,k);
                          printf(" && rm");
                          if (last)
                            printf(" %s.L%d.%d.0.las",root,i,j);
                          for (k = low; k <= hgh; k++)
                            printf(" %s.L%d.%d.%d.las",root,i,j,k);
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
                      printf("LAmerge");
                      if (VON)
                        printf(" -v");
                      if (i == level)
                        printf(" %s.%d",root,j);
                      else
                        printf(" %s.L%d.%d.%d",root,i+1,j,p);
                      for (k = low; k <= hgh; k++)
                        printf(" %s.L%d.%d.%d",root,i,j,k);
                      printf(" && rm");
                      for (k = low; k <= hgh; k++)
                        printf(" %s.L%d.%d.%d.las",root,i,j,k);
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
