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
  { "[-vba] [-t<int>] [-w<int(6)>] [-l<int(1000)>] [-s<int(100)]",
    "       [-M<int>] [-B<int(4)>] [-D<int( 250)>] [-m<track>]+",
    "     ( [-k<int(14)>] [-h<int(35)>] [-e<double(.70)>] [-AI] [-H<int>] |",
    "       [-k<int(20)>] [-h<int(50)>] [-e<double(.85)>]  <ref:db|dam>   )",
    "       <reads:db|dam> [<first:int>[-<last:int>]"
  };

  //  Command Options

static int    DUNIT, BUNIT;
static int    VON, BON, AON, ION, CON;
static int    WINT, TINT, HGAP, HINT, KINT, SINT, LINT, MINT;
static double EREL;
static int    MMAX, MTOP;
static char **MASK;

static int power(int base, int exp)
{ int i, pow;

  pow = 1;
  for (i = 0; i < exp; i++)
    pow *= base;
  return (pow);
}

#define LSF_ALIGN "bsub -q medium -n 4 -o DALIGNER.out -e DALIGNER.err -R span[hosts=1] -J align#%d"
#define LSF_SORT  "bsub -q short -n 12 -o SORT.DAL.out -e SORT.DAL.err -R span[hosts=1] -J sort#%d"
#define LSF_MERGE \
          "bsub -q short -n 12 -o MERGE%d.DAL.out -e MERGE%d.DAL.err -R span[hosts=1] -J merge#%d"

void daligner_script(int argc, char *argv[])
{ int   nblocks;
  int   usepath;
  int   useblock;
  int   fblock, lblock;
#ifdef LSF
  int   jobid;
#endif

  char *pwd, *root;

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
    if (fscanf(dbvis,"blocks = %d\n",&nblocks) != 1 || nblocks == 1)
      { useblock = 0;
        nblocks  = 1;
      }

    usepath = (strcmp(pwd,".") != 0);
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
        useblock = 1;
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
      { file = fopen(Catenate(pwd,"/",root,Numbered_Suffix(".",fblock-1,".las")),"r");
        if (file == NULL)
          { if (usepath)
              fprintf(stderr,"%s: File %s/%s.%d.las should already be present!\n",
                             Prog_Name,pwd,root,fblock-1);
            else
              fprintf(stderr,"%s: File %s.%d.las should already be present!\n",
                             Prog_Name,root,fblock-1);
            exit (1);
          }
        else
          fclose(file);
      }
    if (useblock)
      file = fopen(Catenate(pwd,"/",root,Numbered_Suffix(".",fblock,".las")),"r");
    else
      file = fopen(Catenate(pwd,"/",root,".las"),"r");
    if (file != NULL)
      { if (usepath)
          if (useblock)
            fprintf(stderr,"%s: File %s/%s.%d.las should not yet exist!\n",
                           Prog_Name,pwd,root,fblock);
          else
            fprintf(stderr,"%s: File %s/%s.las should not yet exist!\n",Prog_Name,pwd,root);
        else
          if (useblock)
            fprintf(stderr,"%s: File %s.%d.las should not yet exist!\n",Prog_Name,root,fblock);
          else
            fprintf(stderr,"%s: File %s.las should not yet exist!\n",Prog_Name,root);
        exit (1);
      }
  }

  { int level, njobs;
    int i, j, k;

    //  Produce all necessary daligner jobs ...

    njobs = 0;
    for (i = fblock; i <= lblock; i++)
      njobs += (i-1)/BUNIT+1;

    printf("# Daligner jobs (%d)\n",njobs);

#ifdef LSF
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; i++)
      { int bits;
        int low, hgh;

        bits = (i-1)/BUNIT+1;
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
            if (EREL > 0.)
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
          printf(LSF_SORT,jobid++);
          printf(" \"");
#endif
          printf("LAsort");
          if (VON)
            printf(" -v");
          if (CON)
            printf(" -a");
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
          if (CON)
            printf(" -a");
          if (lblock == 1)
            { if (usepath)
                printf(" %s/%s",pwd,root);
              else
                printf(" %s",root);
              if (useblock)
                printf(".1");
            }
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
#ifdef LSF
          printf("\"");
#endif
          printf("\n");
        }

    //  Check .las files (optional)

    printf("# Check all level 1 .las files (optional but recommended)\n");

    for (i = 1; i <= lblock; i++)
      for (j = (i < fblock ? fblock : 1); j <= lblock; j++)
        { printf("LAcheck -vS");
          if (usepath)
            printf(" %s/%s",pwd,root);
          else
            printf(" %s",root);
          if (lblock == 1)
            { if (usepath)
                printf(" %s/%s",pwd,root);
              else
                printf(" %s",root);
              if (useblock)
                printf(".1");
            }
          else if (i < fblock)
            printf(" L1.%d.%d",i,(j-fblock)+1);
          else
            printf(" L1.%d.%d",i,j);
          printf("\n");
        }

    //  Clean up (optional)

    printf("# Remove initial .las files (optional)\n");

    for (i = 1; i <= lblock; i++)
      for (j = (i < fblock ? fblock : 1); j <= lblock; j++)
        { printf("rm");
          for (k = 0; k < NTHREADS; k++)
            if (useblock)
              { printf(" %s.%d.%s.%d.C%d.las",root,i,root,j,k);
                printf(" %s.%d.%s.%d.N%d.las",root,i,root,j,k);
              }
            else
              { printf(" %s.%s.C%d.las",root,root,k);
                printf(" %s.%s.N%d.las",root,root,k);
              }
          printf("\n");
          printf("rm");
          for (k = 0; k < NTHREADS; k++)
            if (useblock)
              { printf(" %s.%d.%s.%d.C%d.S.las",root,i,root,j,k);
                printf(" %s.%d.%s.%d.N%d.S.las",root,i,root,j,k);
              }
            else
              { printf(" %s.%s.C%d.S.las",root,root,k);
                printf(" %s.%s.N%d.S.las",root,root,k);
              }
          printf("\n");
        }

    //  Higher level merges (if lblock > 1)

    if (lblock > 1)
      { int pow, mway;

        //  Determine most balance mway for merging in ceil(log_mrg lblock) levels

        pow = 1;
        for (level = 0; pow < lblock; level++)
          pow *= DUNIT;

        for (mway = DUNIT; mway >= 3; mway--)
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
                      printf(LSF_MERGE,i,i,jobid++);
                      printf(" \"");
#endif
                      if (last)
                        { if (usepath)
                            printf("mv %s/%s.%d.las L%d.%d.0.las && ",pwd,root,j,i,j);
                          else
                            printf("mv %s.%d.las L%d.%d.0.las && ",root,j,i,j);
                        }
                      low = 1;
                      for (p = 1; p <= dits; p++)
                        { hgh = (dnt*p)/dits;
#ifdef LSF
                          if (p > 1)
                            { printf(LSF_MERGE,i,i,jobid++);
                              printf(" \"");
                            }
#endif
                          printf("LAmerge");
                          if (VON)
                            printf(" -v");
                          if (CON)
                            printf(" -a");
                          if (last)
                            if (usepath)
                              printf(" %s/%s.%d L%d.%d.0",pwd,root,j,i,j);
                            else
                              printf(" %s.%d L%d.%d.0",root,j,i,j);
                          else
                            printf(" L%d.%d.%d",i+1,j,p);
                          for (k = low; k <= hgh; k++)
                            printf(" L%d.%d.%d",i,j,k);
#ifdef LSF
                          printf("\"");
#endif
                          printf("\n");
                          low = hgh+1;
                        }
                    }
                }
              else
                printf("# Level %d jobs (%d)\n",i,bits*((lblock-fblock)+1));

              //  New block merges

              for (j = fblock; j <= lblock; j++) 
                { low = 1;
                  for (p = 1; p <= bits; p++)
                    { hgh = (cnt*p)/bits;
#ifdef LSF
                      printf(LSF_MERGE,i,i,jobid++);
                      printf(" \"");
#endif
                      printf("LAmerge");
                      if (VON)
                        printf(" -v");
                      if (CON)
                        printf(" -a");
                      if (i == level)
                        if (usepath)
                          printf(" %s/%s.%d",pwd,root,j);
                        else
                          printf(" %s.%d",root,j);
                      else
                        printf(" L%d.%d.%d",i+1,j,p);
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

              if (dnt >= 1)
                { int last;

                  last = (dnt == 1 || i == level);
                  for (j = 1; j < fblock; j++)
                    { low = 1;
                      for (p = 1; p <= dits; p++)
                        { hgh = (dnt*p)/dits;
                          printf("LAcheck -vS");
                          if (usepath)
                            printf(" %s/%s",pwd,root);
                          else
                            printf(" %s",root);
                          if (last)
                            if (usepath)
                              printf(" %s/%s.%d",pwd,root,j);
                            else
                              printf(" %s.%d",root,j);
                          else
                            printf(" L%d.%d.%d",i+1,j,p);
                          printf("\n");
                          low = hgh+1;
                        }
                    }
                }

              for (j = fblock; j <= lblock; j++) 
                { low = 1;
                  for (p = 1; p <= bits; p++)
                    { hgh = (cnt*p)/bits;
                      printf("LAcheck -vS");
                      if (usepath)
                        printf(" %s/%s",pwd,root);
                      else
                        printf(" %s",root);
                      if (i == level)
                        if (usepath)
                          printf(" %s/%s.%d",pwd,root,j);
                        else
                          printf(" %s.%d",root,j);
                      else
                        printf(" L%d.%d.%d",i+1,j,p);
                      printf("\n");
                      low = hgh+1;
                    }
                }

              //  Cleanup (optional)

              printf("# Remove level %d .las files (optional)\n",i);

              if (dnt >= 1)
                { int last;

                  last = (dnt == 1 || i == level);
                  for (j = 1; j < fblock; j++)
                    { low = 1;
                      for (p = 1; p <= dits; p++)
                        { hgh = (dnt*p)/dits;
                          printf("rm");
                          if (last)
                            printf(" L%d.%d.0.las",i,j);
                          for (k = low; k <= hgh; k++)
                            printf(" L%d.%d.%d.las",i,j,k);
                          printf("\n");
                          low = hgh+1;
                        }
                    }
                }

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

              if (dnt >= 1)
                { if (dnt > 1)
                    dnt = dits;
                  else
                    dnt = 0;
                }
              cnt  = bits;
            }
        }
    }
  }

  free(root);
  free(pwd);
}

/*********************************************************************************************\
 *
 *  Produce a script to compute overlaps for all block pairs between two DBs, and then sort
 *    and merge *    them into as many .las files as their are blocks of the 1st DB.
 *
 *  Author:  Gene Myers
 *  Date  :  December 31, 2014
 *
 *********************************************************************************************/
 
#define LSF_MALIGN "bsub -q medium -n 4 -o MAPALL.out -e MAPALL.err -R span[hosts=1] -J align#%d"
#define LSF_MSORT  "bsub -q short -n 12 -o SORT.ALL.out -e SORT.ALL.err -R span[hosts=1] -J sort#%d"
#define LSF_MMERGE \
            "bsub -q short -n 12 -o MERGE%d.ALL.out -e MERGE%d.ALL.err -R span[hosts=1] -J merge#%d"

void mapper_script(int argc, char *argv[])
{ int   nblocks1, nblocks2;
  int   useblock1, useblock2;
  int   usepath1, usepath2;
  int   fblock, lblock;
#ifdef LSF
  int   jobid;
#endif

  char *pwd1, *root1;
  char *pwd2, *root2;

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
        if ((*eptr != '\0' && *eptr != '-') || eptr <= argv[3])
          { fprintf(stderr,"%s: final argument '%s' does not start with an integer\n",
                           Prog_Name,argv[3]);
            exit (1);
          }
        useblock2 = 1;
        if (*eptr == '-')
          { lblock = strtol(eptr+1,&fptr,10);
            if (*fptr != '\0' || fptr <= eptr+1)
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

    njobs = nblocks1 * ( (lblock-fblock)/BUNIT + 1);

    printf("# Daligner jobs (%d)\n",njobs);

#ifdef LSF
    jobid = 1;
#endif
    for (i = fblock; i <= lblock; i++)
      { int bits;
        int low, hgh;

        bits = (nblocks1-1)/BUNIT+1;
        low  = 1;
        for (j = 1; j <= bits; j++)
          {
#ifdef LSF
            printf(LSF_MALIGN,jobid++);
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
            if (EREL > 0.)
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
          printf(LSF_MSORT,jobid++);
          printf(" \"");
#endif
          printf("LAsort ");
          if (VON)
            printf("-v ");
          if (CON)
            printf("-a ");
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
            printf("-a ");
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
          pow *= DUNIT;

        for (mway = DUNIT; mway >= 3; mway--)
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
                      printf(LSF_MMERGE,i,i,jobid++);
                      printf(" \"");
#endif
                      printf("LAmerge ");
                      if (VON)
                        printf("-v ");
                      if (CON)
                        printf("-a ");
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

int main(int argc, char *argv[])
{ int    i, j, k;
  int    flags[128];
  char  *eptr;
  int    mapper;

  //  Process options and decide if its a overlap or mapper script

  ARG_INIT("HPC.daligner")

  KINT  = 0;
  HINT  = 0;
  HGAP  = 0;
  EREL  = 0.;

  BUNIT = 4;
  DUNIT = 250;
  TINT  = 0;
  WINT  = 6;
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
          ARG_FLAGS("vbaAI");
          break;
        case 'e':
          ARG_REAL(EREL)
          if (EREL < .7 || EREL >= 1.)
            { fprintf(stderr,"%s: Average correlation must be in [.7,1.) (%g)\n",Prog_Name,EREL);
              exit (1);
            }
          break;
        case 'h':
          ARG_POSITIVE(HINT,"Hit threshold (in bp.s)")
          break;
        case 'k':
          ARG_POSITIVE(KINT,"K-mer length")
          break;
        case 'l':
          ARG_POSITIVE(LINT,"Minimum ovlerap length")
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
        case 's':
          ARG_POSITIVE(SINT,"Trace spacing")
          break;
        case 't':
          ARG_POSITIVE(TINT,"Tuple suppression frequency")
          break;
        case 'w':
          ARG_POSITIVE(WINT,"Log of bin width")
          break;
        case 'B':
          ARG_NON_NEGATIVE(BUNIT,"Blocks per command")
          break;
        case 'D':
          ARG_NON_NEGATIVE(DUNIT,"File per merge")
          if (DUNIT < 3)
            { fprintf(stderr,"%s: Files per merge must be at least 3 (%d)\n",
                             Prog_Name,DUNIT);
              exit (1);
            }
          break;
        case 'H':
          ARG_POSITIVE(HGAP,"HGAP threshold (in bp.s)")
          break;
        case 'M':
          ARG_NON_NEGATIVE(MINT,"Memory allocation (in Gb)")
          break;
      }
    else
      argv[j++] = argv[i];
  argc = j;

  VON = flags['v'];
  BON = flags['b'];
  AON = flags['A'];
  ION = flags['I'];
  CON = flags['a'];

  if (argc < 2 || argc > 4)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
      fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
      fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
      fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[3]);
      fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[4]);
      exit (1);
    }

  if (argc == 2)
    mapper = 0;
  else if (argc == 4)
    mapper = 1;
  else
    { (void) strtol(argv[2],&eptr,10);
      if ((*eptr == '\0' || *eptr == '-') && eptr > argv[2])
        mapper = 0;
      else
        mapper = 1;
    }

  if (mapper)
    { if (AON || ION)
        { fprintf(stderr,"%s: Cannot use -A or -I options in a comparison script\n",Prog_Name);
          exit (1);
        }
      if (HGAP > 0)
        { fprintf(stderr,"%s: Cannot use -H option in a comparison script\n",Prog_Name);
          exit (1);
        }
      if (KINT <= 0)
        KINT = 20;
      if (HINT <= 0)
        HINT = 50;
      if (EREL <= 0.)
        EREL = .85;
    }
  else
    { if (KINT <= 0)
        KINT = 14;
      if (HINT <= 0)
        HINT = 35;
    }

  if (mapper)
    mapper_script(argc,argv);
  else
    daligner_script(argc,argv);

  exit (0);
}
