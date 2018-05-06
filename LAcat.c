/*******************************************************************************************
 *
 *  Merge together in index order, overlap files <XXX>.1.las, <XXX>.2.las, ... into a
 *    single overlap file and output to the standard output
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"
#include "align.h"

static char *Usage = "[-v] <source:las> ... > <target>.las";

#define MEMORY   1000         //  How many megabytes for output buffer

int main(int argc, char *argv[])
{ char     *iblock, *oblock;
  FILE     *input;
  int64     novl, bsize, ovlsize, ptrsize;
  int       tspace, tbytes;
  int       c;

  int       VERBOSE;

  //  Process options

  { int i, j, k;
    int flags[128];

    ARG_INIT("LAcat")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("v") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"    <source>'s may contain a template that is %c-sign optionally\n",
                        BLOCK_SYMBOL);
        fprintf(stderr,"      followed by an integer or integer range\n");
        exit (1);
      }
  }

  ptrsize = sizeof(void *);
  ovlsize = sizeof(Overlap) - ptrsize;
  bsize   = MEMORY * 1000000ll;
  oblock  = (char *) Malloc(bsize,"Allocating output block");
  iblock  = (char *) Malloc(bsize + ptrsize,"Allocating input block");
  if (oblock == NULL || iblock == NULL)
    exit (1);
  iblock += ptrsize;

  novl   = 0;
  tspace = -1;
  for (c = 1; c < argc; c++)
    { Block_Looper *parse;
      FILE *input;

      parse = Parse_Block_Arg(argv[c]);

      while ((input = Next_Block_Arg(parse)) != NULL)
        { int64 povl;
          int   mspace;

          if (fread(&povl,sizeof(int64),1,input) != 1)
            SYSTEM_READ_ERROR
          novl += povl;
          if (fread(&mspace,sizeof(int),1,input) != 1)
            SYSTEM_READ_ERROR
          if (tspace < 0)
            tspace = mspace;
          else if (tspace != mspace)
            { fprintf(stderr,"%s: trace-point spacing conflict between %s and earlier files",
                             Prog_Name,Block_Arg_Root(parse));
              fprintf(stderr," (%d vs %d)\n",tspace,mspace);
              exit (1);
            }

          fclose(input);
        }

      Free_Block_Arg(parse);
    }

  if (tspace <= TRACE_XOVR && tspace != 0)
    tbytes = sizeof(uint8);
  else
    tbytes = sizeof(uint16);
  if (fwrite(&novl,sizeof(int64),1,stdout) != 1)
    SYSTEM_READ_ERROR
  if (fwrite(&tspace,sizeof(int),1,stdout) != 1)
    SYSTEM_READ_ERROR

  { Block_Looper *parse;
    int      c, j;
    Overlap *w;
    int64    tsize, povl;
    int      mspace;
    char    *iptr, *itop;
    char    *optr, *otop;

    optr = oblock;
    otop = oblock + bsize;

    for (c = 1; c < argc; c++)
      { parse = Parse_Block_Arg(argv[c]);

        while ((input = Next_Block_Arg(parse)) != NULL)
          { if (fread(&povl,sizeof(int64),1,input) != 1)
              SYSTEM_READ_ERROR
            if (fread(&mspace,sizeof(int),1,input) != 1)
              SYSTEM_READ_ERROR

            if (VERBOSE)
              { fprintf(stderr,
                    "  Concatenating %s: %lld la\'s\n",Block_Arg_Root(parse),povl);
                fflush(stderr);
              }

            iptr = iblock;
            itop = iblock + fread(iblock,1,bsize,input);

            for (j = 0; j < povl; j++)
              { if (iptr + ovlsize > itop)
                  { int64 remains = itop-iptr;
                    if (remains > 0)
                      memmove(iblock,iptr,remains);
                    iptr  = iblock;
                    itop  = iblock + remains;
                    itop += fread(itop,1,bsize-remains,input);
                  }

                w = (Overlap *) (iptr - ptrsize);
                tsize = w->path.tlen*tbytes;

                if (optr + ovlsize + tsize > otop)
                  { if (fwrite(oblock,1,optr-oblock,stdout) != (size_t) (optr-oblock))
                      SYSTEM_READ_ERROR
                    optr = oblock;
                  }

                memmove(optr,iptr,ovlsize);
                optr += ovlsize;
                iptr += ovlsize;

                if (iptr + tsize > itop)
                  { int64 remains = itop-iptr;
                    if (remains > 0)
                      memmove(iblock,iptr,remains);
                    iptr  = iblock;
                    itop  = iblock + remains;
                    itop += fread(itop,1,bsize-remains,input);
                  }

                memmove(optr,iptr,tsize);
                optr += tsize;
                iptr += tsize;
              }

            fclose(input);
          }

        Free_Block_Arg(parse);
      }

    if (optr > oblock)
      { if (fwrite(oblock,1,optr-oblock,stdout) != (size_t) (optr-oblock))
          SYSTEM_READ_ERROR
      }
  }

  if (VERBOSE)
    { fprintf(stderr,"  Totalling %lld la\'s\n",novl);
      fflush(stderr);
    }

  free(oblock);
  free(iblock-ptrsize);

  exit (0);
}
