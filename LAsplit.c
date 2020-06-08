/*******************************************************************************************
 *
 *  Split an OVL file arriving from the standard input into 'parts' equal sized .las-files
 *    <align>.1.las, <align>.2.las ... or according to a current partitioning of <path>
 *
 *  Author:  Gene Myers
 *  Date  :  June 2014
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

static char *Usage = "-v <target:las> (<parts:int> | <path:db|dam>) < <source>.las";

#define MEMORY   1000   //  How many megabytes for output buffer

int main(int argc, char *argv[])
{ char      *iblock, *oblock;
  FILE      *output;
  DAZZ_STUB *stub;
  int64      novl, bsize, ovlsize, ptrsize;
  int        parts, tspace, tbytes;
  char      *pwd, *root, *root2;

  int        VERBOSE;

  //  Process options

  { int i, j, k;
    int flags[128];

    ARG_INIT("LAsplit")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("v") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"    <target> is a template that must have a single %c-sign in it\n",
                       BLOCK_SYMBOL);
        fprintf(stderr,"    This symbol is replaced by numbers 1 to n = the number of parts\n");
        exit (1);
      }
  }

  { char *eptr;

    parts = strtol(argv[2],&eptr,10);
    if (*eptr != '\0')
      { pwd = PathTo(argv[2]);
        if (strcmp(argv[2]+(strlen(argv[2])-4),".dam") == 0)
          { root = Root(argv[2],".dam");
            stub = Read_DB_Stub(Catenate(pwd,"/",root,".dam"),DB_STUB_BLOCKS);
            parts = stub->nblocks;
          }
        else
          { root = Root(argv[2],".db");
            stub = Read_DB_Stub(Catenate(pwd,"/",root,".db"),DB_STUB_BLOCKS);
            parts = stub->nblocks;
          }
        free(pwd);
        free(root);
      }
    else
      { stub = NULL;
        if (parts <= 0)
          { fprintf(stderr,"%s: Number of parts is not positive\n",Prog_Name);
            exit (1);
          }
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

  pwd   = PathTo(argv[1]);
  root  = Root(argv[1],".las");

  root2 = index(root,BLOCK_SYMBOL);
  if (root2 == NULL)
    { fprintf(stderr,"%s: No %c-sign in source name '%s'\n",Prog_Name,BLOCK_SYMBOL,root);
      exit (1);
    }
  if (index(root2+1,BLOCK_SYMBOL) != NULL)
    { fprintf(stderr,"%s: Two or more occurences of %c-sign in source name '%s'\n",
                     Prog_Name,BLOCK_SYMBOL,root);
      exit (1);
    }
  *root2++ = '\0';

  if (fread(&novl,sizeof(int64),1,stdin) != 1)
    SYSTEM_READ_ERROR
  if (fread(&tspace,sizeof(int),1,stdin) != 1)
    SYSTEM_READ_ERROR
  if (tspace <= TRACE_XOVR && tspace != 0)
    tbytes = sizeof(uint8);
  else
    tbytes = sizeof(uint16);

  if (VERBOSE)
    { printf("  Distributing %lld la\'s\n",novl);
      fflush(stdout);
    }

  { int      i;
    Overlap *w;
    int64    j, low, hgh, last;
    int64    tsize, povl;
    char    *iptr, *itop;
    char    *optr, *otop;

    iptr = iblock;
    itop = iblock + fread(iblock,1,bsize,stdin);

    hgh = 0;
    for (i = 0; i < parts; i++)
      { output = Fopen(Catenate(pwd,"/",Numbered_Suffix(root,i+1,root2),".las"),"w");
        if (output == NULL)
          exit (1);

        low = hgh;
        if (stub != NULL)
          { last = stub->tblocks[i+1];
            hgh  = 0;
          }
        else
          { last = 0;
            hgh  = (novl*(i+1))/parts;
          }

        povl = 0;
        fwrite(&povl,sizeof(int64),1,output);
        fwrite(&tspace,sizeof(int),1,output);

        optr = oblock;
        otop = oblock + bsize;

        for (j = low; j < novl; j++)
          { if (iptr + ovlsize > itop)
              { int64 remains = itop-iptr;
                if (remains > 0)
                  memmove(iblock,iptr,remains);
                iptr  = iblock;
                itop  = iblock + remains;
                itop += fread(itop,1,bsize-remains,stdin);
              }

            w = (Overlap *) (iptr-ptrsize);
            if (stub == NULL)
              { if (j >= hgh && w->aread > last)
                  break;
                last = w->aread;
              }
            else
              { if (w->aread >= last)
                  break;
              }

            tsize = w->path.tlen*tbytes;
            if (optr + ovlsize + tsize > otop)
              { fwrite(oblock,1,optr-oblock,output);
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
                itop += fread(itop,1,bsize-remains,stdin);
              }
	    memmove(optr,iptr,tsize);
            optr += tsize;
            iptr += tsize;
          }
        hgh = j;

        if (optr > oblock)
          fwrite(oblock,1,optr-oblock,output);

        rewind(output);
        povl = hgh-low;
        fwrite(&povl,sizeof(int64),1,output);

        if (VERBOSE)
          { printf("  Split off %s: %lld la\'s\n",Numbered_Suffix(root,i+1,root2),povl);
            fflush(stdout);
          }

        fclose(output);
      }
  }

  free(pwd);
  free(root);
  Free_DB_Stub(stub);
  free(iblock-ptrsize);
  free(oblock);

  exit (0);
}
