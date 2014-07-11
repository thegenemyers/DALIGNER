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

static char *Usage = "<source:las> > <target>.las";

#define MEMORY   1000   //  How many megabytes for output buffer

int main(int argc, char *argv[])
{ void     *iblock, *oblock;
  FILE     *input;
  int64     novl, bsize, ovlsize, ptrsize;
  int       tspace, tbytes;
  char     *pwd, *root;

  Prog_Name = Strdup("LAcat","");

  if (argc <= 1)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
      exit (1);
    }

  ptrsize = sizeof(void *);
  ovlsize = sizeof(Overlap) - ptrsize;
  bsize   = MEMORY * 1000000ll;
  oblock  = Malloc(bsize,"Allocating output block");
  iblock  = Malloc(bsize,"Allocating input block");
  if (oblock == NULL || iblock == NULL)
    exit (1);

  pwd    = PathTo(argv[1]);
  root   = Root(argv[1],".las");

  { int64    povl;
    int      i, mspace;

    novl   = 0;
    tspace = 0;
    for (i = 0; 1; i++)
      { char *name = Catenate(pwd,"/",root,Numbered_Suffix(".",i+1,".las"));
        if ((input = fopen(name,"r")) == NULL) break;

        fread(&povl,sizeof(int64),1,input);
        novl += povl;
        fread(&mspace,sizeof(int),1,input);
        if (i == 0)
          { tspace = mspace;
            if (tspace <= TRACE_XOVR)
              tbytes = sizeof(uint8);
            else
              tbytes = sizeof(uint16);
          }
        else if (tspace != mspace)
          { fprintf(stderr,"%s: PT-point spacing conflict (%d vs %d)\n",Prog_Name,tspace,mspace);
            exit (1);
          }

        fclose(input);
      }
    fwrite(&novl,sizeof(int64),1,stdout);
    fwrite(&tspace,sizeof(int32),1,stdout);
  }

  { int      i, j;
    Overlap *w;
    int64    tsize, povl;
    int      mspace;
    void    *iptr, *itop;
    void    *optr, *otop;

    optr = oblock;
    otop = oblock + bsize;

    for (i = 0; 1; i++)
      { char *name = Catenate(pwd,"/",root,Numbered_Suffix(".",i+1,".las"));
        if ((input = fopen(name,"r")) == NULL) break;

        fread(&povl,sizeof(int64),1,input);
        fread(&mspace,sizeof(int),1,input);

        iptr = iblock;
        itop = iblock + fread(iblock,1,bsize,input);

        for (j = 0; j < povl; j++)
          { if (iptr + ovlsize > itop)
              { int64 remains = itop-iptr;
                if (remains > 0)
                  memcpy(iblock,iptr,remains);
                iptr  = iblock;
                itop  = iblock + remains;
                itop += fread(itop,1,bsize-remains,input);
              }

            w = (Overlap *) (iptr - ptrsize);
            tsize = w->path.tlen*tbytes;

            if (optr + ovlsize + tsize > otop)
              { fwrite(oblock,1,optr-oblock,stdout);
                optr = oblock;
              }

            memcpy(optr,iptr,ovlsize);
            optr += ovlsize;
            iptr += ovlsize;

            if (iptr + tsize > itop)
              { int64 remains = itop-iptr;
                if (remains > 0)
                  memcpy(iblock,iptr,remains);
                iptr  = iblock;
                itop  = iblock + remains;
                itop += fread(itop,1,bsize-remains,input);
              }
            
            memcpy(optr,iptr,tsize);
            optr += tsize;
            iptr += tsize;
          }

        fclose(input);
      }

    if (optr > oblock)
      fwrite(oblock,1,optr-oblock,stdout);
  }

  free(pwd);
  free(root);
  free(oblock);
  free(iblock);

  exit (0);
}
