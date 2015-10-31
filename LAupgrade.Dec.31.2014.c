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
 *  Convert older .las files so that the alen and blen fields are removed.
 *
 *  Author:  Gene Myers
 *  Date  :  Dec 2014
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

typedef struct
  { void     *trace;
    uint16    tlen;
    uint16    diffs;
    uint16    abpos, bbpos;
    uint16    aepos, bepos;
  } PathOld;

typedef struct {
  PathOld path;
  int     aread;
  int     bread;
  uint16  alen;
  uint16  blen;
  int     flags;
} OverlapOld;

static char *Usage = "<source:las> > <target>.las";

#define MEMORY   1000   //  How many megabytes for output buffer

int main(int argc, char *argv[])
{ char     *iblock, *oblock;
  FILE     *input;
  int64     novl, bsize, ovlsize, newsize, ptrsize;
  int       tspace, tbytes;
  char     *pwd, *root;

  Prog_Name = Strdup("Upgrade","");

  if (argc <= 1)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
      exit (1);
    }

  ptrsize = sizeof(void *);
  ovlsize = sizeof(OverlapOld) - ptrsize;
  newsize = sizeof(Overlap   ) - ptrsize;
  bsize   = MEMORY * 1000000ll;
  oblock  = (char *) Malloc(bsize,"Allocating output block");
  iblock  = (char *) Malloc(bsize + ptrsize,"Allocating input block");
  if (oblock == NULL || iblock == NULL)
    exit (1);
  iblock += ptrsize;

  pwd   = PathTo(argv[1]);
  root  = Root(argv[1],".las");
  input = Fopen(Catenate(pwd,"/",root,".las"),"r");
  if (input == NULL)
    exit (1);
  free(pwd);
  free(root);

  if (fread(&novl,sizeof(int64),1,input) != 1)
    SYSTEM_ERROR
  if (fread(&tspace,sizeof(int),1,input) != 1)
    SYSTEM_ERROR
  if (tspace <= TRACE_XOVR)
    tbytes = sizeof(uint8);
  else
    tbytes = sizeof(uint16);

  fwrite(&novl,sizeof(int64),1,stdout);
  fwrite(&tspace,sizeof(int),1,stdout);

  { int         j;
    OverlapOld *w;
    Overlap    *v;
    int64       tsize;
    char       *iptr, *itop;
    char       *optr, *otop;

    optr = oblock;
    otop = oblock + bsize;

    iptr = iblock;
    itop = iblock + fread(iblock,1,bsize,input);

    for (j = 0; j < novl; j++)
      { if (iptr + ovlsize > itop)
          { int64 remains = itop-iptr;
            if (remains > 0)
              memcpy(iblock,iptr,remains);
            iptr  = iblock;
            itop  = iblock + remains;
            itop += fread(itop,1,bsize-remains,input);
          }

        w = (OverlapOld *) (iptr - ptrsize);
        tsize = w->path.tlen*tbytes;

        if (optr + newsize + tsize > otop)
          { fwrite(oblock,1,optr-oblock,stdout);
            optr = oblock;
          }

        v = (Overlap *) (optr - ptrsize);
        v->path.abpos = w->path.abpos;
        v->path.bbpos = w->path.bbpos;
        v->path.aepos = w->path.aepos;
        v->path.bepos = w->path.bepos;
        v->path.diffs = w->path.diffs;
        v->path.tlen  = w->path.tlen;
        v->aread      = w->aread;
        v->bread      = w->bread;
        v->flags      = w->flags;

        optr += newsize;
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
    if (optr > oblock)
      fwrite(oblock,1,optr-oblock,stdout);
  }

  fclose(input);
  free(oblock);
  free(iblock-ptrsize);

  exit (0);
}
