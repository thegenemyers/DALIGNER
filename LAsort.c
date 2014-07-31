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
 *  Load a file U.las of overlaps into memory, sort them all by A,B index,
 *    and then output the result to U.S.las
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

static char *Usage = "[-v] <align:las> ...";

#define MEMORY   1000   //  How many megabytes for output buffer

static void *IBLOCK;

static int SORT_OVL(const void *x, const void *y)
{ int64 l = *((int64 *) x);
  int64 r = *((int64 *) y);

  Overlap *ol, *or;
  int      al, ar;
  int      bl, br;
  int      pl, pr;

  ol = (Overlap *) (IBLOCK+l);
  al = ol->aread;
  bl = ol->bread;

  or = (Overlap *) (IBLOCK+r);
  ar = or->aread;
  br = or->bread;

  if (al != ar)
    return (al-ar);
  if (bl != br)
    return (bl-br);

  pl = ol->path.abpos;
  pr = or->path.abpos;

  return (pl-pr);
}

int main(int argc, char *argv[])
{ void     *iblock, *fblock;
  int64     isize,   osize;
  int64     ovlsize, ptrsize;
  int       tspace, tbytes;
  int       i;

  int       VERBOSE;
 
  //  Process options

  { int j, k;
    int flags[128];

    ARG_INIT("LAsort")

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
        exit (1);
      }

  }

  //  For each file do

  ptrsize = sizeof(void *);
  ovlsize = sizeof(Overlap) - ptrsize;
  isize   = 0;
  iblock  = NULL;
  osize   = MEMORY * 1000000ll;
  fblock  = Malloc(osize,"Allocating LAsort output block");

  for (i = 1; i < argc; i++)
    { int64    *perm;
      FILE     *input, *foutput;
      int64     novl;

      //  Read in the entire file and output header

      { int64  size;
        struct stat info;
        char  *pwd, *root, *name;

        pwd   = PathTo(argv[i]);
        root  = Root(argv[i],".las");
        name  = Catenate(pwd,"/",root,".las");
        input = Fopen(name,"r");
        if (input == NULL)
          exit (1);

        stat(name,&info);
        size = info.st_size;

        fread(&novl,sizeof(int64),1,input);
        fread(&tspace,sizeof(int),1,input);

        if (tspace <= TRACE_XOVR)
          tbytes = sizeof(uint8);
        else
          tbytes = sizeof(uint16);

        if (VERBOSE)
          { printf("  %s: ",root);
            Print_Number(novl,0,stdout);
            printf(" records ");
            Print_Number(size-novl*ovlsize,0,stdout);
            printf(" trace bytes\n");
            fflush(stdout);
          }

        foutput = Fopen(Catenate(pwd,"/",root,".S.las"),"w");
        if (foutput == NULL)
          exit (1);

        fwrite(&novl,sizeof(int64),1,foutput);
        fwrite(&tspace,sizeof(int),1,foutput);

        free(pwd);
        free(root);

        if (size > isize)
          { if (iblock == NULL)
              iblock = Malloc(size+ptrsize,"Allocating LAsort input block");
            else
              iblock = Realloc(iblock-ptrsize,size+ptrsize,"Allocating LAsort input block");
            if (iblock == NULL)
              exit (1);
            iblock += ptrsize;
            isize   = size;
          }
        fread(iblock,1,isize,input);
        fclose(input);
      }

      //  Set up unsorted permutation array

      perm = (int64 *) Malloc(sizeof(int64)*novl,"Allocating LAsort permutation vector");
      if (perm == NULL)
        exit (1);

      { int64 off;
        int   j;

        off = -ptrsize;
        for (j = 0; j < novl; j++)
          { perm[j] = off;
            off += ovlsize + ((Overlap *) (iblock+off))->path.tlen*tbytes;
          }
      }

      //  Sort permutation array of ptrs to records

      IBLOCK = iblock;
      qsort(perm,novl,sizeof(int64),SORT_OVL);

      //  Output the records in sorted order

      { int      j;
        Overlap *w;
        int64    tsize, span;
        void    *fptr, *ftop;

        fptr = fblock;
        ftop = fblock + osize;
        for (j = 0; j < novl; j++)
          { w = (Overlap *) (iblock+perm[j]);
            tsize = w->path.tlen*tbytes;
            span  = ovlsize + tsize;
            if (fptr + span > ftop)
              { fwrite(fblock,1,fptr-fblock,foutput);
                fptr = fblock;
              }
            memcpy(fptr,((void *) w)+ptrsize,ovlsize);
            fptr += ovlsize;
            memcpy(fptr,(void *) (w+1),tsize);
            fptr += tsize;
          }
        if (fptr > fblock)
          fwrite(fblock,1,fptr-fblock,foutput);
      }

      free(perm);
      fclose(foutput);
    }

  if (iblock != NULL)
    free(iblock - ptrsize);
  free(fblock);

  exit (0);
}
