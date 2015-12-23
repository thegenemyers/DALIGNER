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

static char *Usage = "<align:las> (<parts:int> | <path:db|dam>) < <source>.las";

#define MEMORY   1000   //  How many megabytes for output buffer

int main(int argc, char *argv[])
{ char     *iblock, *oblock;
  FILE     *output, *dbvis;
  int64     novl, bsize, ovlsize, ptrsize;
  int       parts, tspace, tbytes;
  int       olast, blast;
  char     *root, *pwd;

  Prog_Name = Strdup("LAsplit","");

  if (argc != 3)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
      exit (1);
    }

  { char *eptr;
    int   nfiles, cutoff, all;
    int64 size;
    char  buffer[2*MAX_NAME+100];

    parts = strtol(argv[2],&eptr,10);
    if (*eptr != '\0')
      { pwd   = PathTo(argv[2]);
        if (strcmp(argv[2]+(strlen(argv[2])-4),".dam") == 0)
          root  = Root(argv[2],".dam");
        else
          root  = Root(argv[2],".db");
        dbvis = fopen(Catenate(pwd,"/",root,".dam"),"r");
        if (dbvis == NULL)
          { dbvis = fopen(Catenate(pwd,"/",root,".db"),"r");
            if (dbvis == NULL)
              { fprintf(stderr,"%s: Second argument '%s' is not an integer or a DB\n",
                               Prog_Name,argv[2]);
                exit (1);
              }
          }
        free(pwd);
        free(root);

        if (fscanf(dbvis,DB_NFILE,&nfiles) != 1)
          SYSTEM_ERROR
        while (nfiles-- > 0)
          if (fgets(buffer,2*MAX_NAME+100,dbvis) == NULL)
            SYSTEM_ERROR
        parts = 0;
        if (fscanf(dbvis,DB_NBLOCK,&parts) != 1)
          { fprintf(stderr,"%s: DB %s has not been partitioned\n",Prog_Name,argv[2]);
            exit (1);
          }
        if (fscanf(dbvis,DB_PARAMS,&size,&cutoff,&all) != 3)
          SYSTEM_ERROR
        if (fscanf(dbvis,DB_BDATA,&olast,&blast) != 2)
          SYSTEM_ERROR
      }
    else
      { dbvis = NULL;
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

  if (fread(&novl,sizeof(int64),1,stdin) != 1)
    SYSTEM_ERROR
  if (fread(&tspace,sizeof(int),1,stdin) != 1)
    SYSTEM_ERROR
  if (tspace <= TRACE_XOVR)
    tbytes = sizeof(uint8);
  else
    tbytes = sizeof(uint16);

  { int      i, j;
    Overlap *w;
    int      low, hgh, last;
    int64    tsize, povl;
    char    *iptr, *itop;
    char    *optr, *otop;

    iptr = iblock;
    itop = iblock + fread(iblock,1,bsize,stdin);

    hgh = 0;
    for (i = 0; i < parts; i++)
      { output = Fopen(Catenate(pwd,"/",root,Numbered_Suffix(".",i+1,".las")),"w");
        if (output == NULL)
          exit (1);

        low = hgh;
        if (dbvis != NULL)
          { if (fscanf(dbvis,DB_BDATA,&olast,&blast) != 2)
              SYSTEM_ERROR
            last = blast-1;
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
                  memcpy(iblock,iptr,remains);
                iptr  = iblock;
                itop  = iblock + remains;
                itop += fread(itop,1,bsize-remains,stdin);
              }

            w = (Overlap *) (iptr-ptrsize);
            if (dbvis == NULL)
              { if (j >= hgh && w->aread > last)
                  break;
                last = w->aread;
              }
            else
              { if (w->aread > last)
                  break;
              }

            tsize = w->path.tlen*tbytes;
            if (optr + ovlsize + tsize > otop)
              { fwrite(oblock,1,optr-oblock,output);
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
                itop += fread(itop,1,bsize-remains,stdin);
              }
	    memcpy(optr,iptr,tsize);
            optr += tsize;
            iptr += tsize;
          }
        hgh = j;

        if (optr > oblock)
          fwrite(oblock,1,optr-oblock,output);

        rewind(output);
        povl = hgh-low;
        fwrite(&povl,sizeof(int64),1,output);

        fclose(output);
      }
  }

  free(pwd);
  free(root);
  free(iblock-ptrsize);
  free(oblock);

  exit (0);
}
