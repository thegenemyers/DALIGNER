/*******************************************************************************************
 *
 *  Check the structural integrity of .las files
 *
 *  Author:  Gene Myers
 *  Date  :  July 2014
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

static char *Usage = "[-vS] [-a:<db>] <align:las> ...";

#define MEMORY   1000   //  How many megabytes for output buffer

int main(int argc, char *argv[])
{ HITS_DB   _db,  *db  = &_db;
  HITS_READ *reads;
  int        VERBOSE;
  int        SORTED;
  int        RANGE;

  //  Process options

  { int i, j, k;
    int flags[128];

    ARG_INIT("LAcheck")

    RANGE = 0;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vS")
            break;
          case 'a':
            RANGE = 1;
            if (argv[i][2] != ':')
              { fprintf(stderr,"%s: Unrecognizable option %s\n",Prog_Name,argv[i]);
                exit (1);
              }
            if (Open_DB(argv[i]+3,db))
              exit (1);
            Trim_DB(db);
            reads = db->reads;
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    SORTED  = flags['S'];

    if (argc < 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  { void  *iblock;
    int64  bsize, ovlsize, ptrsize;
    int    i, j;

    //  Setup IO buffers

    ptrsize = sizeof(void *);
    ovlsize = sizeof(Overlap) - ptrsize;
    bsize   = MEMORY * 1000000ll;
    iblock  = Malloc(bsize,"Allocating input block");
    if (iblock == NULL)
      exit (1);

    //  For each file do

    for (i = 1; i < argc; i++)
      { char     *pwd, *root;
        FILE     *input;
        void     *iptr, *itop;
        Overlap   last;
        int64     novl;
        int       tspace, tbytes;

        //  Establish IO and (novl,tspace) header

        pwd    = PathTo(argv[i]);
        root   = Root(argv[i],".las");
        if ((input = Fopen(Catenate(pwd,"/",root,".las"),"r")) == NULL)
          goto error;

        fread(&novl,sizeof(int64),1,input);
        fread(&tspace,sizeof(int),1,input);
        if (novl < 0)
          { if (VERBOSE)
              fprintf(stderr,"  %s: Number of alignments < 0\n",root);
            goto error;
          }
        if (tspace < 0)
          { if (VERBOSE)
              fprintf(stderr,"  %s: Trace spacing < 0\n",root);
            goto error;
          }

        if (tspace <= TRACE_XOVR)
          tbytes = sizeof(uint8);
        else
          tbytes = sizeof(uint16);

        iptr = iblock;
        itop = iblock + fread(iblock,1,bsize,input);

        //  For each record in file do

        last.aread = -1;
        for (j = 0; j < novl; j++)
          { Overlap ovl;
            int     tsize;
            int     equal;

            //  Fetch next record

            if (iptr + ovlsize > itop)
              { int64 remains = itop-iptr;
                if (remains > 0)
                  memcpy(iblock,iptr,remains);
                iptr  = iblock;
                itop  = iblock + remains;
                itop += fread(itop,1,bsize-remains,input);
                if (iptr + ovlsize > itop)
                  { if (VERBOSE)
                      fprintf(stderr,"  %s: Too few alignment records\n",root);
                    goto error;
                  }
              }

            ovl   = *((Overlap *) (iptr - ptrsize));
            iptr += ovlsize;
            tsize = ovl.path.tlen*tbytes;

            if (iptr + tsize > itop)
              { int64 remains = itop-iptr;
                if (remains > 0)
                  memcpy(iblock,iptr,remains);
                iptr  = iblock;
                itop  = iblock + remains;
                itop += fread(itop,1,bsize-remains,input);
                if (iptr + tsize > itop)
                  { if (VERBOSE)
                      fprintf(stderr,"  %s: Too few alignment records\n",root);
                    goto error;
                  }
              }
            ovl.path.trace = iptr;
            iptr += tsize;

            //  Basic checks

            if (ovl.aread < 0 || ovl.bread < 0)
              { if (VERBOSE)
                  fprintf(stderr,"  %s: Read indices < 0\n",root);
                goto error;
              }
            if (ovl.path.abpos >= ovl.path.aepos || ovl.path.bbpos >= ovl.path.bepos ||
                ovl.path.aepos >  ovl.alen       || ovl.path.bepos >  ovl.blen)
              { if (VERBOSE)
                  fprintf(stderr,"  %s: Non-sense alignment intervals\n",root);
                goto error;
              }

            if (Check_Trace_Points(&ovl,tspace,VERBOSE,root))
              goto error;

            //  Range check is -a set

            if (RANGE)
              { if (ovl.aread >= db->nreads || ovl.bread >= db->nreads)
                  { if (VERBOSE)
                      fprintf(stderr,"  %s: Read indices out of range\n",root);
                    goto error;
                  }
                if (ovl.alen != reads[ovl.aread].end - reads[ovl.aread].beg ||
                    ovl.blen != reads[ovl.bread].end - reads[ovl.bread].beg)
                  { if (VERBOSE)
                      fprintf(stderr,"  %s: Read lengths corrupted\n",root);
                    goto error;
                  }
              }

            //  Duplicate check and sort check if -S set

            equal = 0;
            if (SORTED)
              { if (ovl.aread > last.aread) goto inorder;
                if (ovl.aread == last.aread)
                  { if (ovl.bread > last.bread) goto inorder;
                    if (ovl.bread == last.bread)
                      { if (COMP(ovl.flags) > COMP(last.flags)) goto inorder;
                        if (COMP(ovl.flags) == COMP(last.flags))
                          { if (ovl.path.abpos > last.path.abpos) goto inorder;
                            if (ovl.path.abpos == last.path.abpos)
                              { equal = 1;
                                goto inorder;
                              }
                          }
                      }
                  }
                if (VERBOSE)
                  fprintf(stderr,"  %s: Reads are not sorted\n",root);
                goto error;
              }
            else
              { if (ovl.aread == last.aread && ovl.bread == last.bread &&
                    COMP(ovl.flags) == COMP(last.flags) && ovl.path.abpos == last.path.abpos)
                  equal = 1;
              }
          inorder:
            if (equal)
              { if (ovl.path.aepos == last.path.aepos &&
                    ovl.path.bbpos == last.path.bbpos &&
                    ovl.path.bepos == last.path.bepos)
                  { if (VERBOSE)
                      fprintf(stderr,"  %s: Reads are not sorted\n",root);
                    goto error;
                  }
              }

            last = ovl;
          }

        //  File processing epilog: Check all data read and print OK if -v

        if (iptr < itop)
          { if (VERBOSE)
              fprintf(stderr,"  %s: Too many alignment records\n",root);
            goto error;
          }

        if (VERBOSE)
          { fprintf(stderr,"  %s: ",root);
            Print_Number(novl,0,stderr);
            fprintf(stderr," all OK\n");
          }

      error:
        fclose(input);
        free(pwd);
        free(root);
      }

    free(iblock);
  }

  if (RANGE)
    Close_DB(db);

  exit (0);
}
