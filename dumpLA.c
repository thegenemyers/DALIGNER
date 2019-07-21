#include <stdlib.h>
#include <stdio.h>

#include "DB.h"
#include "align.h"

int main(int argc, char *argv[])
{ char    code, which;
  int64   total;
  int     aread, bread;
  char    orient, chain;
  int     alen, blen;
  int     ab, ae, bb, be;
  int     diffs;
  int     len;
  int     tspace, small, tbytes;
  uint8  *tbuffer = NULL;
  uint16 *sbuffer;
  int     hasNext, haveC, haveT, haveD;
  int     i;
  int64   novls;
  Overlap _ovl, *ovl = &_ovl;
  Path   *path = & (ovl->path);
  FILE   *output;
  char   *pwd, *root;

  //  Process arguments

  if (argc != 2)
    { fprintf(stderr,"Usage: DumpLA <align:las> < (ascii dump)\n");
      exit (1);
    }

  pwd = PathTo(argv[1]);
  root = Root(argv[1],".las");
  if ((output = fopen(Catenate(pwd,"/",root,".las"),"w")) == NULL)
    { fprintf(stderr,"DumpLA: Cannot open %s for writing\n",argv[1]);
      exit (1);
    }
  free(root);
  free(pwd);

  hasNext = 0;
  while (scanf(" %c",&code) == 1)       //  Header lines
    if (code == '@' || code == '+' || code == '%')
      { scanf(" %c %lld",&which,&total);
        if (code == '@' && which == 'T')
          { tbuffer = (uint8 *) malloc(2*total*sizeof(uint16));
            sbuffer = (uint16 *) tbuffer;
          }
      }
    else
      { if (tbuffer == NULL)
          { fprintf(stderr,"DumpLA: .las dump must contain trace header lines\n");
            exit (1);
          }
        if (code != 'X')
          { fprintf(stderr,"DumpLA: .las dump must have an X-line after header\n");
            exit (1);
          }
        scanf(" %d",&tspace);
        if (tspace <= TRACE_XOVR && tspace != 0)
          { small = 1;
            tbytes = 1;
          }
        else
          { small = 0;
            tbytes = 2;
          }
        if (scanf(" %c",&code) == 1)
          { if (code != 'P')
              { fprintf(stderr,"DumpLA: .las dump data must being with a P-line\n");
                exit (1);
              }
            hasNext = 1;
          }
        break;
      }

  novls = 0;
  fwrite(&novls,sizeof(int64),1,output);
  fwrite(&tspace,sizeof(int),1,output);

  while (hasNext)       //  For each data line do
    { scanf(" %d %d %c %c",&aread,&bread,&orient,&chain);
      haveC = haveT = haveD = 0;
      hasNext = 0;
      while (scanf(" %c",&code) == 1)       //  For each data line do
        if (code == 'P')
          { hasNext = 1;
            break;
          }
        else
          switch (code)
          { case 'L':                         //  Read lengths
              scanf(" %d %d",&alen,&blen);
              break;
            case 'C':                         //  Coordinate intervals
              scanf(" %d %d %d %d",&ab,&ae,&bb,&be);
              haveC = 1;
              break;
            case 'D':                         //  Differences
              scanf(" %d",&diffs);
              haveD = 1;
              break;
            case 'T':                         //  Mask
              haveT = 1;
              scanf(" %d",&len);
              len *= 2;
              if (small)
                { for (int i = 0; i < len; i += 2)
                    scanf(" %hhd %hhd",tbuffer+i,tbuffer+(i+1));
                }
              else
                { for (int i = 0; i < len; i += 2)
                    scanf(" %hd %hd",sbuffer+i,sbuffer+(i+1));
                }
              break;
            default:
              fprintf(stderr,"DumpLA: Unrecognized line type '%c'\n",code);
              exit (1);
          }
      if (!haveC)
        { fprintf(stderr,"DumpLA: Alignment record does not have a C-line\n");
            exit (1);
        }
      if (!haveT)
        { fprintf(stderr,"DumpLA: Alignment record does not have a T-line\n");
            exit (1);
        }
      if (!haveD)
        { diffs = 0;
          for (i = 0; i < len; i += 2)
            diffs += tbuffer[i];
        }

      novls += 1;
      ovl->aread = aread-1;
      ovl->bread = bread-1;
      ovl->flags = 0;
      if (orient == 'c')
        ovl->flags |= COMP_FLAG;
      if (chain == '-')
        ovl->flags |= NEXT_FLAG;
      else if (chain == '>')
        ovl->flags |= BEST_FLAG;
      else if (chain == '+')
        ovl->flags |= START_FLAG;
      path->abpos = ab;
      path->aepos = ae;
      path->bbpos = bb;
      path->bepos = be;
      path->diffs = diffs;
      path->tlen  = len;
      path->trace = (void *) tbuffer;

      Write_Overlap(output,ovl,tbytes);
    }

  rewind(output);
  fwrite(&novls,sizeof(int64),1,output);

  fclose(output);

  exit (0);
}
