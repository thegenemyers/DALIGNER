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
  int     tspace, small;
  uint8  *tbuffer = NULL;
  uint16 *sbuffer = NULL;

  (void) argv;

  //  Process arguments

  if (argc > 1)
    { fprintf(stderr,"Usage: LAa2b <(ascii) >(binary)\n");
      exit (1);
    }

  if (fread(&code,sizeof(char),1,stdin) == 0)
    code = 0;

  while (code == '@' || code == '+' || code == '%')
    { fread(&which,sizeof(char),1,stdin);
      fread(&total,sizeof(int64),1,stdin);
      printf("%c %c %lld\n",code,which,total);
      if (code == '@')
        { tbuffer = (uint8 *) malloc(2*total*sizeof(uint16));
          sbuffer = (uint16 *) tbuffer;
        }

      if (fread(&code,sizeof(char),1,stdin) == 0)
       code = 0;
    }

  small = 0;
  if (tbuffer != NULL && code != 0)
    { if (code != 'X')
        { fprintf(stderr,"LAb2a: .las dump has traces but no X-info\n");
          exit (1);
        }
      fread(&tspace,sizeof(int),1,stdin);
      small = (tspace <= TRACE_XOVR && tspace != 0);
      printf("X %d\n",tspace);
    }

  while (code != 0)                        //  For each data item do
    { switch (code)
      { case 'P':                         //  Alignment pair
          fread(&aread,sizeof(int),1,stdin);
          fread(&bread,sizeof(int),1,stdin);
          fread(&orient,sizeof(char),1,stdin);
          fread(&chain,sizeof(char),1,stdin);
          printf("%c %d %d %c %c\n",code,aread,bread,orient,chain);
          break;
        case 'L':                         //  Read lengths
          scanf(" %d %d",&alen,&blen);
          fread(&len,sizeof(int),1,stdin);
          fread(&blen,sizeof(int),1,stdin);
          printf("%c %d %d\n",code,alen,blen);
          break;
        case 'C':                         //  Coordinate intervals
          fread(&ab,sizeof(int),1,stdin);
          fread(&ae,sizeof(int),1,stdin);
          fread(&bb,sizeof(int),1,stdin);
          fread(&be,sizeof(int),1,stdin);
          printf("%c %d %d %d %d\n",code,ab,ae,bb,be);
          break;
        case 'D':                         //  Differences
          fread(&diffs,sizeof(int),1,stdin);
          printf("%c %d\n",code,diffs);
          break;
        case 'T':                         //  Mask
          if (tbuffer == NULL)
            { fprintf(stderr,"LAb2a: .las dump has traces but no @ T-info\n");
              exit (1);
            }
          fread(&len,sizeof(int),1,stdin);
          printf("%c %d\n",code,len);
          len *= 2;
          if (small)
            { fread(tbuffer,sizeof(uint8),len,stdin);
              for (int i = 0; i < len; i += 2)
                printf(" %d %d\n",tbuffer[i],tbuffer[i+1]);
            }
          else
            { fread(sbuffer,sizeof(uint16),len,stdin);
              for (int i = 0; i < len; i += 2)
                printf(" %d %d\n",sbuffer[i],sbuffer[i+1]);
            }
      }

      if (fread(&code,sizeof(char),1,stdin) == 0)
       code = 0;
    }

  exit (0);
}
