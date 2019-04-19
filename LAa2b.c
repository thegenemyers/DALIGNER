#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{ char   code, which;
  long long int total;
  int    aread, bread;
  char   orient, chain;
  int    alen, blen;
  int    ab, ae, bb, be;
  int    diffs;
  int    len;
  unsigned char *tbuffer;

  (void) argv;

  //  Process arguments

  if (argc > 1)
    { fprintf(stderr,"Usage: LAa2b <(ascii) >(binary)\n");
      exit (1);
    }

  while (scanf(" %c",&code) == 1)       //  Header lines
    if (code == '@' || code == '+' || code == '%')
      { scanf(" %c %lld",&which,&total);
        fwrite(&code,sizeof(char),1,stdout);
        fwrite(&which,sizeof(char),1,stdout);
        fwrite(&total,sizeof(long long int),1,stdout);
        if (code == '@')
          tbuffer = malloc(2*total*sizeof(unsigned char));
      }
    else
      { ungetc(code,stdin);
        break;
      }

  while (scanf(" %c",&code) == 1)       //  For each data line do
    { fwrite(&code,sizeof(char),1,stdout);
      switch (code)
      { case 'P':                         //  Alignment pair
          scanf(" %d %d %c %c",&aread,&bread,&orient,&chain);
          fwrite(&aread,sizeof(int),1,stdout);
          fwrite(&bread,sizeof(int),1,stdout);
          fwrite(&orient,sizeof(char),1,stdout);
          fwrite(&chain,sizeof(char),1,stdout);
          break;
        case 'L':                         //  Read lengths
          scanf(" %d %d",&alen,&blen);
          fwrite(&len,sizeof(int),1,stdout);
          fwrite(&blen,sizeof(int),1,stdout);
          break;
        case 'C':                         //  Coordinate intervals
          scanf(" %d %d %d %d",&ab,&ae,&bb,&be);
          fwrite(&ab,sizeof(int),1,stdout);
          fwrite(&ae,sizeof(int),1,stdout);
          fwrite(&bb,sizeof(int),1,stdout);
          fwrite(&be,sizeof(int),1,stdout);
          break;
        case 'D':                         //  Differences
          scanf(" %d",&diffs);
          fwrite(&diffs,sizeof(int),1,stdout);
          break;
        case 'T':                         //  Mask
          scanf(" %d",&len);
          fwrite(&len,sizeof(int),1,stdout);
          len *= 2;
          for (int i = 0; i < len; i += 2)
            scanf(" %hhd %hhd",tbuffer+i,tbuffer+(i+1));
          fwrite(tbuffer,sizeof(unsigned char),len,stdout);
      }
    }

  exit (0);
}
