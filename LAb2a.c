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
  unsigned char *tbuffer = NULL;

  (void) argv;

  //  Process arguments

  if (argc > 1)
    { fprintf(stderr,"Usage: LAa2b <(ascii) >(binary)\n");
      exit (1);
    }

  if (fread(&code,sizeof(char),1,stdin) == 0)
    code = 0;

  while (code == '@' || code == '+' || code == '%')
    { fread(&which,sizeof(char),1,stdout);
      fread(&total,sizeof(long long int),1,stdout);
      printf("%c %c %lld\n",code,which,total);
      if (code == '@')
        tbuffer = malloc(2*total*sizeof(unsigned char));

      if (fread(&code,sizeof(char),1,stdin) == 0)
       code = 0;
    }

  while (code != 0)                        //  For each data item do
    { switch (code)
      { case 'P':                         //  Alignment pair
          fread(&aread,sizeof(int),1,stdout);
          fread(&bread,sizeof(int),1,stdout);
          fread(&orient,sizeof(char),1,stdout);
          fread(&chain,sizeof(char),1,stdout);
          printf("%c %d %d %c %c\n",code,aread,bread,orient,chain);
          break;
        case 'L':                         //  Read lengths
          scanf(" %d %d",&alen,&blen);
          fread(&len,sizeof(int),1,stdout);
          fread(&blen,sizeof(int),1,stdout);
          printf("%c %d %d\n",code,alen,blen);
          break;
        case 'C':                         //  Coordinate intervals
          fread(&ab,sizeof(int),1,stdout);
          fread(&ae,sizeof(int),1,stdout);
          fread(&bb,sizeof(int),1,stdout);
          fread(&be,sizeof(int),1,stdout);
          printf("%c %d %d %d %d\n",code,ab,ae,bb,be);
          break;
        case 'D':                         //  Differences
          fread(&diffs,sizeof(int),1,stdout);
          printf("%c %d\n",code,diffs);
          break;
        case 'T':                         //  Mask
          fread(&len,sizeof(int),1,stdout);
          printf("%c %d\n",code,len);
          len *= 2;
          fread(tbuffer,sizeof(unsigned char),len,stdout);
          for (int i = 0; i < len; i += 2)
            printf("   %d %d\n",tbuffer[i],tbuffer[i+1]);
      }

      if (fread(&code,sizeof(char),1,stdin) == 0)
       code = 0;
    }

  exit (0);
}
