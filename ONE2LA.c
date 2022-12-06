#include <stdlib.h>
#include <stdio.h>

#include "DB.h"
#include "align.h"
#include "ONElib.h"

static char *One_Schema =
  "P 3 dal\n"
  "D X 1 3 INT\n"             //  Data prolog: trace spacing
  "O P 2 3 INT 8 INT_LIST\n"  //  A-read and B-read list
                              //    All per B-read
  "D O 1 6 STRING\n"          //       orientation [+-]
  "D C 1 6 STRING\n"          //       chain directive [>+-.]
  "D A 1 8 INT_LIST\n"        //       (ab,ae)
  "D B 1 8 INT_LIST\n"        //       (be,be)
  "D L 2 3 INT 8 INT_LIST\n"  //       la and then each lb
  "D D 1 8 INT_LIST\n"        //       diff
                              //    One line per B-read
  "D T 1 8 INT_LIST\n"        //       trace segment length
  "D Q 1 8 INT_LIST\n";       //       trace segment diffs

int main(int argc, char *argv[])
{ int64    novls, Pmax, Tmax, psize;
  Overlap *ovls, *otop, *o;
  void    *trace, *ttop;
  int      aread;
  int      has[128];
  int      tspace, small, tbytes;

  OneFile   *file1;
  OneSchema *schema;
  int64     *list;
  char      *string;

  int t, i, j, k;

  //  Process arguments

  if (argc != 2)
    { fprintf(stderr,"Usage: ONE2LA <align:dal> > (.las)\n");
      exit (1);
    }

  { char *pwd, *root, *path;
    FILE *output;

    pwd = PathTo(argv[1]);
    root = Root(argv[1],".dal");
    path = Catenate(pwd,"/",root,".dal");
    if ((output = fopen(path,"r")) == NULL)
      { free(root);
        root = Root(argv[1],".1dal");
        path = Catenate(pwd,"/",root,".1dal");
        if ((output = fopen(path,"r")) == NULL)
          { fprintf(stderr,"ONE2LA: Cannot open %s for reading\n",argv[1]);
            exit (1);
          }
      }
    fclose(output);
    free(root);
    free(pwd);

    schema = oneSchemaCreateFromText(One_Schema);

    file1 = oneFileOpenRead(path,schema,"dal",1);
  }

  if (file1->info['A']->given.count == 0)
    { fprintf(stderr,"ONE2LA: .dal file does not contatin coordinate information\n");
      exit (1);
    }
  if (file1->info['T']->given.count == 0)
    { fprintf(stderr,"ONE2LA: .dal file does not contatin trace information\n");
      exit (1);
    }

  t = oneReadLine(file1);
  if (t == 0 || t != 'X')
    { fprintf(stderr,"ONE2LA: .dal data segment does not begine with an 'X'-line\n");
      exit (1);
    }

  tspace = oneInt(file1,0);
  if (tspace <= TRACE_XOVR && tspace != 0)
    { small  = 1;
      tbytes = 1;
    }
  else
    { small  = 0;
      tbytes = 2;
    }

  Pmax = file1->info['P']->given.max;
  Tmax = file1->info['T']->given.max;

  trace  = Malloc(sizeof(uint16)*Tmax*Pmax,"Allocating trace buffer");
  ovls   = Malloc(sizeof(Overlap)*Pmax,"Allocating overlap vector");
  list   = Malloc(sizeof(int64)*(Pmax+Tmax),"Allocating integer list");
  string = Malloc(Pmax+1,"Allocating max string");

  novls = file1->info['P']->given.total;
  fwrite(&novls,sizeof(int64),1,stdout);
  fwrite(&tspace,sizeof(int),1,stdout);

  while ((t = oneReadLine(file1)) != 0)
    { if (t != 'P')
        { fprintf(stderr,"ONE2LA: Pile data does not begin with a P-line\n");
          exit(1);
        } 
      psize = oneLen(file1);
      list  = oneIntList(file1);
      aread = oneInt(file1,0)-1;
      for (i = 0; i < psize; i++)
        { ovls[i].aread = aread;
          ovls[i].bread = list[i]-1;
        }

      ttop = trace;
      otop = ovls + psize;
      has['O'] = has['C'] = has['A'] = has['B'] = has['L'] = has['D'] = has['T'] = has['Q'] = 0;
      for (o = ovls; o < otop; o++)
        { o->flags = 0;
          o->path.tlen = -1;
        }


      for (j = 0; j < 6+2*psize; j++)
        { t = oneReadLine(file1);

          if (t == 0)
            { fprintf(stderr,"ONE2LA: Pile object not followed by sufficient auxilliary lines\n");
              exit (1);
            }
          if (has[t] > 0 && t != 'T' && t != 'Q')
            { fprintf(stderr,"ONE2LA: Pile has more than one '%c' line\n",t);
              exit (1);
            }
          has[t] += 1;
          if (t == 'A' || t == 'B')
            { if (oneLen(file1) != 2*psize)
                { fprintf(stderr,"ONE2LA: %c-line has incorrect list length\n",t);
                  exit (1);
                }
            }
          else if (t != 'T' && t != 'Q')
            { if (oneLen(file1) != psize)
                { fprintf(stderr,"ONE2LA: %c-line has incorrect list length\n",t);
                  exit (1);
                }
            }
          else
            { if (has[t] > psize)
                { fprintf(stderr,"ONE2LA: Too many %c-lines for pile\n",t);
                  exit (1);
                }
            }

          switch (t)
          { case 'O':
              string = oneString(file1);
              i = 0;
              for (o = ovls; o < otop; o++)
                if (string[i++] == 'c')
                  o->flags |= COMP_FLAG;
              break;
            case 'C':
              string = oneString(file1);
              i = 0;
              for (o = ovls; o < otop; o++)
                if (string[i] == '-')
                  o->flags |= NEXT_FLAG;
                else if (string[i] == '>')
                  o->flags |= BEST_FLAG;
                else if (string[i] == '+')
                  o->flags |= START_FLAG;
              break;
            case 'A':
              list = oneIntList(file1);
              i = 0;
              for (o = ovls; o < otop; o++)
                { o->path.abpos = list[i++];
                  o->path.aepos = list[i++];
                }
              break;
            case 'B':
              list = oneIntList(file1);
              i = 0;
              for (o = ovls; o < otop; o++)
                { o->path.bbpos = list[i++];
                  o->path.bepos = list[i++];
                }
              break;
            case 'L':
              break;
            case 'D':
              list = oneIntList(file1);
              i = 0;
              for (o = ovls; o < otop; o++)
                o->path.diffs = list[i++];
              break;
            case 'T':
            case 'Q':
              list = oneIntList(file1);
              o = ovls + (has[t]-1);
              if (o->path.tlen >= 0)
                { if (o->path.tlen != 2*oneLen(file1))
                    { fprintf(stderr,"LA2ONE: T and Q line lengths do not correspond\n");
                      exit (1);
                    }
                }
              else
                { o->path.tlen  = 2*oneLen(file1);
                  o->path.trace = ttop;
                  ttop += o->path.tlen*tbytes;
                }
              if (t == 'Q')
                k = 0;
              else
                k = 1;
              if (tbytes == 1)
                { uint8 *t8 = (uint8 *) o->path.trace;
                  for (i = 0; k < o->path.tlen; k += 2)
                    t8[k] = list[i++];
                }
              else
                { uint16 *t16 = (uint16 *) o->path.trace;
                  for (i = 0; k < o->path.tlen; k += 2)
                    t16[k] = list[i++];
                }
              break;
            default:
              fprintf(stderr,"LA2ONE: Unrecognized line type '%c'\n",t);
              exit (1);
          }
        }

      if (has['T'] != psize || has['Q'] != psize)
        { fprintf(stderr,"ONE2LA: # of pile traces != pile size\n");
          exit (1);
        }

      for (o = ovls; o < otop; o++)
        Write_Overlap(stdout,o,tbytes);
    }

  oneSchemaDestroy(schema);

  exit (0);
}
