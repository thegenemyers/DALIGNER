/*********************************************************************************************\
 *
 *  Find all local alignment between long, noisy DNA reads:
 *    Compare sequences in 'subject' database against those in the list of 'target' databases
 *    searching for local alignments of 1000bp or more (defined constant MIN_OVERLAP in
 *    filter.c).  Subject is compared in both orientations againt each target.  An output
 *    stream of 'Overlap' records (see align.h) is written in binary to the standard output,
 *    each encoding a given found local alignment between two of the sequences.  The -v
 *    option turns on a verbose reporting mode that gives statistics on each major stage.
 *
 *    There cannot be more than 65,535 reads in a given db, and each read must be less than
 *    66,535 characters long.
 *
 *    The filter operates by looking for a pair of diagonal bands of width 2^'s' that contain
 *    a collection of exact matching 'k'-mers between the two sequences, such that the total
 *    number of bases covered by 'k'-mer hits is 'h'.  k cannot be larger than 15 in the
 *    current implementation.
 *
 *    Some k-mers are significantly over-represented (e.g. homopolymer runs).  These are
 *    suppressed as seed hits, with the parameter 'm' -- any k-mer that occurs more than
 *    'm' times in either the subject or target is not counted as a seed hit.  If the -m
 *    option is absent then no k-mer is suppressed.
 *
 *    For each subject, target pair, say XXX and YYY, the program outputs a file containing
 *    overlaps of the form XXX.YYY.[C|N]#.las where C implies that the reads in XXX were
 *    complemented and N implies they were not (both comparisons are performed), and # is
 *    the thread that detected and wrote out the collection of overlaps.  For example, if
 *    NTHREAD in the program is 4, then 8 files are output for each subject, target pair.
 *
 *  Author:  Gene Myers
 *  Date  :  June 1, 2014
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"
#include "filter.h"


static char *Usage[] =
  { "[-vbd] [-k<int(14)>] [-w<int(6)>] [-h<int(35)>] [-t<int>] [-H<int>]",
    "       [-e<double(.70)] [-m<double(.55)>] [-l<int(1000)>] [-s<int(100)>]",
    "       <subject:file> <target:file> ...",
  };

int VERBOSE;   //   Globally visible to filter.c
int BIASED;
int MINOVER;
int HGAP_MIN;


static HITS_DB *read_DB(char *name, int dust)
{ static HITS_DB  block;
  HITS_TRACK     *dtrack;

  if (Open_DB(name,&block))
    exit (1);

  if (dust)
    dtrack = Load_Track(&block,"dust");
  else
    dtrack = NULL;

  Trim_DB(&block);

  if (block.totlen > 0x7fffffffll)
    { fprintf(stderr,"File (%s) is too large\n",name);
      exit (1);
    }
  if (block.nreads > 0xffff)
    { fprintf(stderr,"There are more than %d reads in file (%s)\n",0xffff,name);
      exit (1);
    }
  if (block.maxlen > 0xffff)
    { fprintf(stderr,"Reads are over %d bases long in file (%s)\n",0xffff,name);
      exit (1);
    }

  Read_All_Sequences(&block,0);

  if (dtrack != NULL)
    { int *anno = (int *) (dtrack->anno); 
      int  i;

      for (i = 0; i <= block.nreads; i++)
        anno[i] /= sizeof(int);
    }

  return (&block);
}

static void complement(char *s, int len)
{ char *t;
  int   c;

  t = s + (len-1);
  while (s < t)
    { c = *s;
      *s = 3-*t;
      *t = 3-c;
      s += 1;
      t -= 1;
    }
  if (s == t)
    *s = 3-*s;
}

static HITS_DB *complement_DB(HITS_DB *block)
{ static HITS_DB cblock;
  int            i, nreads;
  HITS_READ     *reads;
  char          *seq;
  float          x;
  
  nreads = block->nreads;
  reads  = block->reads;
  seq    = (char *) Malloc(block->reads[nreads].boff+1,"Allocating dazzler sequence block");
  if (seq == NULL)
    exit (1);
  *seq++ = 4;
  memcpy(seq,block->bases,block->reads[nreads].boff);

  for (i = 0; i < nreads; i++)
    complement(seq+reads[i].boff,reads[i].end-reads[i].beg);

  cblock = *block;
  cblock.bases = (void *) seq;

  x = cblock.freq[0];
  cblock.freq[0] = cblock.freq[3];
  cblock.freq[3] = x;

  x = cblock.freq[1];
  cblock.freq[1] = cblock.freq[2];
  cblock.freq[2] = x;

  { HITS_TRACK *t, *dust;
    int        *data, *tano, *tata;
    int         j, p, rlen;

    for (t = block->tracks; t != NULL; t++)
      if (strcmp(t->name,"dust") == 0)
        break;
    if (t != NULL)
      { tano = (int *) t->anno;
        tata = (int *) t->data;

        data = (int *) Malloc(sizeof(int)*tano[nreads],"Allocating dazzler .dust index");
        dust = (HITS_TRACK *) Malloc(sizeof(HITS_TRACK),"Allocating dazzler .dust track");
        if (data == NULL || dust == NULL)
          exit (1);

        dust->next = NULL;
        dust->name = t->name;
        dust->size = 4;
        dust->anno = (void *) tano;
        dust->data = (void *) data;
        cblock.tracks = dust;

        p = 0;
        for (i = 0; i < nreads; i++)
          { rlen = (reads[i].end-reads[i].beg)-1;
            for (j = tano[i+1]-1; j >= tano[i]; j--)
              data[p++] = rlen - tata[j];
          }
      }
  }

  return (&cblock);
}

int main(int argc, char *argv[])
{ HITS_DB     ablock,  bblock, cblock;
  char       *afile,  *bfile;
  char       *aroot,  *broot;
  Align_Spec *asettings;

  int    DUSTED;
  int    KMER_LEN;
  int    BIN_SHIFT;
  int    MAX_REPS;
  int    HIT_MIN;
  double AVE_ERROR;
  double MIN_ERROR;
  int    SPACING;

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("daligner")

    KMER_LEN  = 14;
    HIT_MIN   = 35;
    BIN_SHIFT = 6;
    MAX_REPS  = 0;
    HGAP_MIN  = 0;
    AVE_ERROR = .70;
    MIN_ERROR = .55;
    SPACING   = 100;
    MINOVER   = 1000;    //   Globally visible to filter.c

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vbd")
            break;
          case 'k':
            ARG_POSITIVE(KMER_LEN,"K-mer length")
            break;
          case 'w':
            ARG_POSITIVE(BIN_SHIFT,"Log of bin width")
            break;
          case 'h':
            ARG_POSITIVE(HIT_MIN,"Hit threshold (in bp.s)")
            break;
          case 't':
            ARG_POSITIVE(MAX_REPS,"Tuple supression frequency")
            break;
          case 'H':
            ARG_POSITIVE(HGAP_MIN,"HGAP threshold (in bp.s)")
            break;
          case 'e':
            ARG_REAL(AVE_ERROR)
            if (AVE_ERROR < .7 || AVE_ERROR >= 1.)
              { fprintf(stderr,"%s: Average correlation must be in [.7,1.) (%g)\n",
                               Prog_Name,AVE_ERROR);
                exit (1);
              }
            break;
          case 'm':
            ARG_REAL(MIN_ERROR)
            if (MIN_ERROR < .55 || MIN_ERROR >= 1.)
              { fprintf(stderr,"%s: Minimum correlation must be in [.55,1.) (%g)\n",
                               Prog_Name,MIN_ERROR);
                exit (1);
              }
            break;
          case 's':
            ARG_POSITIVE(SPACING,"Trace spacing")
            break;
          case 'l':
            ARG_POSITIVE(MINOVER,"Minimum alignment length")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];   //  Globally declared in filter.h
    BIASED  = flags['b'];   //  Globally declared in filter.h
    DUSTED  = flags['d'];

    if (argc <= 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        exit (1);
      }
  }

  MINOVER *= 2;
  if (Set_Filter_Params(KMER_LEN,BIN_SHIFT,MAX_REPS,HIT_MIN))
    { fprintf(stderr,"Illegal combination of filter parameters\n");
      exit (1);
    }

  /* Read in the reads in A */

  afile  = argv[1];
  aroot  = Root(afile,".db");
  ablock = *read_DB(afile,DUSTED);
  cblock = *complement_DB(&ablock);

  if (ablock.cutoff >= HGAP_MIN)
    HGAP_MIN = ablock.cutoff;

  asettings = New_Align_Spec( AVE_ERROR, MIN_ERROR, SPACING, ablock.freq);

  /* Build indices for A and A complement */

  Build_Table(&ablock);
  Build_Table(&cblock);

  /* Compare against reads in B in both orientations */

  { int i;

    for (i = 2; i < argc; i++)
      { bfile = argv[i];
        broot = Root(bfile,".db");
        if (strcmp(afile,bfile) == 0)
          { Match_Filter(aroot,&ablock,broot,&ablock,1,0,asettings);
            Match_Filter(aroot,&cblock,broot,&ablock,1,1,asettings);
          }
        else
          { bblock = *read_DB(bfile,DUSTED);
            Match_Filter(aroot,&ablock,broot,&bblock,0,0,asettings);
            Match_Filter(aroot,&cblock,broot,&bblock,0,1,asettings);
            Close_DB(&bblock);
          }
      }
  }

  exit (0);
}
