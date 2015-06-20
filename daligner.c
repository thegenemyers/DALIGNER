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
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <sys/param.h>
#if defined(BSD)
#include <sys/sysctl.h>
#endif

#include "DB.h"
#include "filter.h"

static char *Usage[] =
  { "[-vbAI] [-k<int(14)>] [-w<int(6)>] [-h<int(35)>] [-t<int>] [-M<int>]",
    "        [-e<double(.70)] [-l<int(1000)>] [-s<int(100)>] [-H<int>]",
    "        [-m<track>]+ <subject:db|dam> <target:db|dam> ...",
  };

int     VERBOSE;   //   Globally visible to filter.c
int     BIASED;
int     MINOVER;
int     HGAP_MIN;
int     SYMMETRIC;
int     IDENTITY;
uint64  MEM_LIMIT;
uint64  MEM_PHYSICAL;

/*  Adapted from code by David Robert Nadeau (http://NadeauSoftware.com) licensed under
 *     "Creative Commons Attribution 3.0 Unported License"
 *          (http://creativecommons.org/licenses/by/3.0/deed.en_US)
 *
 *   I removed Windows options, reformated, and return int64 instead of size_t
 */

static int64 getMemorySize( )
{
#if defined(CTL_HW) && (defined(HW_MEMSIZE) || defined(HW_PHYSMEM64))

  // OSX, NetBSD, OpenBSD

  int     mib[2];
  size_t  size = 0;
  size_t  len = sizeof( size );

  mib[0] = CTL_HW;
#if defined(HW_MEMSIZE)
  mib[1] = HW_MEMSIZE;            // OSX
#elif defined(HW_PHYSMEM64)
  mib[1] = HW_PHYSMEM64;          // NetBSD, OpenBSD
#endif
  if (sysctl(mib,2,&size,&len,NULL,0) == 0)
    return ((size_t) size);
  return (0);

#elif defined(_SC_AIX_REALMEM)

  // AIX

  return ((size_t) sysconf( _SC_AIX_REALMEM ) * ((size_t) 1024L));

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)

  // FreeBSD, Linux, OpenBSD, & Solaris

  size_t  size = 0;

  size = (size_t) sysconf(_SC_PHYS_PAGES);
  return (size * ((size_t) sysconf(_SC_PAGESIZE)));

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGE_SIZE)

  // ? Legacy ?

  size_t  size = 0;

  size = (size_t) sysconf(_SC_PHYS_PAGES);
  return (size * ((size_t) sysconf(_SC_PAGE_SIZE)));

#elif defined(CTL_HW) && (defined(HW_PHYSMEM) || defined(HW_REALMEM))

  // DragonFly BSD, FreeBSD, NetBSD, OpenBSD, and OSX

  int          mib[2];
  unsigned int size = 0;
  size_t       len  = sizeof( size );

  mib[0] = CTL_HW;
#if defined(HW_REALMEM)
  mib[1] = HW_REALMEM;		// FreeBSD
#elif defined(HW_PYSMEM)
  mib[1] = HW_PHYSMEM;		// Others
#endif
  if (sysctl(mib,2,&size,&len,NULL,0) == 0)
    return (size_t)size;
  return (0);

#else

  return (0);

#endif
}

typedef struct
  { int *ano;
    int *end;
    int  idx;
    int  out;
  } Event;

static void reheap(int s, Event **heap, int hsize)
{ int      c, l, r;
  Event   *hs, *hr, *hl;

  c  = s;
  hs = heap[s];
  while ((l = 2*c) <= hsize)
    { r  = l+1;
      hl = heap[l];
      hr = heap[r];
      if (hr->idx > hl->idx)
        { if (hs->idx > hl->idx)
            { heap[c] = hl;
              c = l;
            }
          else
            break;
        }
      else
        { if (hs->idx > hr->idx)
            { heap[c] = hr;
              c = r;
            }
          else
            break;
        }
    }
  if (c != s)
    heap[c] = hs;
}

int64 Merge_Size(HITS_DB *block, int mtop)
{ Event       ev[mtop+1];
  Event      *heap[mtop+2];
  int         r, mhalf;
  int64       nsize;

  { HITS_TRACK *track;
    int         i;

    track = block->tracks;
    for (i = 0; i < mtop; i++)
      { ev[i].ano = ((int *) (track->data)) + ((int64 *) (track->anno))[0];
        ev[i].out = 1;
        heap[i+1] = ev+i;
        track = track->next;
      }
    ev[mtop].idx = INT32_MAX;
    heap[mtop+1] = ev+mtop;
  }

  mhalf = mtop/2;

  nsize = 0;
  for (r = 0; r < block->nreads; r++)
    { int         i, level, hsize;
      HITS_TRACK *track;

      track = block->tracks;
      for (i = 0; i < mtop; i++)
        { ev[i].end = ((int *) (track->data)) + ((int64 *) (track->anno))[r+1];
          if (ev[i].ano < ev[i].end)
            ev[i].idx = *(ev[i].ano);
          else
            ev[i].idx = INT32_MAX;
          track = track->next;
        }
      hsize = mtop;

      for (i = mhalf; i > 1; i--)
        reheap(i,heap,hsize);

      level = 0;
      while (1)
        { Event *p;

          reheap(1,heap,hsize);

          p = heap[1];
          if (p->idx == INT32_MAX) break;

          p->out = 1-p->out;
          if (p->out)
            { level -= 1;
              if (level == 0)
                nsize += 1;
            }
          else
            { if (level == 0)
                nsize += 1;
              level += 1;
            }
          p->ano += 1;
          if (p->ano >= p->end)
            p->idx = INT32_MAX;
          else
            p->idx = *(p->ano);
        }
    }

  return (nsize);
}

HITS_TRACK *Merge_Tracks(HITS_DB *block, int mtop, int64 nsize)
{ HITS_TRACK *ntrack;
  Event       ev[mtop+1];
  Event      *heap[mtop+2];
  int         r, mhalf;
  int64      *anno;
  int        *data;

  ntrack = (HITS_TRACK *) Malloc(sizeof(HITS_TRACK),"Allocating merged track");
  if (ntrack == NULL)
    exit (1);
  ntrack->name = Strdup("merge","Allocating merged track");
  ntrack->anno = anno = (int64 *) Malloc(sizeof(int64)*(block->nreads+1),"Allocating merged track");
  ntrack->data = data = (int *) Malloc(sizeof(int)*nsize,"Allocating merged track");
  ntrack->size = sizeof(int);
  ntrack->next = NULL;
  if (anno == NULL || data == NULL || ntrack->name == NULL)
    exit (1);

  { HITS_TRACK *track;
    int         i;

    track = block->tracks;
    for (i = 0; i < mtop; i++)
      { ev[i].ano = ((int *) (track->data)) + ((int64 *) (track->anno))[0];
        ev[i].out = 1;
        heap[i+1] = ev+i;
        track = track->next;
      }
    ev[mtop].idx = INT32_MAX;
    heap[mtop+1] = ev+mtop;
  }

  mhalf = mtop/2;

  nsize = 0;
  for (r = 0; r < block->nreads; r++)
    { int         i, level, hsize;
      HITS_TRACK *track;

      anno[r] = nsize;

      track = block->tracks;
      for (i = 0; i < mtop; i++)
        { ev[i].end = ((int *) (track->data)) + ((int64 *) (track->anno))[r+1];
          if (ev[i].ano < ev[i].end)
            ev[i].idx = *(ev[i].ano);
          else
            ev[i].idx = INT32_MAX;
          track = track->next;
        }
      hsize = mtop;

      for (i = mhalf; i > 1; i--)
        reheap(i,heap,hsize);

      level = 0;
      while (1)
        { Event *p;

          reheap(1,heap,hsize);

          p = heap[1];
          if (p->idx == INT32_MAX) break;

          p->out = 1-p->out;
          if (p->out)
            { level -= 1;
              if (level == 0)
                data[nsize++] = p->idx;
            }
          else
            { if (level == 0)
                data[nsize++] = p->idx;
              level += 1;
            }
          p->ano += 1;
          if (p->ano >= p->end)
            p->idx = INT32_MAX;
          else
            p->idx = *(p->ano);
        }
    }
  anno[r] = nsize;

  return (ntrack);
}

static HITS_DB *read_DB(char *name, char **mask, int *mstat, int mtop, int kmer, int *isdam)
{ static HITS_DB  block;
  int             i, status, stop;

  status = Open_DB(name,&block);
  if (status < 0)
    exit (1);
  *isdam = status;

  for (i = 0; i < mtop; i++)
    { status = Check_Track(&block,mask[i]);
      if (status > mstat[i])
        mstat[i] = status;
      if (status == 0)
        Load_Track(&block,mask[i]);
    }

  Trim_DB(&block);

  stop = 0;
  for (i = 0; i < mtop; i++)
    { HITS_TRACK *track;
      int64      *anno;
      int         j;

      status = Check_Track(&block,mask[i]);
      if (status < 0)
        continue;
      stop += 1;
      track = Load_Track(&block,mask[i]);

      anno = (int64 *) (track->anno); 
      for (j = 0; j <= block.nreads; j++)
        anno[j] /= sizeof(int);
    }

  if (stop > 1)
    { int64       nsize;
      HITS_TRACK *track;

      nsize = Merge_Size(&block,stop);
      track = Merge_Tracks(&block,stop,nsize);

      while (block.tracks != NULL)
        Close_Track(&block,block.tracks->name);

      block.tracks = track;
    }

  if (block.cutoff < kmer)
    { for (i = 0; i < block.nreads; i++)
        if (block.reads[i].rlen < kmer)
          { fprintf(stderr,"%s: Block %s contains reads < %dbp long !  Run DBsplit.\n",
                           Prog_Name,name,kmer);
            exit (1);
          }
    }

  Read_All_Sequences(&block,0);

  return (&block);
}

static void complement(char *s, int len)
{ char *t;
  int   c;

  t = s + (len-1);
  while (s < t)
    { c = *s;
      *s = (char) (3-*t);
      *t = (char) (3-c);
      s += 1;
      t -= 1;
    }
  if (s == t)
    *s = (char) (3-*s);
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
    complement(seq+reads[i].boff,reads[i].rlen);

  cblock = *block;
  cblock.bases = (void *) seq;

  x = cblock.freq[0];
  cblock.freq[0] = cblock.freq[3];
  cblock.freq[3] = x;

  x = cblock.freq[1];
  cblock.freq[1] = cblock.freq[2];
  cblock.freq[2] = x;

  { HITS_TRACK *t, *dust;
    int        *data, *tata;
    int         p, rlen;
    int64       j, *tano;

    for (t = block->tracks; t != NULL; t = t->next)
      { tano = (int64 *) t->anno;
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
          { rlen = reads[i].rlen;
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
  int         isdam;
  int         MMAX, MTOP, *MSTAT;
  char      **MASK;

  int    KMER_LEN;
  int    BIN_SHIFT;
  int    MAX_REPS;
  int    HIT_MIN;
  double AVE_ERROR;
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
    SPACING   = 100;
    MINOVER   = 1000;    //   Globally visible to filter.c

    MEM_PHYSICAL = getMemorySize();
    MEM_LIMIT    = MEM_PHYSICAL;
    if (MEM_PHYSICAL == 0)
      { fprintf(stderr,"\nWarning: Could not get physical memory size\n");
        fflush(stderr);
      }

    MTOP  = 0;
    MMAX  = 10;
    MASK  = (char **) Malloc(MMAX*sizeof(char *),"Allocating mask track array");
    MSTAT = (int *) Malloc(MMAX*sizeof(int),"Allocating mask status array");
    if (MASK == NULL || MSTAT == NULL)
      exit (1);

    j    = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vbAI")
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
          case 'l':
            ARG_POSITIVE(MINOVER,"Minimum alignment length")
            break;
          case 's':
            ARG_POSITIVE(SPACING,"Trace spacing")
            break;
          case 'M':
            { int limit;

              ARG_NON_NEGATIVE(limit,"Memory allocation (in Gb)")
              MEM_LIMIT = limit * 0x40000000ll;
              break;
            }
          case 'm':
            if (MTOP >= MMAX)
              { MMAX  = 1.2*MTOP + 10;
                MASK  = (char **) Realloc(MASK,MMAX*sizeof(char *),"Reallocating mask track array");
                MSTAT = (int *) Realloc(MSTAT,MMAX*sizeof(int),"Reallocating mask status array");
                if (MASK == NULL || MSTAT == NULL)
                  exit (1);
              }
            MSTAT[MTOP]  = -2;
            MASK[MTOP++] = argv[i]+2;
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE   = flags['v'];   //  Globally declared in filter.h
    BIASED    = flags['b'];   //  Globally declared in filter.h
    SYMMETRIC = 1-flags['A'];
    IDENTITY  = flags['I'];

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
  ablock = *read_DB(afile,MASK,MSTAT,MTOP,KMER_LEN,&isdam);
  cblock = *complement_DB(&ablock);
  if (isdam)
    aroot = Root(afile,".dam");
  else
    aroot = Root(afile,".db");

  if (ablock.cutoff >= HGAP_MIN)
    HGAP_MIN = ablock.cutoff;

  asettings = New_Align_Spec( AVE_ERROR, SPACING, ablock.freq);

  /* Compare against reads in B in both orientations */

  { int i, j;

    broot = NULL;
    for (i = 2; i < argc; i++)
      { bfile = argv[i];
        if (strcmp(afile,bfile) != 0)
          { bblock = *read_DB(bfile,MASK,MSTAT,MTOP,KMER_LEN,&isdam);
            if (isdam)
              broot = Root(bfile,".dam");
            else
              broot = Root(bfile,".db");
          }

        if (i == 2)
          { for (j = 0; j < MTOP; j++)
              { if (MSTAT[j] == -2)
                  printf("%s: Warning: -m%s option given but no track found.\n",Prog_Name,MASK[i]);
                else if (MSTAT[j] == -1)
                  printf("%s: Warning: %s track not sync'd with relevant db.\n",Prog_Name,MASK[i]);
              }

            if (VERBOSE)
              printf("\nBuilding index for %s\n",aroot);
            Build_Table(&ablock);

            if (VERBOSE)
              printf("\nBuilding index for c(%s)\n",aroot);
            Build_Table(&cblock);
          }
      
        if (strcmp(afile,bfile) == 0)
          { Match_Filter(aroot,&ablock,aroot,&ablock,1,0,asettings);
            Match_Filter(aroot,&cblock,aroot,&ablock,1,1,asettings);
          }
        else
          { Match_Filter(aroot,&ablock,broot,&bblock,0,0,asettings);
            Match_Filter(aroot,&cblock,broot,&bblock,0,1,asettings);
            Close_DB(&bblock);
            free(broot);
          }
      }
  }

  exit (0);
}
