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
 *  Compressed data base module.  Auxiliary routines to open and manipulate a data base for
 *    which the sequence and read information are separated into two separate files, and the
 *    sequence is compressed into 2-bits for each base.  Support for tracks of additional
 *    information, and trimming according to the current partition.  Eventually will also
 *    support compressed quality information.
 *
 *  Author :  Gene Myers
 *  Date   :  July 2013
 *  Revised:  April 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <dirent.h>

#include "DB.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif


/*******************************************************************************************
 *
 *  GENERAL UTILITIES
 *
 ********************************************************************************************/

char *Prog_Name;

void *Malloc(int64 size, char *mesg)
{ void *p;

  if ((p = malloc(size)) == NULL)
    { if (mesg == NULL)
        fprintf(stderr,"%s: Out of memory\n",Prog_Name);
      else
        fprintf(stderr,"%s: Out of memory (%s)\n",Prog_Name,mesg);
    }
  return (p);
}

void *Realloc(void *p, int64 size, char *mesg)
{ if ((p = realloc(p,size)) == NULL)
    { if (mesg == NULL)
        fprintf(stderr,"%s: Out of memory\n",Prog_Name);
      else
        fprintf(stderr,"%s: Out of memory (%s)\n",Prog_Name,mesg);
    }
  return (p);
}

char *Strdup(char *name, char *mesg)
{ char *s;

  if (name == NULL)
    return (NULL);
  if ((s = strdup(name)) == NULL)
    { if (mesg == NULL)
        fprintf(stderr,"%s: Out of memory\n",Prog_Name);
      else
        fprintf(stderr,"%s: Out of memory (%s)\n",Prog_Name,mesg);
    }
  return (s);
}

FILE *Fopen(char *name, char *mode)
{ FILE *f;

  if (name == NULL || mode == NULL)
    return (NULL);
  if ((f = fopen(name,mode)) == NULL)
    fprintf(stderr,"%s: Cannot open %s for '%s'\n",Prog_Name,name,mode);
  return (f);
}

char *PathTo(char *name)
{ char *path, *find;

  if (name == NULL)
    return (NULL);
  if ((find = rindex(name,'/')) != NULL)
    { *find = '\0';
      path = Strdup(name,"Extracting path from");
      *find = '/';
    }
  else
    path = Strdup(".","Allocating default path");
  return (path);
}

char *Root(char *name, char *suffix)
{ char *path, *find, *dot;
  int   epos;

  if (name == NULL)
    return (NULL);
  find = rindex(name,'/');
  if (find == NULL)
    find = name;
  else
    find += 1;
  if (suffix == NULL)
    { dot = strchr(find,'.');
      if (dot != NULL)
        *dot = '\0';
      path = Strdup(find,"Extracting root from");
      if (dot != NULL)
        *dot = '.';
    }
  else
    { epos = strlen(find) - strlen(suffix);
      if (epos > 0 && strcasecmp(find+epos,suffix) == 0)
        { find[epos] = '\0';
          path = Strdup(find,"Extracting root from");
          find[epos] = suffix[0];
        }
      else
        path = Strdup(find,"Allocating root");
    }
  return (path);
}

char *Catenate(char *path, char *sep, char *root, char *suffix)
{ static char *cat = NULL;
  static int   max = -1;
  int len;

  if (path == NULL || root == NULL || sep == NULL || suffix == NULL)
    return (NULL);
  len = strlen(path) + strlen(sep) + strlen(root) + strlen(suffix);
  if (len > max)
    { max = 1.2*len + 100;
      if ((cat = (char *) realloc(cat,max+1)) == NULL)
        { fprintf(stderr,"%s: Out of memory (Making path name for %s)\n",Prog_Name,root);
          return (NULL);
        }
    }
  sprintf(cat,"%s%s%s%s",path,sep,root,suffix);
  return (cat);
}

char *Numbered_Suffix(char *left, int num, char *right)
{ static char *suffix = NULL;
  static int   max = -1;
  int len;

  if (left == NULL || right == NULL)
    return (NULL);
  len = strlen(left) + strlen(right) + 40;
  if (len > max)
    { max = 1.2*len + 100;
      if ((suffix = (char *) realloc(suffix,max+1)) == NULL)
        { fprintf(stderr,"%s: Out of memory (Making number suffix for %d)\n",Prog_Name,num);
          return (NULL);
        }
    }
  sprintf(suffix,"%s%d%s",left,num,right);
  return (suffix);
}


#define  COMMA  ','

//  Print big integers with commas/periods for better readability

void Print_Number(int64 num, int width, FILE *out)
{ if (width == 0)
    { if (num < 1000ll)
        fprintf(out,"%lld",num);
      else if (num < 1000000ll)
        fprintf(out,"%lld%c%03lld",num/1000ll,COMMA,num%1000ll);
      else if (num < 1000000000ll)
        fprintf(out,"%lld%c%03lld%c%03lld",num/1000000ll,
                                           COMMA,(num%1000000ll)/1000ll,COMMA,num%1000ll);
      else
        fprintf(out,"%lld%c%03lld%c%03lld%c%03lld",num/1000000000ll,
                                                  COMMA,(num%1000000000ll)/1000000ll,
                                                  COMMA,(num%1000000ll)/1000ll,COMMA,num%1000ll);
    }
  else
    { if (num < 1000ll)
        fprintf(out,"%*lld",width,num);
      else if (num < 1000000ll)
        fprintf(out,"%*lld%c%03lld",width-4,num/1000ll,COMMA,num%1000ll);
      else if (num < 1000000000ll)
        fprintf(out,"%*lld%c%03lld%c%03lld",width-8,num/1000000ll,COMMA,(num%1000000ll)/1000ll,
                                            COMMA,num%1000ll);
      else
        fprintf(out,"%*lld%c%03lld%c%03lld%c%03lld",width-12,num/1000000000ll,COMMA,
                                                        (num%1000000000ll)/1000000ll,COMMA,
                                                        (num%1000000ll)/1000ll,COMMA,num%1000ll);
    }
}


/*******************************************************************************************
 *
 *  READ COMPRESSION/DECOMPRESSION UTILITIES
 *
 ********************************************************************************************/

//  Compress read into 2-bits per base (from [0-3] per byte representation

void Compress_Read(int len, char *s)
{ int   i, c, d;
  char *s0, *s1, *s2, *s3;

  s0 = s;
  s1 = s0+1;
  s2 = s1+1;
  s3 = s2+1;

  c = s1[len];
  d = s2[len];
  s0[len] = s1[len] = s2[len] = 0;

  for (i = 0; i < len; i += 4)
    *s++ = (s0[i] << 6) | (s1[i] << 4) | (s2[i] << 2) | s3[i];

  s1[len] = c;
  s2[len] = d;
}

//  Uncompress read form 2-bits per base into [0-3] per byte representation

void Uncompress_Read(int len, char *s)
{ int   i, tlen, byte;
  char *s0, *s1, *s2, *s3;
  char *t;

  s0 = s;
  s1 = s0+1;
  s2 = s1+1;
  s3 = s2+1;

  tlen = (len-1)/4;

  t = s+tlen;
  for (i = tlen*4; i >= 0; i -= 4)
    { byte = *t--;
      s0[i] = ((byte >> 6) & 0x3);
      s1[i] = ((byte >> 4) & 0x3);
      s2[i] = ((byte >> 2) & 0x3);
      s3[i] = (byte & 0x3);
    }
  s[len] = 4;
}

//  Convert read in [0-3] representation to ascii representation (end with '\n')

void Lower_Read(char *s)
{ static char letter[4] = { 'a', 'c', 'g', 't' };

  for ( ; *s != 4; s++)
    *s = letter[(int) *s];
  *s = '\0';
}

void Upper_Read(char *s)
{ static char letter[4] = { 'A', 'C', 'G', 'T' };

  for ( ; *s != 4; s++)
    *s = letter[(int) *s];
  *s = '\0';
}

//  Convert read in ascii representation to [0-3] representation (end with 4)

void Number_Read(char *s)
{ static char number[128] =
    { 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
    };

  for ( ; *s != '\0'; s++)
    *s = number[(int) *s];
  *s = 4;
}


/*******************************************************************************************
 *
 *  DB OPEN, TRIM & CLOSE ROUTINES
 *
 ********************************************************************************************/

// Open the given database "root" into the supplied HITS_DB record "db"
//   The index array is allocated and read in, the 'bases' file is open for reading.

int Open_DB(char* path, HITS_DB *db)
{ char *root, *pwd, *bptr, *fptr;
  int   nreads;
  FILE *index, *dbvis;
  int   status;
  int   part, cutoff, all;
  int   ofirst, bfirst, olast;

  status = 0;

  root = Root(path,".db");
  pwd  = PathTo(path);

  bptr = rindex(root,'.');
  if (bptr != NULL && bptr[1] != '\0' && bptr[1] != '-')
    { part = strtol(bptr+1,&fptr,10);
      if (*fptr != '\0' || part == 0)
        part = 0;
      else
        *bptr = '\0';
    }
  else
    part = 0;

  if (db->cutoff < 0 && part > 0)
    { fprintf(stderr,"%s: DB %s has not yet been partitioned, cannot request a block !\n",
                     Prog_Name,root);
      goto exit;
    }
  
  if ((dbvis = Fopen(Catenate(pwd,"/",root,".db"),"r")) == NULL)
    { status = 1;
      goto exit;
    }

  if ((index = Fopen(Catenate(pwd,PATHSEP,root,".idx"),"r")) == NULL)
    { status = 1;
      goto exit1;
    }
  fread(db,sizeof(HITS_DB),1,index);
  nreads = db->oreads;

  { int   p, nblocks, nfiles, blast;
    int64 size;
    char buffer[2*MAX_NAME+100];

    nblocks = 0;
    fscanf(dbvis,DB_NFILE,&nfiles);
    for (p = 0; p < nfiles; p++)
      fgets(buffer,2*MAX_NAME+100,dbvis);
    if (fscanf(dbvis,DB_NBLOCK,&nblocks) != 1 || part > nblocks)
      if (part > 0)
        { status = 1;
          if (nblocks == 0)
            fprintf(stderr,"%s: DB has not been partitioned\n",Prog_Name);
          else
            fprintf(stderr,"%s: DB has only %d blocks\n",Prog_Name,nblocks);
          goto exit2;
        }
      else
        { cutoff = 0;
          all    = 1;
        }
    else
      fscanf(dbvis,DB_PARAMS,&size,&cutoff,&all);

    if (part > 0)
      { for (p = 1; p <= part; p++)
          fscanf(dbvis,DB_BDATA,&ofirst,&bfirst);
        fscanf(dbvis,DB_BDATA,&olast,&blast);
      }
    else
      { ofirst = bfirst = 0;
        olast  = nreads;
      }
  }

  db->trimmed = 0;
  db->tracks  = NULL;
  db->part    = part;
  db->cutoff  = cutoff;
  db->all     = all;
  db->ofirst  = ofirst;
  db->bfirst  = bfirst;

  if (part <= 0)
    { db->reads = (HITS_READ *) Malloc(sizeof(HITS_READ)*(nreads+1),"Allocating Open_DB index");
      fread(db->reads,sizeof(HITS_READ),nreads,index);
    }
  else
    { HITS_READ *reads;
      int        i, r, maxlen;
      int64      totlen;

      nreads = olast-ofirst;
      reads  = (HITS_READ *) Malloc(sizeof(HITS_READ)*(nreads+1),"Allocating Open_DB index");

      fseeko(index,sizeof(HITS_READ)*ofirst,SEEK_CUR);
      fread(reads,sizeof(HITS_READ),nreads,index);

      totlen = 0;
      maxlen = 0;
      for (i = 0; i < nreads; i++)
        { r = reads[i].end - reads[i].beg;
          totlen += r;
          if (r > maxlen)
            maxlen = r;
        }

      db->maxlen = maxlen;
      db->totlen = totlen;
      db->reads  = reads;
    }

  db->nreads = nreads;
  db->path   = Strdup(Catenate(pwd,PATHSEP,root,""),"Allocating Open_DB path");
  db->bases  = NULL;
  db->loaded = 0;

exit2:
  fclose(index);
exit1:
  fclose(dbvis);
exit:
  if (bptr != NULL)
    *bptr = '.';

  return (status);
}


// Trim the DB or part thereof and all loaded tracks according to the cuttof and all settings
//   of the current DB partition.  Reallocate smaller memory blocks for the information kept
//   for the retained reads.

void Trim_DB(HITS_DB *db)
{ int         i, j, r;
  int         allflag, cutoff;
  int64       totlen;
  int         maxlen, nreads;
  HITS_TRACK *record;
  HITS_READ  *reads;

  if (db->trimmed) return;

  if (db->cutoff <= 0 && db->all) return;

  cutoff = db->cutoff;
  if (db->all)
    allflag = 0;
  else
    allflag = DB_BEST;

  reads  = db->reads;
  nreads = db->nreads;

  for (record = db->tracks; record != NULL; record = record->next)
    if (strcmp(record->name,".@qvs") == 0)
      { uint16 *table = ((HITS_QV *) record)->table;

        j = 0;
        for (i = 0; i < db->nreads; i++)
          if ((reads[i].flags & DB_BEST) >= allflag && reads[i].end - reads[i].beg >= cutoff)
            table[j++] = table[i];
      }
    else
      { int   *anno4, size;
        int64 *anno8;
        void  *anno, *data;

        size = record->size;
        data = record->data; 
        if (data == NULL)
          { anno = record->anno;
            j = 0;
            for (i = r = 0; i < db->nreads; i++, r += size)
              if ((reads[i].flags & DB_BEST) >= allflag && reads[i].end - reads[i].beg >= cutoff)
                { memmove(anno+j,anno+r,size);
                  j += size;
                }
            memmove(anno+j,anno+r,size);
          }
        else if (size == 4)
          { int ai;

            anno4 = (int *) (record->anno);
            j = anno4[0] = 0;
            for (i = 0; i < db->nreads; i++)
              if ((reads[i].flags & DB_BEST) >= allflag && reads[i].end - reads[i].beg >= cutoff)
                { ai = anno4[i];
                  anno4[j+1] = anno4[j] + (anno4[i+1]-ai);
                  memmove(data+anno4[j],data+ai,anno4[i+1]-ai);
                  j += 1;
                }
            record->data = Realloc(record->data,anno4[j],NULL);
          }
        else // size == 8
          { int64 ai;

            anno8 = (int64 *) (record->anno);
            j = anno8[0] = 0;
            for (i = 0; i < db->nreads; i++)
              if ((reads[i].flags & DB_BEST) >= allflag && reads[i].end - reads[i].beg >= cutoff)
                { ai = anno8[i];
                  anno8[j+1] = anno8[j] + (anno8[i+1]-ai);
                  memmove(data+anno8[j],data+ai,anno8[i+1]-ai);
                  j += 1;
                }
            record->data = Realloc(record->data,anno8[j],NULL);
          }
        record->anno = Realloc(record->anno,record->size*(j+1),NULL);
      }

  totlen = maxlen = 0;
  for (j = i = 0; i < nreads; i++)
    { r = reads[i].end - reads[i].beg;
      if ((reads[i].flags & DB_BEST) >= allflag && r >= cutoff)
        { totlen += r;
          if (r > maxlen)
            maxlen = r;
          reads[j++] = reads[i];
        }
    }
  
  db->totlen  = totlen;
  db->maxlen  = maxlen;
  db->nreads  = j;
  db->trimmed = 1;

  if (j < nreads)
    reads = Realloc(reads,sizeof(HITS_READ)*(j+1),NULL);
}

// Shut down an open 'db' by freeing all associated space, including tracks and QV structures, 
//   and any open file pointers.  The record pointed at by db however remains (the user
//   supplied it and so should free it).

void Close_DB(HITS_DB *db)
{ HITS_TRACK *t, *p;

  if (db->loaded)
    free(((char *) (db->bases)) - 1);
  else if (db->bases != NULL)
    fclose((FILE *) db->bases);
  free(db->reads);
  free(db->path);

  Close_QVs(db);

  for (t = db->tracks; t != NULL; t = p)
    { p = t->next;
      free(t->anno);
      free(t->data);
      free(t);
    }
}


/*******************************************************************************************
 *
 *  QV LOAD & CLOSE ROUTINES
 *
 ********************************************************************************************/

HITS_DB *Active_DB = NULL;  //  Last db/qv used by "Load_QVentry"
HITS_QV *Active_QV;         //    Becomes invalid after closing

void Load_QVs(HITS_DB *db)
{ FILE        *quiva, *istub, *indx;
  uint16      *table;
  QVcoding    *coding;
  HITS_QV     *qvtrk;

  if (db->tracks != NULL && strcmp(db->tracks->name,".@qvs") == 0)
    return;

  if (db->trimmed)
    { fprintf(stderr,"%s: Cannot load QVs after trimming the DB\n",Prog_Name);
      exit (1);
    }

  if (db->reads[db->nreads-1].coff == 0)
    { fprintf(stderr,"%s: The requested QVs have not been added to the DB!\n",Prog_Name);
      exit (1);
    }

  //  Open .qvs, .idx, and .db files

  quiva = Fopen(Catenate(db->path,"","",".qvs"),"r");
  if (quiva == NULL)
    return;
  if (db->part > 0)
    { indx = Fopen(Catenate(db->path,"","",".idx"),"r");
      if (indx == NULL)
        exit (1);
    }

  { char *x;

    x = rindex(db->path,'/');
    *x = '\0';
    istub = Fopen(Catenate(db->path,"/",x+2,".db"),"r");
    if (istub == NULL)
      exit (1);
    *x = '/';
  }

  { int   first, last, nfiles, ncodes;
    char  prolog[MAX_NAME], fname[MAX_NAME];
    int   i, j;

    table = (uint16 *) Malloc(sizeof(uint16)*db->nreads,"Allocating QV table indices");

    fscanf(istub,DB_NFILE,&nfiles);

    if (db->part > 0)
      { int pfirst, plast;
        int fbeg, fend;
        int n, k;

        //  Determine first how many and which files span the block (fbeg to fend)

        pfirst = db->ofirst;
        plast  = pfirst + db->nreads;

        first = 0;
        for (fbeg = 0; fbeg < nfiles; fbeg++)
          { fscanf(istub,DB_FDATA,&last,fname,prolog);
            if (last > pfirst)
              break;
            first = last;
          }
        for (fend = fbeg+1; fend <= nfiles; fend++)
          { if (last >= plast)
              break;
            fscanf(istub,DB_FDATA,&last,fname,prolog);
            first = last;
          }

        ncodes = fend-fbeg;
        coding = (QVcoding *) Malloc(sizeof(QVcoding)*ncodes,"Allocating coding schemes");

        //  Carefully get the first coding scheme (its offset is most likely in a HITS_RECORD
        //    in .idx that is *not* in memory).  Get all the other coding schemes normally and
        //    assign the tables # for each read in the block in "tables".

        rewind(istub);
        fscanf(istub,DB_NFILE,&nfiles);

        first = 0;
        for (n = 0; n < fbeg; n++)
          { fscanf(istub,DB_FDATA,&last,fname,prolog);
            first = last;
          }

        for (n = fbeg; n < fend; n++)
          { fscanf(istub,DB_FDATA,&last,fname,prolog);

            i = n-fbeg;
            if (first < pfirst)
              { HITS_READ read;

                fseeko(indx,sizeof(HITS_DB) + sizeof(HITS_READ)*first,SEEK_SET);
                fread(&read,sizeof(HITS_READ),1,indx);
                fseeko(quiva,read.coff,SEEK_SET);
                coding[i] = *Read_QVcoding(quiva);
              }
            else
              { fseeko(quiva,db->reads[first-pfirst].coff,SEEK_SET);
                coding[i] = *Read_QVcoding(quiva);
                db->reads[first-pfirst].coff = ftello(quiva);
              }

            j = first-pfirst;
            if (j < 0)
              j = 0;
            k = last-pfirst;
            if (k > db->nreads)
              k = db->nreads;
            while (j < k)
              table[j++] = i;

            first = last;
          }

        fclose(indx);
      }

    else
      { //  Load in coding scheme for each file, adjust .coff of first read in the file, and
        //    record which table each read uses

        ncodes = nfiles;
        coding = (QVcoding *) Malloc(sizeof(QVcoding)*nfiles,"Allocating coding schemes");
  
        first = 0;
        for (i = 0; i < nfiles; i++)
          { fscanf(istub,DB_FDATA,&last,fname,prolog);
  
            fseeko(quiva,db->reads[first].coff,SEEK_SET);
            coding[i] = *Read_QVcoding(quiva);
	    db->reads[first].coff = ftello(quiva);

            for (j = first; j < last; j++)
              table[j] = i;

            first = last;
          }
      }

    //  Allocate and fill in the HITS_QV record and add it to the front of the
    //    track list

    qvtrk = (HITS_QV *) Malloc(sizeof(HITS_QV),"Allocating QV pseudo-track");
    qvtrk->next   = db->tracks;
    db->tracks    = (HITS_TRACK *) qvtrk;
    qvtrk->name   = Strdup(".@qvs","Allocating QV pseudo-track name");
    qvtrk->ncodes = ncodes;
    qvtrk->table  = table;
    qvtrk->coding = coding;
    qvtrk->quiva  = quiva;
  }

  fclose(istub);
  fclose(indx);
}

// Close the QV stream, free the QV pseudo track and all associated memory

void Close_QVs(HITS_DB *db)
{ HITS_TRACK *track;
  HITS_QV    *qvtrk;
  int         i;

  Active_DB = NULL;

  track = db->tracks;
  if (track != NULL && strcmp(track->name,".@qvs") == 0)
    { qvtrk = (HITS_QV *) track;
      for (i = 0; i < qvtrk->ncodes; i++)
        Free_QVcoding(qvtrk->coding+i);
      free(qvtrk->coding);
      free(qvtrk->table);
      fclose(qvtrk->quiva);
      db->tracks = track->next;
      free(track);
    }
  return;
}


/*******************************************************************************************
 *
 *  TRACK LOAD & CLOSE ROUTINES
 *
 ********************************************************************************************/

// If track is not already in the db's track list, then allocate all the storage for it,
//   read it in from the appropriate file, add it to the track list, and return a pointer
//   to the newly created HITS_TRACK record.  If the track does not exist or cannot be
//   opened for some reason, then NULL is returned.

HITS_TRACK *Load_Track(HITS_DB *db, char *track)
{ FILE       *afile, *dfile;
  int         tracklen, size;
  int         nreads;
  void       *anno;
  void       *data;
  HITS_TRACK *record;

  if (track[0] == '.')
    { fprintf(stderr,"Track names cannot begin with a .\n");
      exit (1);
    }

  for (record = db->tracks; record != NULL; record = record->next)
    if (strcmp(record->name,track) == 0)
      return (record);

  afile = Fopen(Catenate(db->path,".",track,".anno"),"r");
  if (afile == NULL)
    return (NULL);
  dfile = Fopen(Catenate(db->path,".",track,".data"),"r");

  fread(&tracklen,sizeof(int),1,afile);
  fread(&size,sizeof(int),1,afile);

  if (db->trimmed)
    { if (tracklen != db->breads)
        { fprintf(stderr,"%s: Track %s not same size as database !\n",Prog_Name,track);
          exit (1);
        }
      if (db->part > 0)
        fseeko(afile,size*db->bfirst,SEEK_CUR);
    }
  else
    { if (tracklen != db->oreads)
        { fprintf(stderr,"%s: Track %s not same size as database !\n",Prog_Name,track);
          exit (1);
        }
      if (db->part > 0)
        fseeko(afile,size*db->ofirst,SEEK_CUR);
    }
  nreads = db->nreads;

  anno = (void *) Malloc(size*(nreads+1),"Allocating Track Anno Vector");

  fread(anno,size,nreads+1,afile);

  if (dfile != NULL)
    { int64 *anno8, off8, dlen;
      int   *anno4, off4;
      int    i;

      if (size == 4)
        { anno4 = (int *) anno;
          off4  = anno4[0];
          if (off4 != 0)
            { for (i = 0; i <= nreads; i++)
                anno4[i] -= off4;
              fseeko(dfile,off4,SEEK_SET);
            }
          dlen = anno4[nreads];
          data = (void *) Malloc(dlen,"Allocating Track Data Vector");
        }
      else
        { anno8 = (int64 *) anno;
          off8  = anno8[0];
          if (off8 != 0)
            { for (i = 0; i <= nreads; i++)
                anno8[i] -= off8;
              fseeko(dfile,off8,SEEK_SET);
            }
          dlen = anno8[nreads];
          data = (void *) Malloc(dlen,"Allocating Track Data Vector");
        }
      fread(data,dlen,1,dfile);
      fclose(dfile);
    }
  else
    data = NULL;

  fclose(afile);

  record = (HITS_TRACK *) Malloc(sizeof(HITS_TRACK),"Allocating Track Record");
  record->name = Strdup(track,"Allocating Track Name");
  record->data = data;
  record->anno = anno;
  record->size = size;

  if (db->tracks != NULL && strcmp(db->tracks->name,".@qvs") == 0)
    { record->next     = db->tracks->next;
      db->tracks->next = record;
    }
  else
    { record->next = db->tracks;
      db->tracks   = record;
    }

  return (record);
}

void Close_Track(HITS_DB *db, char *track)
{ HITS_TRACK *record, *prev;

  prev = NULL;
  for (record = db->tracks; record != NULL; record = record->next)
    { if (strcmp(record->name,track) == 0)
        { free(record->anno);
          free(record->data);
          free(record->name);
          if (prev == NULL)
            db->tracks = record->next;
          else
            prev->next = record->next;
          free(record);
          return;
        }
      prev = record;
    }
  return;
}


/*******************************************************************************************
 *
 *  READ BUFFER ALLOCATION AND READ ACCESS
 *
 ********************************************************************************************/

// Allocate and return a buffer big enough for the largest read in 'db', leaving room
//   for an initial delimiter character

char *New_Read_Buffer(HITS_DB *db)
{ char *read;

  read = (char *) Malloc(db->maxlen+4,"Allocating New Read Buffer");
  if (read == NULL)
    exit (1);
  return (read+1);
}

// Load into 'read' the i'th read in 'db'.  As an upper case ASCII string if ascii is 2, as a
//   lower-case ASCII string is ascii is 1, and as a numeric string over 0(A), 1(C), 2(G), and
//   3(T) otherwise.
//
// **NB**, the byte before read will be set to a delimiter character!

void Load_Read(HITS_DB *db, int i, char *read, int ascii)
{ FILE      *bases  = (FILE *) db->bases;
  int64      off;
  int        len;
  HITS_READ *r = db->reads;

  if (bases == NULL)
    { db->bases = (void *) (bases = Fopen(Catenate(db->path,"","",".bps"),"r"));
      if (bases == NULL)
        exit (1);
    }
  if (i >= db->nreads)
    { fprintf(stderr,"%s: Index out of bounds (Load_Read)\n",Prog_Name);
      exit (1);
    }

  off = r[i].boff;
  len = r[i].end - r[i].beg;

  if (ftello(bases) != off)
    fseeko(bases,off,SEEK_SET);
  fread(read,1,COMPRESSED_LEN(len),bases);
  Uncompress_Read(len,read);
  if (ascii == 1)
    { Lower_Read(read);
      read[-1] = '\0';
    }
  else if (ascii == 2)
    { Upper_Read(read);
      read[-1] = '\0';
    }
  else
    read[-1] = 4;
}


/*******************************************************************************************
 *
 *  QV BUFFER ALLOCATION QV READ ACCESS
 *
 ********************************************************************************************/

// Allocate and return a buffer of 5 vectors big enough for the largest read in 'db'

char **New_QV_Buffer(HITS_DB *db)
{ char **entry;
  char  *qvs;
  int    i;

  qvs   = (char *) Malloc(db->maxlen*5,"Allocating New QV Buffer");
  entry = (char **) Malloc(sizeof(char *)*5,"Allocating New QV Buffer");
  if (qvs == NULL || entry == NULL)
    exit (1);
  for (i = 0; i < 5; i++)
    entry[i] = qvs + i*db->maxlen;
  return (entry);
}

// Load into entry the QV streams for the i'th read from db.  The parameter ascii applies to
//  the DELTAG stream as described for Load_Read.

void Load_QVentry(HITS_DB *db, int i, char **entry, int ascii)
{ HITS_READ *reads;
  FILE      *quiva;
  int        rlen;

  if (db != Active_DB)
    { if (db->tracks == NULL || strcmp(db->tracks->name,".@qvs") != 0)
        { fprintf(stderr,"%s: QV's are not loaded!\n",Prog_Name);
          exit (1);
        }
      Active_QV = (HITS_QV *) db->tracks;
      Active_DB = db;
    }
  if (i >= db->nreads)
    { fprintf(stderr,"%s: Index out of bounds (Load_QVentry)\n",Prog_Name);
      exit (1);
    }

  reads = db->reads;
  quiva = Active_QV->quiva;
  rlen  = reads[i].end-reads[i].beg;

  fseeko(quiva,reads[i].coff,SEEK_SET);
  Uncompress_Next_QVentry(quiva,entry,Active_QV->coding+Active_QV->table[i],rlen);

  if (ascii != 1)
    { char *deltag = entry[1];

      if (ascii != 2)
        { int x = deltag[rlen];
          deltag[rlen] = '\0';
          Number_Read(deltag);
          deltag[rlen] = x;
        }
      else
        { int i;
          int u = 'A'-'a';

          for (i = 0; i < rlen; i++)
            deltag[i] += u;
        }
    }
}


/*******************************************************************************************
 *
 *  BLOCK LOAD OF ALL READS (PRIMARILY FOR DALIGNER)
 *
 ********************************************************************************************/

// Allocate a block big enough for all the uncompressed sequences, read them into it,
//   reset the 'off' in each read record to be its in-memory offset, and set the
//   bases pointer to point at the block after closing the bases file.  If ascii is
//   non-zero then the reads are converted to ACGT ascii, otherwise the reads are left
//   as numeric strings over 0(A), 1(C), 2(G), and 3(T).

void Read_All_Sequences(HITS_DB *db, int ascii)
{ FILE      *bases  = (FILE *) db->bases;
  int        nreads = db->nreads;
  HITS_READ *reads = db->reads;
  void     (*translate)(char *s);

  char  *seq;
  int64  o, off;
  int    i, len;

  if (bases == NULL)
    db->bases = (void *) (bases = Fopen(Catenate(db->path,"","",".bps"),"r"));
  else
    rewind(bases);

  seq = (char *) Malloc(db->totlen+nreads+4,"Allocating All Sequence Reads");

  *seq++ = 4;

  if (ascii == 1)
    translate = Lower_Read;
  else
    translate = Upper_Read;

  o = 0;
  for (i = 0; i < nreads; i++)
    { len = reads[i].end - reads[i].beg;
      off = reads[i].boff;
      if (ftello(bases) != off)
        fseeko(bases,off,SEEK_SET);
      fread(seq+o,1,COMPRESSED_LEN(len),bases);
      Uncompress_Read(len,seq+o);
      if (ascii)
        translate(seq+o);
      reads[i].boff = o;
      o += (len+1);
    }
  reads[nreads].boff = o;

  fclose(bases);

  db->bases  = (void *) seq;
  db->loaded = 1;
}

int List_DB_Files(char *path, void foreach(char *path, char *extension))
{ int            status, rlen, dlen;
  char          *root, *pwd, *name;
  DIR           *dirp;
  struct dirent *dp;

  status = 0;
  pwd    = PathTo(path);
  root   = Root(path,".db");
  rlen   = strlen(root);

  if (root == NULL || pwd == NULL)
    { status = 1;
      goto exit;
    }

  if ((dirp = opendir(pwd)) == NULL)
    { status = 1;
      goto exit;
    }

  while ((dp = readdir(dirp)) != NULL)     //   Get case dependent root name (if necessary)
    { name = dp->d_name;
      if (strcmp(name,Catenate("","",root,".db")) == 0)
        break;
      if (strcasecmp(name,Catenate("","",root,".db")) == 0)
        { strncpy(root,name,rlen);
          break;
        }
    }
  if (dp == NULL)
    { status = 1;
      closedir(dirp);
      goto exit;
    }

  foreach(Catenate(pwd,"/",root,".db"),"db");

  rewinddir(dirp);                         //   Report each auxiliary file
  while ((dp = readdir(dirp)) != NULL)
    { name = dp->d_name;
      dlen = strlen(name);
#ifdef HIDE_FILES
      if (name[0] != '.')
        continue;
      dlen -= 1;
      name += 1;
#endif
      if (dlen < rlen+1)
        continue;
      if (name[rlen] != '.')
        continue;
      if (strncmp(name,root,rlen) != 0)
        continue;
      foreach(Catenate(pwd,PATHSEP,name,""),name+(rlen+1));
    }
  closedir(dirp);

exit:
  free(pwd);
  free(root);
  return (status);
}

void Print_Read(char *s, int width)
{ int i;

  if (s[0] < 4)
    { for (i = 0; s[i] != 4; i++)
        { if (i%width == 0 && i != 0)
            printf("\n");
          printf("%d",s[i]);
        }
      printf("\n");
    }
  else
    { for (i = 0; s[i] != '\0'; i++)
        { if (i%width == 0 && i != 0)
            printf("\n");
          printf("%c",s[i]);
        }
      printf("\n");
    }
}
