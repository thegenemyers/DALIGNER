/*******************************************************************************************
 *
 *  Filter interface for the dazzler.
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
 *
 ********************************************************************************************/

#ifndef _FILTER

#define _FILTER

#include "DB.h"
#include "align.h"

extern int BIASED;
extern int VERBOSE;
extern int MINOVER;
extern int HGAP_MIN;

#define NTHREADS  4    //  Must be a power of 2
#define NSHIFT    2    //  log_2 NTHREADS

int Set_Filter_Params(int kmer, int binshift, int suppress, int hitmin); 

void Build_Table(HITS_DB *block);

void Match_Filter(char *aname, HITS_DB *ablock, char *bname, HITS_DB *bblock,
                  int self, int comp, Align_Spec *asettings);

#endif
