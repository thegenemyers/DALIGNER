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

#undef PROFILE

extern int    VERBOSE;      //  -v flag is set?
extern int    MINOVER;      //  minimum overlap (-l)
extern int    HGAP_MIN;     //  HGap minimum (-H)
extern int    SYMMETRIC;    //  output both A vs B and B vs A? ( ! -A)
extern int    IDENTITY;     //  compare reads against themselves?  (-I)
extern int    BRIDGE;       //  bridge consecutive, chainable alignments  (-B)
extern char  *SORT_PATH;    //  where to place temporary files (-P)

extern uint64 MEM_LIMIT;    //  memory limit (-M)
extern uint64 MEM_PHYSICAL;

void Set_Filter_Params(int kmer, int mod, int binshift, int suppress, int hitmin, int nthreads); 

void *Sort_Kmers(DAZZ_DB *block, int *len);

void Match_Filter(char *aname, DAZZ_DB *ablock, char *bname, DAZZ_DB *bblock,
                  void *atable, int alen, void *btable, int blen, Align_Spec *asettings);

void Clean_Exit(int val);

#endif
