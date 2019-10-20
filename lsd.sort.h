#ifndef LSD_SORT
#define LSD_SORT

void Set_LSD_Params(int nthread, int verbose);

void *LSD_Sort(long long len, void *src, void *trg, int rsize, int dsize, int *bytes);

#endif // LSD_SORT
