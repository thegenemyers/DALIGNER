#ifndef RADIX_SORT
#define RADIX_SORT

void Set_Radix_Params(int nthread, int verbose);

void *LSD_Sort(long long len, void *src, void *trg, int rsize, int dsize, int *bytes);

#endif // RADIX_SORT
