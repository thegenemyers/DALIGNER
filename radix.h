#ifndef RADIX_SORT
#define RADIX_SORT

void Set_Radix_Params(int nthread, int verbose);

void *Radix_Sort(long long len, void *src, void *trg, int *bytes);

#endif // RADIX_SORT
