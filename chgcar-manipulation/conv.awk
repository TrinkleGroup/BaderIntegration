/E/{
  N=NF/3;
  for (n=0; n<N; ++n) {
    printf(" %s%s", $(1+3*n), $(2+3*n));
    e = $(3+3*n)+1;
    if (e>=0) printf("+");
    if (e<0) {printf("-"); e=-e}
    printf("%02d", e);
  }
  printf("\n");
  next;}
{print}
