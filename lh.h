
#ifndef _lh_h
#define _lh_h

#define LH_MAX_DIM 1048576
#define SPN_LEN 4096
#define MAX_SPE_LEN 16384

typedef struct {
  char *name;
  float *h;
  int dim;
  int cnt;
  int outlyers;
} lh;

lh *lh_init(char *name, int dim);
lh *lh_init_ia(char *name, int *ia, int dim);
void lh_incr(lh *x, int ch);
void lh2spe(lh *x);
void lh2spn(lh *x);

#endif  /* _rad_h */
