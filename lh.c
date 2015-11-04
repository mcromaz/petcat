/*** lh: lightweight histogrammer */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "lh.h"
#include "wrap.h"

struct spePrefix {
  int reclA;
  unsigned int titleA;
  unsigned int titleB;
  int dim;
  int a1;
  int a2;
  int a3;
  int reclB;
} spePre = {24, 0, 0, 0, 1, 1, 1, 24};

lh *lh_init(char *name, int dim) {
  lh *x;

  assert(name != 0);			/* pre */
  assert(dim > 0 && dim < LH_MAX_DIM);

  x = Calloc(1, sizeof(lh));
  x->name = name;
  x->h = Calloc(dim, sizeof(float));
  x->dim = dim;
  x->cnt = 0, x->outlyers = 0;
  return x;
}

lh *lh_init_ia(char *name, int *h, int dim) {

  lh *x;
  int i;

  x = Calloc(1, sizeof(lh));
  x->name = name;
  x->h = Calloc(dim, sizeof(float));
  x->dim = dim;
  for (i = 0; i < x->dim; i++) {
    x->h[i] = (float) h[i];
    x->cnt += h[i];
  }
  x->outlyers = 0;
  return x;
}

void lh_incr(lh *x, int ch) {
  
  assert(x != 0);	/* pre */

  if (ch >= 0 && ch < x->dim) {
    x->h[ch] += 1.0;
    x->cnt += 1;
  }
  else {
    x->outlyers += 1;
  }
  return;
}

void lh2spe(lh *x) {

  FILE *fp_spectra;
  char *filename, *rootname;
  int recl;
  int i;

  assert(x != 0);		/* pre */
  assert(x->dim > 0 && x->dim <= 16384);

  rootname = x->name;
  filename = Malloc(strlen(rootname) + 4);
  sprintf(filename, "%s.spe", rootname);
  
  fp_spectra = fopen(filename, "w");
  if (fp_spectra == NULL) {
    fprintf(stderr, "unable to open file %s for spectra", filename);
    return;
  }

  spePre.dim = x->dim;
  recl = sizeof(float) * x->dim;

  (void) fwrite(&spePre, sizeof(struct spePrefix), 1, fp_spectra);
  (void) fwrite(&recl, sizeof(int), 1, fp_spectra);
  (void) fwrite(x->h, sizeof(float), x->dim, fp_spectra);
  (void) fwrite(&recl, sizeof(int), 1, fp_spectra);

  fclose(fp_spectra);
}


void lh2spn(lh *x)  {
  ;
}
