/*** log2: routines supporting diagnostics for mode2 data */

#include <stdio.h>
#include <assert.h>
#include "petcat.h"

static FILE *flog;

void log2intpts(struct crys_intpts *a) {

  static int init = 0;
  int i;

  if (!init) {
    flog = fopen("intpts.txt", "w");
    assert(flog != 0);
    init = 1;
  }

  for (i = 0; i < a->num; i++) {
    fprintf(flog, "%f, %f, %f, %f\n", (a->intpts)[i].x, (a->intpts)[i].y, (a->intpts)[i].z, (a->intpts)[i].e);
  }
  //fprintf(flog, "\n");
}

void logintpt(float x, float y, float z, float e, FILE *flog) {

  fprintf(flog, "%f, %f, %f, %f\n", x, y, z, e); 
}
