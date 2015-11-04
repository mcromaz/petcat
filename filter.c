/*** filter.c: routines for using the detector filter file */

#include <assert.h>
#include <stdio.h>
#include "eb.h"
#include "wrap.h"

struct filter *read_filterfile(char *filename) {

  FILE *fin;
  struct filter *fltr;
  char s[80];
  int num;

  assert(filename != 0);	/* pre */

  fin = fopen(filename, "r");
  if (fin == 0)
    return 0;

  fltr = Malloc(sizeof(struct filter));

  num = fscanf(fin, "%s %f %f\n", s, &(fltr->emin), &(fltr->emax));
  assert(num == 3);
  assert(strcmp(s, "egate") == 0);
  assert(fltr->emin >= 0 && fltr->emax < 16384 && fltr->emin <= fltr->emax);

  num = fscanf(fin, "%s %d %d\n", s, &(fltr->xmin), &(fltr->xmax));
  assert(num == 3);
  assert(strcmp(s, "xgate") == 0);
  assert(fltr->xmin >= 0 && fltr->xmax < 16384 && fltr->xmin <= fltr->xmax);

  num = fscanf(fin, "%s %d %d\n", s, &(fltr->segmin), &(fltr->segmax)); 
  assert(num == 3);
  assert(strcmp(s, "seg") == 0);
  assert(fltr->segmin >= 0 && fltr->segmax < 16384 && fltr->segmin <= fltr->segmax);

  fclose(fin);
  return fltr;
}
