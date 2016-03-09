/*** lkupPts: output trace of closest basisPts */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "petcat.h"
#include "neigh.h"

struct cood {
  int run;
  float x, y, z;
} ref_cood[2] = {
  { 111, -27.4, -3.3, 28.5 },
  { 113, -9.6, 0.4, 29.7 }};

int main(int argc, char **argv) {

  FILE *fin, *fou;
  Basis_Point *b, *b_min, *basis;
  char hdr[256];
  struct cood *c;
  char *seg_labels[] = {"n", "l", "r", "u", "d"};
  double dist, dist_min;
  int gridPts, gridSegs, timeSteps, srad, sphi, sz;
  int i, r, n, num;

  if (argc != 2) { fprintf(stderr, "usage: basisPts basisFileName\n"); exit(1); }

  fin = fopen(argv[1], "r");
  if (fin == 0) { fprintf(stderr, "cannot open %s\n", argv[1]); exit(1); }
  fou = fopen("bpt_tr.csv", "w");
  if (fou == 0) { fprintf(stderr, "cannot create output file bpt_tr.csv\n"); exit(1); }

  /* read in header to determine # of basis points - from David's read_basis() */
  num = fread(hdr, 256, 1, fin);
  assert(num == 1);

  num = sscanf(hdr,
          "GRETA basis signals; Cylindrical; version 1.2\n"
          "%d basis point, %d grid segments, %d time steps\n"
          "grid_pos_lu_rpz SRAD SPHI SZ: %d %d %d\n",
           &gridPts, &gridSegs, &timeSteps, &srad, &sphi, &sz);
  assert(num == 6);

  fprintf(stdout, "%d basis points\n", gridPts);
  basis = malloc(sizeof(Basis_Point) * gridPts);
  num = fread(basis, sizeof(Basis_Point), gridPts, fin);
  assert(num == gridPts);

  fprintf(fou, "run,seg,
  ti,v\n");
  for (r = 0; r < 2; r++) {
    c = ref_cood + r;
    dist_min = 100000;
    for (i = 0; i < gridPts; i++) {
      b = basis + i;
      dist = sqrt((b->x - c->x) * (b->x - c->x) +
             (b->y - c->y) * (b->y - c->y) +
             (b->z - c->z) * (b->z - c->z));
      if (dist < dist_min) {
        dist_min = dist;
        b_min = b;
      }
    }
    fprintf(stdout, "closest pt: seg = %d, x=%5.2f, y= %5.2f, z=%5.2f\n",
      b_min->iseg, b_min->x, b_min->y, b_min->z);
    fprintf(stdout, "dx = %5.2f, dy = %5.2f, dz = %5.2f (%f)\n", c->x - b_min->x,
      c->y - b_min->y, c->z - b_min->z, dist_min);
    fflush(stdout);

    for (i = 0; i < 50; i++) {
      fprintf(fou, "%d,%s,%d,%f\n", c->run, seg_labels[0], i,
        b_min->signal[b_min->iseg][i]);
    }
    for (n = 0; n < 4; n++) {
      for (i = 0; i < 50; i++) {
        fprintf(fou, "%d,%s,%d,%f\n", c->run, seg_labels[n + 1], i,
          b_min->signal[neigh[b_min->iseg][n]][i]);
      }
    }

  }

  fclose(fin);
  fclose(fou);

  return 0;
}
