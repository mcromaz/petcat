/*** basisPts: output cood of basis pts to file */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "petcat.h"
#include "neigh.h"

int main(int argc, char **argv) {

  FILE *fin, *fou;
  Basis_Point bpt;
  char hdr[256];
  int gridPts, gridSegs, timeSteps, srad, sphi, sz;
  int localIndex = 1;	/* cnt from 1 */
  int i, j, num;
  int cnt[36], sum = 0;

  (void) memset(cnt, 0, sizeof(cnt));

  if (argc != 2) { fprintf(stderr, "usage: basisPts basisFileName\n"); exit(1); }

  fin = fopen(argv[1], "r");
  if (fin == 0) { fprintf(stderr, "cannot open %s\n", argv[1]); exit(1); }
  fou = fopen("basisPts.csv", "w");
  if (fou == 0) { fprintf(stderr, "cannot create output file basisPts.csv\n"); exit(1); }

  /* read in header to determine # of basis points - from David's read_basis() */
  num = fread(hdr, 256, 1, fin);
  assert(num == 1);

  num = sscanf(hdr,
          "GRETA basis signals; Cylindrical; version 1.2\n"
          "%d basis point, %d grid segments, %d time steps\n"
          "grid_pos_lu_rpz SRAD SPHI SZ: %d %d %d\n",
           &gridPts, &gridSegs, &timeSteps, &srad, &sphi, &sz);
  assert(num == 6);

  fprintf(fou, "i, x, y, z, iseg, ir, ip, iz\n");
  for (i = 0; i < gridPts; i++) {
    num = fread(&bpt, sizeof(Basis_Point), 1, fin);
    assert(num == 1);
    if (bpt.iseg >= 0 && bpt.iseg < 36) { cnt[bpt.iseg]++; }
    fprintf(fou, "%d, %f, %f, %f, %d, %d, %d, %d\n", localIndex, bpt.x, bpt.y, bpt.z, bpt.iseg, bpt.ir, bpt.ip, bpt.iz);
    localIndex++;
  }
  fclose(fin);
  fclose(fou);

  /* consistency check */
  for (i = 0; i < 36; i++) { sum += cnt[i]; }
  assert(gridPts == sum);

  fprintf(stdout, "%s: %d grid pts, (%d, %d, %d)\n", argv[1], gridPts, srad, sphi, sz);
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      fprintf(stdout, "%8d", cnt[6 * i + j]);
    }
    fprintf(stdout, "\n");
  }
  return 0;
}
