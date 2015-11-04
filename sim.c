/*** sim: run simulated data through decomp w/o preprocessing */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "petcat.h"


int numGridPoints(char *name) {

  FILE *fin;
  char hdr[256];
  int grid_pts, grid_segs, time_steps, srad, sphi, sz;
  int num;

  fin = fopen(name, "r");
  assert(fin != 0);
  fread(hdr, 256, 1, fin);
  num = sscanf(hdr, "GRETA basis signals; Cylindrical; version 1.2\n" "%d basis point, %d grid segments, %d time steps\n" "grid_pos_lu_rpz SRAD SPHI SZ: %d %d %d\n", &grid_pts, &grid_segs, &time_steps, &srad, &sphi, &sz);
  assert(num == 6);
  fclose(fin);

  return grid_pts;
}

Basis_Point *nearestPt(float x0, float y0, float z0, int gridPts) {

  float dx, dy, dz;
  float d2, d2min;
  int i, imin;
 
  d2min = 100000.;
  imin = 0;
  for (i = 0; i < gridPts; i++) {
    dx = basis[i].x - x0;    
    dy = basis[i].y - y0;    
    dz = basis[i].z - z0;    
    d2 = dx * dx + dy * dy + dz * dz;
    if (d2 < d2min) {
      d2min = d2;
      imin = i;
    }
  }

  return basis + imin;
}

int main() {

  Event_Signal e;
  struct decomp_state *a;
  struct crys_intpts *x;
  struct postprocCnt *postCnt;
  char *basisFile = "basis_xt.dat"; 
  int gridPts;
  Basis_Point *b1, *b2;
  Basis_Point *b[16];
  float eint[16];
  FILE *fsim, *fsimlog;
  float e0, x0, y0, z0;
  int id, num, numInt, i, j, k;
  char s[256];

  /* read in interacton poistions */
  /* lookup signals in basis */
  /* load asig structure */

  gridPts = numGridPoints(basisFile);

  a = dl_decomp_init(basisFile, 1); // 1 to suppress diag info
  assert(a != 0);

  fsim = fopen("outputA.txt", "r");
  assert(fsim != 0);

  fsimlog = fopen("simpts.txt", "w");
  assert(fsimlog != 0);

  /* read in header */
  numInt = 0;
  if (fgets(s, 80, fsim) == 0) {
      exit(0);
  } 
  while (1) { /* loop over events */
    printf("new event!\n");
    while (1) { /* loop over int pts */ 
      if (fgets(s, 80, fsim) == 0) {
        exit(0);
      } 
      num = sscanf(s, "%d %f %*f %*f %*f %f %f %f %*d", &id, &e0, &x0, &y0, &z0);
      if (id == -1) {	/* end evt, process */
        /* load basis points into e */
        assert(numInt >= 1 && numInt < 16); /* pre */
        memset(&e, 0, sizeof(e));
        for (i = 0; i < numInt; i++) {
          e.seg_energy[b[i]->iseg] += eint[i];
          e.total_energy += eint[i];
        }
        printf("numInt = %d\n", numInt);
        for (i = 0; i < numInt; i++) {
          for (j = 0; j < 37; j++) {
            for (k = 0; k < 50; k++) {
              e.signal[j][k] += (e.seg_energy[b[i]->iseg] / e.total_energy)  * b[i]->signal[j][k];
            }
          }
        }
        /* decompose e */
        x = dl_decomp(a, &e, postCnt);
        log2intpts(x);

        /* reset */
        numInt = 0;
        break;
      }
      assert(numInt <= 16);
      eint[numInt] = e0;
      b[numInt] = nearestPt(x0, y0, z0, gridPts);
      numInt++;
      logintpt(x0, y0, z0, e0, fsimlog);
      printf("id = %d, eint = %f, x = %f, y = %f, z = %f\n", id, e0, x0, y0, z0);
    }
  }

/*
  b1 = nearestPt(-18.0, 5.0, 33.0, gridPts);
  printf("seg = %d, x = %f, y = %f, z = %f\n", b1->iseg, b1->x, b1->y, b1->z);
  b2 = nearestPt(-15.0, 5.0, 36.0, gridPts);
  printf("seg = %d, x = %f, y = %f, z = %f\n", b2->iseg, b2->x, b2->y, b2->z);
*/

  /* cp b[] to e */
 memset(&e, 0, sizeof(e));
 e.total_energy = 1000.;
 e.seg_energy[b1->iseg] = 1000.;
 for (i = 0; i < 37; i++) {
    for (j = 0; j < 50; j++) {
      e.signal[i][j] = 0.666 * b1->signal[i][j] + 0.333 * b2->signal[i][j];
    }
  }
 x = dl_decomp(a, &e, postCnt);
  log2intpts(x);
}
