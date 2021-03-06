#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdlib.h>

#include "petcat.h"
#include "wrap.h"

typedef struct {
  char  iseg, ir, ip, iz;             /* integer cylindrical coordinates of grid point */
  float x, y, z;                      /* cartesian coordinates of grid point */
  float signal[MAX_GRID_SEGS][50];        /* actual basis signals */
  int   lo_time[MAX_GRID_SEGS], hi_time[MAX_GRID_SEGS];     /* limits for non-zero signal */
} Basis_Point_50;

int read_basis(char *basisfile)
{
/* routine to read decomposition basis signal files from
   .dat unformatted binary file BASIS_FILE (defined in gdecomp.h)
   This file would normally have been created by
   the program convert_basis_sig2dat

   returns 0 on success, 1 on failure
   modifies:       Basis_Point basis[GRID_PTS];
                   int         grid_pos_lu[SX][SY][SZ];
      defined in gdecomp.h.

   Author:  D.C. Radford    Aug 2004
*/
  char  header[256];
  int   i, ii, j, k, t, s, nscanned;
  int r,p,z;
  FILE  *file;
  Basis_Point_50 *tmp;
  Basis_Point *bpt;
  int num;
  static int no_bad_segs[1] = {-1};

  tmp = malloc(sizeof(Basis_Point_50));

  bad_segs = no_bad_segs;
  /* if (set_bad_segs()) bad_segs = no_bad_segs; */

  if (basis) free(basis);


  printf("Reading basis signals from %s\n", basisfile);
  if (!(file=fopen(basisfile, "r"))) {
    printf("ERROR -- Cannot open file %s for reading.\n", basisfile);
    return 1;
  }

  fread(header, 256, 1, file);
   nscanned = sscanf(header,
          "GRETA basis signals; Cylindrical; version 1.2\n"
          "%d basis point, %d grid segments, %d time steps\n"
          "grid_pos_lu_rpz SRAD SPHI SZ: %d %d %d\n",
           &GRID_PTS, &GRID_SEGS, &TIME_STEPS, &SRAD, &SPHI, &SZZZ);

   if (nscanned != 6) {
      printf("ERROR reading basis file parameters\n");
      fclose(file);
      return 1;
   }
   /* Note that GRID SEGS is not truly variable. If it must be, then the basis
      read will have to expanded similarly (but with more complexity) to the
      lookup table read.
    */

#ifdef SAMPLE25
  printf("*** SAMPLE25 defined ***\n");
  TIME_STEPS = 25; // redefine value in basis header
#endif

   if (GRID_PTS > MAX_GRID_PTS || GRID_SEGS != MAX_GRID_SEGS ||
       TIME_STEPS > MAX_TIME_STEPS || SRAD > MAX_SRAD ||
       SPHI > MAX_SPHI || SZZZ > MAX_SZZZ) {
      printf("ERROR: basis file parameter out of range\n");
      fclose(file);
      return 1;
   }

   /* malloc the space for the basis */
   if (!(basis = (Basis_Point*) malloc(sizeof(Basis_Point) * GRID_PTS))) {
     printf("\nERROR  -  cannot malloc basis!\n\n");
     exit(-1);
   }

   /* read basis itself */
#ifdef SAMPLE25
    /* read to intermediate 50 sample basis point and compress */
    for (i = 0; i < GRID_PTS; i++) {
      num = fread(tmp, sizeof(Basis_Point_50), 1, file);
      if (num != 1) {
        printf("Something's wrong with the basis data file!\n");
        fclose(file);
        return 1;
      }
      bpt = basis + i;
      memset(bpt, 0, sizeof(Basis_Point));
      bpt->iseg = tmp->iseg;
      bpt->ir = tmp->ir;
      bpt->ip = tmp->ip;
      bpt->iz = tmp->iz;
      bpt->x = tmp->x;
      bpt->y = tmp->y;
      bpt->z = tmp->z;
      for (j = 0; j < MAX_GRID_SEGS; j++) {
        bpt->lo_time[j] = tmp->lo_time[j] / 2;
        bpt->hi_time[j] = tmp->hi_time[j] / 2;
        for (k = 0; k < 25; k++) {
          bpt->signal[j][k] = (tmp->signal[j][2 * k] + tmp->signal[j][2 * k + 1]) / 2;
      }
    }
  }
#else
    if (fread(basis, sizeof(Basis_Point), GRID_PTS, file) != GRID_PTS) {
    /* something's wrong */
    printf("Something's wrong with the basis data file!\n");
    fclose(file);
    return 1;
  }
#endif

  /* read lookup table  */
  for (s = 0; s < SSEG; s++) {
     for (r = 0; r < SRAD; r++) {
         for (p = 0; p < SPHI; p++) {
            for (z = 0; z < SZZZ; z++ ) {
                if (fread(&grid_pos_lu[s][r][p][z],sizeof(int),1,file)!=1) {
                   printf("Error reading lookup table at %d %d %d %d\n",
                                s, r, p, z);
                   fclose(file);
                   return 1;
                }
            }
         }
      }
   }

  fclose(file);

  /* set bad segment signals to zero */
  for (i=0; bad_segs[i]>=0; i++) {
    ii = bad_segs[i];
    for (j=0; j<GRID_PTS; j++) {
      for (t=0; t<TIME_STEPS; t++) basis[j].signal[ii][t] = 0.0f;
      basis[j].lo_time[ii] = TIME_STEPS;
      basis[j].hi_time[ii] = 0;
    }
  }

  /* find the maximum values of ir, ip and iz for each segment */
  for (i=0; i<SSEG; i++) {
    maxir[i] = maxip[i] = maxiz[i] = 0;
  }
  for (i=0; i<GRID_PTS; i++) {
    s = basis[i].iseg;
    if (maxir[s] < basis[i].ir) maxir[s] = basis[i].ir;
    if (maxip[s] < basis[i].ip) maxip[s] = basis[i].ip;
    if (maxiz[s] < basis[i].iz) maxiz[s] = basis[i].iz;
  }
  if (!quiet) {
    printf(" seg  max.ir  max.ip  max.iz\n");
    for (i=0; i<SSEG; i++) {
      printf(" %3d %7d %7d %7d\n", i, maxir[i], maxip[i], maxiz[i]);
    }
    printf(" All %7d %7d %7d\n\n", SRAD-1, SPHI-1, SZZZ-1);
  }

  return 0;
}
