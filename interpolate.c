/* program to decompose event signals for GRETINA / GRETA.

   Author:  D.C. Radford    Aug 2004
*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "petcat.h"

int nearest_grid_points(int seg, double rin, double pin, double zin,
			float *rdiff, float *pdiff, float *zdiff, int pos_out[8])
{
  /* returns 1 on failure */
  int   ir, ip, iz, i, iter;

  if (rin < -1.0           || pin < -1.0           || zin < -1.0 ||
      rin > maxir[seg] + 1 || pin > maxip[seg] + 1 || zin > maxiz[seg] + 1) {
    /* printf("=== %d %f %f %f\n", seg, rin, pin, zin); */
    return 1;
  }

  ir = rin + 0.0001;
  ip = pin + 0.0001;
  iz = zin + 0.0001;
  if (ir<0) ir = 0;
  if (ip<0) ip = 0;
  if (iz<0) iz = 0;
  if (ir>maxir[seg] - 1) ir = maxir[seg] - 1;
  if (ip>maxip[seg] - 1) ip = maxip[seg] - 1;
  if (iz>maxiz[seg] - 1) iz = maxiz[seg] - 1;

  for (iter=0; iter<10; iter++) {
    pos_out[0] = grid_pos_lu[seg][ir  ][ip  ][iz  ];
    pos_out[1] = grid_pos_lu[seg][ir  ][ip  ][iz+1];
    pos_out[2] = grid_pos_lu[seg][ir  ][ip+1][iz  ];
    pos_out[3] = grid_pos_lu[seg][ir  ][ip+1][iz+1];
    pos_out[4] = grid_pos_lu[seg][ir+1][ip  ][iz  ];
    pos_out[5] = grid_pos_lu[seg][ir+1][ip  ][iz+1];
    pos_out[6] = grid_pos_lu[seg][ir+1][ip+1][iz  ];
    pos_out[7] = grid_pos_lu[seg][ir+1][ip+1][iz+1];

    for (i=0; i<8; i++) {
      if (pos_out[i] < 0) {
	if (iter == 2) return 1;
	/* printf("*>*>* %d %d %d %d <<< %d %f %f %f ", i, ir, ip, iz, seg, rin, pin, zin); */
	if (ir > maxir[seg]/2) {
	  ir-= 2;
	} else {
	  ir+= 2;
	}
	if (iz > maxiz[seg]/2) {
	  iz-= 2;
	} else {
	  iz+= 2;
	}
	/* printf(">>> %d %d %d %d\n", i, ir, ip, iz); */
	break;
      }
    }
    if (i==8) break;
  }
  *rdiff = rin - (double) ir;
  *pdiff = pin - (double) ip;
  *zdiff = zin - (double) iz;

  return 0;
} /* nearest_grid_points */

/* ================================================================ */

int interpolate(int seg, double rin, double pin, double zin, Basis_Point *signal,
	        Basis_Point sigderiv[3], int calc_deriv, int *ssel)
{
  /*
    basis-point-interpolation routine
    Given a position (rin, pin, zin), obtain a set of eight basis
    point positions pos[8], which are the nearest neighbors.
    Use them to calculate the interpolated signal *sig,
    and, if (calc_deriv), the partial derivatives *sigderiv[3]
    with respect to r,p,z.
    (rdiff, pdiff, zdiff) is the position of the point relative
    to pos[0].
  */
  /*
    returns 1 on failure
     - for a Cartesian basis, this is usually because there are not
       enough valid neighbors to interpolate or calculate the derivatives
  */

  /* double sum; */
  float  rdiff, pdiff, zdiff;
  float  d, ddr, ddp, ddz, dist[8], dr[8], dp[8], dz[8];
  int    i, j, k, nneigh, id[8], pos[8], posi;
  float  *sig, *sd1, *sd2, *sd3, *bsig;

  /* get indices of eight nearest neighbor points */
  if (nearest_grid_points(seg, rin, pin, zin, &rdiff, &pdiff, &zdiff, pos))
    return 1;

  sig = signal->signal[0];
  sd1 = sigderiv[0].signal[0];
  sd2 = sigderiv[1].signal[0];
  sd3 = sigderiv[2].signal[0];

  /* weighting factors to calculate sig from the 8 signals */
  dist[0] = (1.0f-rdiff)*(1.0f-pdiff)*(1.0f-zdiff);
  dist[1] = (1.0f-rdiff)*(1.0f-pdiff)*      zdiff;
  dist[2] = (1.0f-rdiff)*      pdiff *(1.0f-zdiff);
  dist[3] = (1.0f-rdiff)*      pdiff *      zdiff;
  dist[4] =       rdiff *(1.0f-pdiff)*(1.0f-zdiff);
  dist[5] =       rdiff *(1.0f-pdiff)*      zdiff;
  dist[6] =       rdiff *      pdiff *(1.0f-zdiff);
  dist[7] =       rdiff *      pdiff *      zdiff;

  if (calc_deriv) {
    /* weighting factors to calculate sigderiv from the 8 signals */
    dr[0]   = -(1.0f-pdiff)*(1.0f-zdiff);
    dr[1]   = -(1.0f-pdiff)*      zdiff;
    dr[2]   =       -pdiff *(1.0f-zdiff);
    dr[3]   =       -pdiff *      zdiff;
    dr[4]   =  (1.0f-pdiff)*(1.0f-zdiff);
    dr[5]   =  (1.0f-pdiff)*      zdiff;
    dr[6]   =        pdiff *(1.0f-zdiff);
    dr[7]   =        pdiff *      zdiff;

    dp[0]   = -(1.0f-rdiff)*(1.0f-zdiff);
    dp[1]   = -(1.0f-rdiff)*      zdiff;
    dp[2]   =  (1.0f-rdiff)*(1.0f-zdiff);
    dp[3]   =  (1.0f-rdiff)*      zdiff;
    dp[4]   =       -rdiff *(1.0f-zdiff);
    dp[5]   =       -rdiff *      zdiff;
    dp[6]   =        rdiff *(1.0f-zdiff);
    dp[7]   =        rdiff *      zdiff;

    dz[0]   = -(1.0f-rdiff)*(1.0f-pdiff);
    dz[1]   =  (1.0f-rdiff)*(1.0f-pdiff);
    dz[2]   = -(1.0f-rdiff)*      pdiff;
    dz[3]   =  (1.0f-rdiff)*      pdiff;
    dz[4]   =       -rdiff *(1.0f-pdiff);
    dz[5]   =        rdiff *(1.0f-pdiff);
    dz[6]   =       -rdiff *      pdiff;
    dz[7]   =        rdiff *      pdiff;
  }

  /* count the actual number of VALID neighbors, and compute
     the minimum and maximum interesting time range for each segment */
  for (j=0; j<GRID_SEGS; j++) {
    signal->lo_time[j] = TIME_STEPS;
    signal->hi_time[j] = 0;
  }
  nneigh = 0;
  for (i=0; i<8; i++) {
    if ((posi = pos[i]) >= 0) {
      id[nneigh++] = i;
      for (j=0; j<GRID_SEGS; j++) {
	if (signal->lo_time[j] > basis[posi].lo_time[j])
	  signal->lo_time[j] = basis[posi].lo_time[j];
	if (signal->hi_time[j] < basis[posi].hi_time[j])
	  signal->hi_time[j] = basis[posi].hi_time[j];
      }
    }
  }

  /* if there are not enough valid neighbors, we need to quit */
  /* if (nneigh == 0 || (calc_deriv>0 && nneigh<7)) { */
  if (nneigh == 0) {
    printf("*** No neighbors! seg r p z: %d %.2f %.2f %.2f\n", seg, rin, pin, zin);
    if (!calc_deriv) {
      for (i=0; i<GRID_SEGS*TIME_STEPS; i++) {
	sig[i] = 0.0f;
      }
    } else {
      for (i=0; i<GRID_SEGS*TIME_STEPS; i++) {
	sig[i] = 0.0f;
	sd1[i] = 0.0f;
	sd2[i] = 0.0f;
	sd3[i] = 0.0f;
      }
    }
    return 1;
  }

  /*
    sum = dist[id[0]);
    for (i=1; i<nneigh; i++) {
    sum = sum + dist[id[i]);
    }
    if (sum==0.0) return 1;
  */

  if (!calc_deriv) {
    /* compute only the interpolated signal */
    for (i=0; i<GRID_SEGS*TIME_STEPS; i++) {
      sig[i] = 0.0f;
    }
    for (k=0; k<nneigh; k++) {
      posi = pos[id[k]];
      d = dist[id[k]];
      for (i=0; i<GRID_SEGS; i++) {
	if (!ssel[i]) continue;
	sig = signal->signal[i];
	bsig = basis[posi].signal[i];
	for (j=signal->lo_time[i]; j<=signal->hi_time[i]; j++) {
	  sig[j] += bsig[j]*d;
	}
      }
    }
    /*
      for (i=0; i<GRID_SEGS; i++) {
      for (j=signal->lo_time[i]; j<=signal->hi_time[i]; j++) {
      sig->signal[i][j] /= sum;
      }
      }
    */

  } else {
    /* compute both the interpolated signal
       and the partial derivatives */
    for (i=0; i<GRID_SEGS*TIME_STEPS; i++) {
      sig[i] = 0.0f;
      sd1[i] = 0.0f;
      sd2[i] = 0.0f;
      sd3[i] = 0.0f;
    }
    for (k=0; k<nneigh; k++) {
      posi = pos[id[k]];
      d = dist[id[k]];
      ddr = dr[id[k]];
      ddp = dp[id[k]];
      ddz = dz[id[k]];
      for (i=0; i<GRID_SEGS; i++) {
	if (!ssel[i]) continue;
	sig = signal->signal[i];
	sd1 = sigderiv[0].signal[i];
	sd2 = sigderiv[1].signal[i];
	sd3 = sigderiv[2].signal[i];
	bsig = basis[posi].signal[i];
	for (j=signal->lo_time[i]; j<=signal->hi_time[i]; j++) {
	  sig[j] += bsig[j]*d;
	  sd1[j] += bsig[j]*ddr;
	  sd2[j] += bsig[j]*ddp;
	  sd3[j] += bsig[j]*ddz;
	}
      }
    }
    /*
      for (i=0; i<GRID_SEGS; i++) {
      if (!ssel[i]) continue;
      for (j=signal->lo_time[i]; j<=signal->hi_time[i]; j++) {
      sig->signal[i][j] /= sum;
      sd1->signal[i][j] /= sum;
      sd2->signal[i][j] /= sum;
      sd3->signal[i][j] /= sum;
      }
      }
    */
  }

  return 0;
} /* interpolate */

/* ================================================================ */

int cyl_to_cart(int seg, double *pars, double*x, double *y, double*z, double *e)
{
  float  rdiff, pdiff, zdiff, dist[8];
  int    j, ir, ip, iz, pos[8];

  *x = *y = *z = 0;

  /* get indices of eight nearest neighbor points */
  if (!nearest_grid_points(seg, pars[0], pars[1], pars[2], &rdiff, &pdiff, &zdiff, pos)) {
    /* weighting factors */
    dist[0] = (1.0f-rdiff)*(1.0f-pdiff)*(1.0f-zdiff);
    dist[1] = (1.0f-rdiff)*(1.0f-pdiff)*      zdiff;
    dist[2] = (1.0f-rdiff)*      pdiff *(1.0f-zdiff);
    dist[3] = (1.0f-rdiff)*      pdiff *      zdiff;
    dist[4] =       rdiff *(1.0f-pdiff)*(1.0f-zdiff);
    dist[5] =       rdiff *(1.0f-pdiff)*      zdiff;
    dist[6] =       rdiff *      pdiff *(1.0f-zdiff);
    dist[7] =       rdiff *      pdiff *      zdiff;
    for (j=0; j<8; j++) {
      if (pos[j] >= 0) {
	*x += dist[j]*basis[pos[j]].x;
	*y += dist[j]*basis[pos[j]].y;
	*z += dist[j]*basis[pos[j]].z;
      }
    }
  } else {
    ir = pars[0] + 0.0001;
    ip = pars[1] + 0.0001;
    iz = pars[2] + 0.0001;
    if (ir<0) ir = 0;
    if (ip<0) ip = 0;
    if (iz<0) iz = 0;
    if (ir>maxir[seg]) ir = maxir[seg];
    if (ip>maxip[seg]) ip = maxip[seg];
    if (iz>maxiz[seg]) iz = maxiz[seg];
    j = grid_pos_lu[seg][ir][ip][iz];
    if (j < 0) return 1;
    *x = basis[j].x;
    *y = basis[j].y;
    *z = basis[j].z;
  }
  *e = pars[3];
  return 0;
} /* cyl_to_cart */
