/* program to decompose event signals for GRETINA / GRETA.

   Author:  D.C. Radford    Aug 2004
*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "petcat.h"


int eval(const Event_Signal *asig,  /* observed signals */
	 Event_Signal *bsig, /* fitted signals */
	 Interaction *ints, double t0, int nints,
	 double *chisq_out, double beta[MAX_PARS],
	 double alpha[MAX_PARS][MAX_PARS], int calc_deriv, int *ssel)
{
  /* returns 1 on failure */

  double rp, pp, zp, e[MAX_PARS/4], chisq;
  float  dif, *sigd[MAX_PARS], *sig1, *sig2;
  int    i, j, k, ii, kk, m, n, tt, npars;

  Basis_Point sigderiv[MAX_PARS], bsig1;


  npars = 4 * nints;
  for (i=0; i<MEAS_SEGS; i++) {
    for (j=0; j<TIME_STEPS; j++) {
      bsig->signal[i][j] = 0.0;
    }
  }
  tt = t0;
  t0 -= tt;

  /* loop over interactions; each one uses four parameters */
  for (i=0; i<nints; i++) {
    ii   = 4*i;
    rp   = ints[i].r;
    pp   = ints[i].p;
    zp   = ints[i].z;
    e[i] = ints[i].e;

    /* find the nearest neighbor grid points, and interpolate to get
       the signal and derivatives at (rp, pp, zp) */
    if (interpolate(ints[i].seg, rp, pp, zp, &bsig1, sigderiv+ii, calc_deriv, ssel)) return 1;
    /* printf("eval:"); for (j=0;j<8;j++) printf(" %d", pos[j]); printf("\n"); */
    
    /* add up the total signal and derivatives */
    if (!calc_deriv) {
      for (j=0; j<GRID_SEGS; j++) {
	if (!ssel[j]) continue;
	for (k=0; k<TIME_STEPS; k++) {
	  bsig->signal[j][k] += bsig1.signal[j][k]*e[i];
	}
      }
    } else {
      for (j=0; j<GRID_SEGS; j++) {
	if (!ssel[j]) continue;
	for (m=ii; m<ii+4; m++) {
	  sigd[m] = sigderiv[m].signal[j];
	}
	sigderiv[ii].lo_time[j] = bsig1.lo_time[j];
	sigderiv[ii].hi_time[j] = bsig1.hi_time[j];
	for (k=0; k<TIME_STEPS; k++) {
	  bsig->signal[j][k] += bsig1.signal[j][k]*e[i];
	  sigd[ii  ][k] *= e[i];
	  sigd[ii+1][k] *= e[i];
	  sigd[ii+2][k] *= e[i];
	  sigd[ii+3][k] = bsig1.signal[j][k];
	}
      }
    }
  }

  /* adjust signals and derivatives for t-zero */
  if (tt > 0) {
    for (j=0; j<GRID_SEGS; j++) {
      if (!ssel[j]) continue;
      sig2 = bsig->signal[j];
      for (k=TIME_STEPS-1; k>=tt; k--) {
	sig2[k] = sig2[k-tt];
      }
      for (k=0; k<tt; k++) {
	sig2[k] = 0;
      }
      if (calc_deriv) {
	for (i=0; i<npars; i++) {
	  sig2 = sigderiv[i].signal[j];
	  /* sigderiv[i].lo_time[j] -= tt;  */ /*CHECKME*/
	  /* if (sigderiv[i].lo_time[j] < 0) sigderiv[i].lo_time[j] = 0;  */ /*CHECKME*/
	  sigderiv[i].lo_time[j] += tt;  /* CHECKME */
	  sigderiv[i].hi_time[j] += tt;  /* CHECKME */
	  if (sigderiv[i].hi_time[j] > TIME_STEPS-1)
	    sigderiv[i].hi_time[j] = TIME_STEPS-1;
	  for (k=TIME_STEPS-1; k>=tt; k--) {
	    sig2[k] = sig2[k-tt];
	  }
	  for (k=0; k<tt; k++) {
	    sig2[k] = 0;
	  }
	}
      }
    }
  }
  for (j=0; j<GRID_SEGS; j++) {
    if (!ssel[j]) continue;
    /* calculate derivative w.r.t. t-zero */
    sigderiv[npars].lo_time[j] = 0;
    sigderiv[npars].hi_time[j] = TIME_STEPS-1;
    sig1 = sigderiv[npars].signal[j];
    sig2 = bsig->signal[j];
    sig1[0] = 0;
    if (t0 < 0.0001) {
      for (k=TIME_STEPS-1; k>0; k--) {
	sig1[k] = sig2[k-1] - sig2[k];
      }
    } else {
      /* interpolate signals and derivs for fractional value of t-zero */
      for (k=TIME_STEPS-1; k>0; k--) {
	sig1[k] = sig2[k-1] - sig2[k];
	sig2[k] += t0 * sig1[k];
      }
      if (calc_deriv) {
	for (i=0; i<npars; i++) {
	  sig2 = sigderiv[i].signal[j];
	  for (k=TIME_STEPS-1; k>0; k--) {
	    sig2[k] += t0 * (sig2[k-1]-sig2[k]);
	  }
	}
      }
    }
  }

  for (j=0; j<MEAS_SEGS; j++) {
    for (k=0; k<TIME_STEPS; k++) {
      // printf("basis: %d %d = %f %f\n", j, k, asig->signal[j][k], bsig->signal[j][k]);
    }
  }
  // exit(0);

  /* calculate chisq */
  chisq = 0.0;
  for (i=0; i<MEAS_SEGS; i++) {
    if (!ssel[i]) continue;
    for (j=0; j<TIME_STEPS; j++) {
      dif = bsig->signal[i][j] - asig->signal[i][j];
      chisq += (double) (dif*dif);
    }
  }

  /* if necessary, calculate the beta and alpha matrices */
  if (calc_deriv) {
    for (k=0; k<npars; k++) {
      beta[k] = 0.0;
    }
    for (i=0; i<MEAS_SEGS; i++) {
      if (!ssel[i]) continue;
      for (m=0; m<npars+1; m++) {
	ii = 4*(m/4);
	sigd[m] = sigderiv[m].signal[i];
	kk = sigderiv[ii].hi_time[i];
	for (k=sigderiv[ii].lo_time[i]; k<=kk; k++) {
	  /* dif  = bsig->signal[i][k] - asig->signal[i][k]; */
	  dif  = asig->signal[i][k] - bsig->signal[i][k];
	  beta[m] = beta[m] + (double) (dif * sigd[m][k]);
	}
	for (n=0; n<=m; n++) {
	  for (k=sigderiv[ii].lo_time[i]; k<=kk; k++) {
	    alpha[m][n] += (double) (sigd[m][k] * sigd[n][k]);
	  }
	}
      }
    }
  }

  *chisq_out = chisq;
  return 0;
} /* eval */

/* ================================================================ */

double eval_int_pos(const Event_Signal *asig,  /* observed signals */
		    Event_Signal *bsig, /* fitted signals */
		    Interaction *ints, double t0, int nints)
{
  /* evaluate fit and chi-squared, using integer positions
     for the interactions when possible.
     NOTE that this eval does NOT evaluate the derivatives */

  /* returns chisq; 1000.0 on failure */

  static int all[37] =
    {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

  double rp, pp, zp, e, chisq = 0.0;
  float  dif, *sig, *sig1, *sig2;
  int    ir, ip, iz, i, j, k, ii, kk, pos, tt, npars, iseg;

  Basis_Point sigderiv, bsig1;

  npars = 4 * nints;
  for (i=0; i<MEAS_SEGS; i++) {
    for (j=0; j<TIME_STEPS; j++) {
      bsig->signal[i][j] = 0.0;
    }
  }

  /* loop over interactions; each one uses four parameters */
  for (i=0; i<npars/4; i++) {
    ii = 4*i;
    rp = ints[i].r;
    pp = ints[i].p;
    zp = ints[i].z;
    e  = ints[i].e;
    ir = rp + 0.5;
    ip = pp + 0.5;
    iz = zp + 0.5;
    iseg  = ints[i].seg;
    if (ir>=0        && ip>=0        && iz>=0 &&
	ir<=maxir[iseg] && ip<=maxip[iseg] && iz<=maxiz[iseg] &&
	fabs(rp - (double) ir) < 0.01 &&
	fabs(pp - (double) ip) < 0.01 &&
	fabs(zp - (double) iz) < 0.01 &&
	(pos = grid_pos_lu[iseg][ir][ip][iz]) >= 0) {

      /* (r,p,z) are integers; no need to interpolate */
      /* add up the total signal */
      for (j=0; j<GRID_SEGS; j++) {
	sig = bsig->signal[j];
	sig1 = basis[pos].signal[j];
	kk = basis[pos].hi_time[j];
	for (k=basis[pos].lo_time[j]; k<=kk; k++) {
	  sig[k] += sig1[k]*e;
	}
      }

    } else {

      /* find the nearest neighbor grid points, and interpolate to get
	 the signal and derivatives at (rp, pp, zp) */
      if (interpolate(iseg, rp, pp, zp, &bsig1, &sigderiv, 0, all)) return 1000.0;

      /* add up the total signal */
      for (j=0; j<GRID_SEGS; j++) {
	sig = bsig->signal[j];
	sig1 = bsig1.signal[j];
	kk = bsig1.hi_time[j];
	for (k=bsig1.lo_time[j]; k<=kk; k++) {
	  sig[k] += sig1[k]*e;
	}
      }
    }
  }

  /* adjust signals for t-zero */
  tt = t0;
  t0 -= tt;
  if (tt > 0) {
    for (j=0; j<GRID_SEGS; j++) {
      sig2 = bsig->signal[j];
      for (k=TIME_STEPS-1; k>=tt; k--) {
	sig2[k] = sig2[k-tt];
      }
      for (k=0; k<tt; k++) {
	sig2[k] = 0;
      }
    }
  }
  if (t0 >= 0.0001) {
    for (j=0; j<GRID_SEGS; j++) {
      sig2 = bsig->signal[j];
      for (k=TIME_STEPS-1; k>0; k--) {
	sig2[k] += t0 * (sig2[k-1]-sig2[k]);
      }
    }
  }

  /* calculate chisq */
  chisq = 0.0;
  for (i=0; i<MEAS_SEGS; i++) {
    for (j=0; j<TIME_STEPS; j++) {
      dif = bsig->signal[i][j] - asig->signal[i][j];
      chisq += (double) (dif*dif);
    }
  }

  return chisq;
} /* eval_int_pos */
