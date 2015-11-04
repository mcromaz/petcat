/* program to decompose event signals for GRETINA / GRETA.

   Author:  D.C. Radford    Aug 2004
*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "petcat.h"

/* identifies near neighbors for hit segments */
static int sel[36][37] = {
  {1,1,0,0,0,1, 1,1,0,0,0,1, 1,1,0,0,0,1, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {1,1,1,0,0,0, 1,1,1,0,0,0, 1,1,1,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {0,1,1,1,0,0, 0,1,1,1,0,0, 0,1,1,1,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {0,0,1,1,1,0, 0,0,1,1,1,0, 0,0,1,1,1,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {0,0,0,1,1,1, 0,0,0,1,1,1, 0,0,0,1,1,1, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {1,0,0,0,1,1, 1,0,0,0,1,1, 1,0,0,0,1,1, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {1,1,0,0,0,1, 1,1,0,0,0,1, 1,1,0,0,0,1, 1,1,0,0,0,1, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {1,1,1,0,0,0, 1,1,1,0,0,0, 1,1,1,0,0,0, 1,1,1,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {0,1,1,1,0,0, 0,1,1,1,0,0, 0,1,1,1,0,0, 0,1,1,1,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {0,0,1,1,1,0, 0,0,1,1,1,0, 0,0,1,1,1,0, 0,0,1,1,1,0, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {0,0,0,1,1,1, 0,0,0,1,1,1, 0,0,0,1,1,1, 0,0,0,1,1,1, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {1,0,0,0,1,1, 1,0,0,0,1,1, 1,0,0,0,1,1, 1,0,0,0,1,1, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {1,1,0,0,0,1, 1,1,0,0,0,1, 1,1,0,0,0,1, 1,1,0,0,0,1, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {1,1,1,0,0,0, 1,1,1,0,0,0, 1,1,1,0,0,0, 1,1,1,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {0,1,1,1,0,0, 0,1,1,1,0,0, 0,1,1,1,0,0, 0,1,1,1,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {0,0,1,1,1,0, 0,0,1,1,1,0, 0,0,1,1,1,0, 0,0,1,1,1,0, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {0,0,0,1,1,1, 0,0,0,1,1,1, 0,0,0,1,1,1, 0,0,0,1,1,1, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {1,0,0,0,1,1, 1,0,0,0,1,1, 1,0,0,0,1,1, 1,0,0,0,1,1, 0,0,0,0,0,0, 0,0,0,0,0,0 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 1,1,0,0,0,1, 1,1,0,0,0,1, 1,1,0,0,0,1, 0,0,0,0,0,0 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 1,1,1,0,0,0, 1,1,1,0,0,0, 1,1,1,0,0,0, 0,0,0,0,0,0 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,1,1,1,0,0, 0,1,1,1,0,0, 0,1,1,1,0,0, 0,0,0,0,0,0 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,1,1,1,0, 0,0,1,1,1,0, 0,0,1,1,1,0, 0,0,0,0,0,0 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,1,1,1, 0,0,0,1,1,1, 0,0,0,1,1,1, 0,0,0,0,0,0 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 1,0,0,0,1,1, 1,0,0,0,1,1, 1,0,0,0,1,1, 0,0,0,0,0,0 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 1,1,0,0,0,1, 1,1,0,0,0,1, 1,1,0,0,0,1 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 1,1,1,0,0,0, 1,1,1,0,0,0, 1,1,1,0,0,0 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,1,1,1,0,0, 0,1,1,1,0,0, 0,1,1,1,0,0 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,1,1,1,0, 0,0,1,1,1,0, 0,0,1,1,1,0 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,1,1,1, 0,0,0,1,1,1, 0,0,0,1,1,1 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 1,0,0,0,1,1, 1,0,0,0,1,1, 1,0,0,0,1,1 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 1,1,0,0,0,1, 1,1,0,0,0,1 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 1,1,1,0,0,0, 1,1,1,0,0,0 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,1,1,1,0,0, 0,1,1,1,0,0 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,1,1,1,0, 0,0,1,1,1,0 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,1,1,1, 0,0,0,1,1,1 ,1},
  {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 1,0,0,0,1,1, 1,0,0,0,1,1 ,1}};


double fitter(const Event_Signal *asig,  /* observed signals */
	      Event_Signal *bsig, /* fitted signals */
	      Interaction *ints, double *t0, int nints, int final) {

  /****************************************************************************
     This subroutine is a modified version of 'CURFIT', in Bevington....
     designed originally for use with gf2
     D.C. Radford                Mar 1984
     Modified for tcovergeantc   May 2003  DCR
     Translated to C             Aug 2004  DCR
     Cylindrical basis           Sep 2006  DCR
     T-zero fitting              Nov 2006  DCR

     Returns chi-squared
  ****************************************************************************/

  float  efrac = 0.02;
  double beta[MAX_PARS], delta[MAX_PARS], pars[MAX_PARS], pars_save[MAX_PARS];
  double array[MAX_PARS][MAX_PARS], alpha[MAX_PARS][MAX_PARS];
  double sqa[MAX_PARS][MAX_PARS], errs[MAX_PARS];
  double convfact, bb[MAX_PARS], xtry[MAX_PARS], maxflamda, minflamda;
  double chisq, chisq1, flamda, f, f2, step, max;
  int    nits, maxits, fail, i, j, k, ii, jj, kk, iseg, npars;
  int    nip, nextp[MAX_PARS], conv, fixed[MAX_PARS], ssel[37];

  /* Initialise for fitting */
  for (j=0; j<GRID_SEGS; j++) {
    ssel[j] = 0; /* neighbouring segment map array */
  }
  npars = nints*4;   /* npars = number of parameters */
  for (i=0; i<nints; i++) {  /* nints = number of fitted interactions */
    ii = 4*i;
    pars[ii]   = ints[i].r;
    pars[ii+1] = ints[i].p;
    pars[ii+2] = ints[i].z;
    pars[ii+3] = ints[i].e;
    for (j=0; j<GRID_SEGS; j++) {
      ssel[j] |= sel[ints[i].seg][j];  /* neighbouring segments; have information for fitting */
    }
  }
  pars[npars] = *t0;
  if (final > 0) {
    convfact = 10.0;
    maxflamda = 50.0;
    maxits = 20;
  } else {
    convfact = 1.0;
    maxflamda = 0.5;
    maxits = 6;
  }
  flamda = 0.1;
  minflamda = 0.0;
  nits   = 0;
  conv   = 0;
  nip    = npars+1;
  chisq  = 0.0;
  chisq1 = 0.0;
  step = 0.999;
  if (!quiet) {
    printf("nits nip  flamda   chisq   chisq1 :"
	   "    r      p      z      r      p      z      e1     e2  seg seg   t0\n");
    printf("%4d %3d %8.4f %7.4f %7.4f :"
	   " %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %2d %2d %5.2f\n",
	   nits, nip, flamda, chisq, chisq1, pars[0], pars[1],
	   pars[2], pars[4], pars[5], pars[6], pars[3], pars[7],
	   ints[0].seg, ints[1].seg, pars[npars]);
    for (j=8; j<npars; j+=8) {
      if (j+4 >= npars) {
	printf("%42.3f %6.3f %6.3f %27.3f %9d\n",
	       pars[j+0], pars[j+1], pars[j+2], pars[j+3], ints[j/4].seg);
      } else {
	printf("%42.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %2d %2d\n",
	       pars[j+0], pars[j+1], pars[j+2], pars[j+4],
	       pars[j+5], pars[j+6], pars[j+3], pars[j+7],
	       ints[j/4].seg, ints[1+j/4].seg);
      }
    }
  }

  while (!conv && nits<maxits) {

    /* evaluate fit, alpha & beta matrices, & chisq */
    for (j=0; j<npars+1; j++) {
      pars_save[j] = pars[j];
      beta[j] = 0.0;
      for (k=0; k<=j; k++) {
	alpha[j][k] = 0.0;
      }
      fixed[j] = 0;
    }

    for (i=0; i<nints; i++) {
      ii   = 4*i;
      ints[i].r = pars[ii];
      ints[i].p = pars[ii+1];
      ints[i].z = pars[ii+2];
      ints[i].e = pars[ii+3];
      ints[i].dr = errs[ii];
      ints[i].dp = errs[ii+1];
      ints[i].dz = errs[ii+2];
      ints[i].de = errs[ii+3];
    }
    *t0 = pars[npars];
    if (eval(asig, bsig, ints, *t0, nints, &chisq1, beta, alpha, 1, ssel)) {
      if (!quiet) printf("*** ACK!!! eval failure! ***\n");
      return 1000.0;
    }
    for (j=0; j<npars+1; j++) {
      bb[j] = beta[j];
      for (k=0; k<=j; k++) {
	if (alpha[j][j]*alpha[k][k] == 0.0) {
	  printf("Cannot - diag. element %d or %d eq. to zero.\n", j, k);
	  return 1000.0;
	}
	sqa[j][k] = sqrt(alpha[j][j]*alpha[k][k]);
	sqa[k][j] = sqa[j][k];
      }
    }

    /* alpha[j][k] = SUM of derivs[j]*derivs[k] for all points */

    /* invert modified curvature matrix to find new parameters */
  INVERT_MATRIX:
    nip = 0;
    for (i=0; i<npars+1; i++) {
      delta[i] = 0.0;
      if (!fixed[i]) nextp[nip++] = i;
    }

    array[0][0] = 1.0 + flamda;
    for (j=1; j<nip; j++) {
      jj = nextp[j];
      for (k=0; k<j; k++) {
	kk = nextp[k];
	array[j][k] = alpha[jj][kk]/sqa[jj][kk];
	array[k][j] = array[j][k];
      }
      array[j][j] = 1.0 + flamda;
    }

    matinv(array[0], nip, MAX_PARS);
    for (j=0; j<nip; j++) {
      if (alpha[j][j]*alpha[j][j] == 0.0) {
	printf("Cannot - diag. element %d eq. to zero.\n", j);
	return 1000.0;
      }
      jj = nextp[j];
      for (k=0; k<nip; k++) {
	kk = nextp[k];
	delta[jj] = delta[jj] + bb[kk]*array[j][k]/sqa[jj][kk];
      }
    }

    /* calculate new par. values */
    for (j=0; j<npars+1; j++) {
      xtry[j] = pars[j] + delta[j]*step;
    }

    /*
      check constraints on r, phi, z, energy-fraction and t0.
      if constraints are violated, fixed the worst one to be violated
      at its limiting value, and ajust parameters and derivatives,
      then fix the constrained parameter and reinvert matrix.
    */

    k = 0;
    f = 1.0;
    for (i=0; i<nip; i++) {
      ii = nextp[i];
      j = ii/4;
      iseg = ints[j].seg;
      j = ii - 4*j;
      f2 = 1.0;
      if (j < 3) {  /* radius, phi or z */
	if (ii == npars) {   /* t0 */
	  max = TIME_STEPS/2;
	} else if (j == 0) { /* r */
	  max = maxir[iseg];
	} else if (j == 1) { /* phi */
	  max = maxip[iseg];
	} else {             /* z */
	  max = maxiz[iseg];
	}
	if (xtry[ii] < 0.0) {
	  f2 = pars[ii]/(pars[ii] - xtry[ii]);
	} else if (xtry[ii] > max) {
	  f2 = (pars[ii] - max)/(pars[ii] - xtry[ii]);
	}
	if (f2 < f) {
	  k = ii;
	  f = f2;
	}
      } else if (j == 3 && xtry[ii] < efrac) {  /* energy fraction */
	f2 = (pars[ii] - efrac)/(pars[ii] - xtry[ii]);
	if (f2 < f) {
	  k = ii;
	  f = f2;
	}
      }
    }

    /* if we have fixed any parameters, go back and reinvert matrix */
    if (f < 1.0) {
      fixed[k] = 1;
      if (f > 0.01) {
	for (j=0; j<npars+1; j++) {
	  pars[j] += delta[j] * step * f;
	  bb[j] = bb[j] * (1.0-f);
	}
      }
      if (f < 0.9) goto INVERT_MATRIX;
    } else {
      for (j=0; j<npars+1; j++) {
	pars[j] = xtry[j];
      }
    }

    /* if chisq increased, increase flamda and try again */
    for (i=0; i<nints; i++) {
      ii   = 4*i;
      ints[i].r = pars[ii];
      ints[i].p = pars[ii+1];
      ints[i].z = pars[ii+2];
      ints[i].e = pars[ii+3];
      ints[i].dr = errs[ii];
      ints[i].dp = errs[ii+1];
      ints[i].dz = errs[ii+2];
      ints[i].de = errs[ii+3];
    }
    *t0 = pars[npars];
    fail = eval(asig, bsig, ints, *t0, nints, &chisq, beta, alpha, 0, ssel);
    if (!quiet) {
      printf("%4d %3d %8.4f %7.4f %7.4f : "
	     "%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %2d %2d %5.2f\n",
	     nits, nip, flamda, chisq, chisq1, pars[0], pars[1],
	     pars[2], pars[4], pars[5], pars[6], pars[3], pars[7],
	     ints[0].seg, ints[1].seg, pars[npars]);
      for (j=8; j<npars; j+=8) {
	if (j+4 >= npars) {
	  printf("%42.3f %6.3f %6.3f %27.3f %9d\n",
		 pars[j+0], pars[j+1], pars[j+2], pars[j+3], ints[j/4].seg);
	} else {
	  printf("%42.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %2d %2d\n",
		 pars[j+0], pars[j+1], pars[j+2], pars[j+4],
		 pars[j+5], pars[j+6], pars[j+3], pars[j+7],
		 ints[j/4].seg, ints[1+j/4].seg);
	}
      }
    }

    if (fail || chisq > chisq1) {
      if (fail && !quiet) printf("... eval failed!\n");
      for (j=0; j<npars+1; j++) {
	bb[j] = beta[j];
	pars[j] = pars_save[j];
	fixed[j] = 0;
      }
      if (flamda < maxflamda) {
	/* flamda *= sqrt(10.0); */
	flamda *= 10.0;
	if (flamda < 0.1) flamda = 0.1;
	goto INVERT_MATRIX;
      }
      conv = 0;
      break;
    }

    /*  evaluate parameters and errors, and test for convergence */
    conv = 1;
    /* f = sqrt(1.0 + flamda) * asig->total_energy/2.0; */
    f = asig->total_energy/2.0;
    for (j=0; j<nip; j++) {
      jj = nextp[j];
      if (array[j][j]<0.0) array[j][j]=0.0;
      errs[jj] = sqrt(array[j][j]/alpha[jj][jj])/f;
      if (fabs(delta[jj]) >= (errs[jj]/convfact)) conv = 0;
    }

    if (flamda > minflamda) flamda /= 10.0;
    nits++;
    if (!quiet) printf("\n");
  }

  /* list data and exit */
  if (conv) {
    chisq = eval_int_pos(asig, bsig, ints, *t0, nints);
    if (!quiet) printf("         Converged after %d iterations,  Chisq = %f\n",
		       nits, chisq);
  } else {
    /* if (nits == 0) chisq = chisq1; */
    for (i=0; i<nints; i++) {
      ii   = 4*i;
      ints[i].r = pars[ii];
      ints[i].p = pars[ii+1];
      ints[i].z = pars[ii+2];
      ints[i].e = pars[ii+3];
      ints[i].dr = errs[ii];
      ints[i].dp = errs[ii+1];
      ints[i].dz = errs[ii+2];
      ints[i].de = errs[ii+3];
    }
    *t0 = pars[npars];
    chisq = eval_int_pos(asig, bsig, ints, *t0, nints);
    if (!quiet) {
      printf("Failed to converge after %d iterations,  Chisq = %f\n",
	     nits, chisq);
    }
  }
  return chisq;
}
