/* program to decompose event signals for GRETINA / GRETA.

   Author:  D.C. Radford    Aug 2004
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "petcat.h"

/* identifies near neighbors for hit segments */
static int  sel[36][37] = {
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


int grid_init(void)
{
  /*
    initializes the adaptive grid search algorithm
    by precalculating sums of products of basis signals
    ags1[] is for one hit segment; uses sel[][] to restrict calculated sums

    returns 1 on error
  */

  /* the six different segment rows have different volumes
     this ags_step gives the size of the course grid for each row
     2x2x2 grid for the front 3 rows, 3x3x3 grid for the back 3 rows */
  int    ags_step[6] = {2,2,3,3,3,3};

  double sum;
  float  *sig1, *sig2;
  int    i, j, k, point1, point2, pos, pos1, pos2, pair;
  int    s1, s2, ss, seg;

  ss = 0;
  for (seg=0; seg<SSEG; seg++) { /* loop over segment */
    ags1[seg].npts = 0;
      
    /* create the course grid for each of the segments */
    /* start from 1, not 0, because we're taking only some of the basis points */
    for (i=1; i<SRAD; i+=ags_step[seg/6]) {
      for (j=1; j<SPHI; j+=ags_step[seg/6]) {
	for (k=1; k<SZZZ; k+=ags_step[seg/6]) {
	  if ((pos = grid_pos_lu[seg][i][j][k]) >= 0) {
	    ags1[seg].grid_pos[ags1[seg].npts++] = pos;
	    if (ags1[seg].npts > MAX_AGS) {
	      printf("ERROR -- MAX_AGS = %d is too small for segment seg no. %d!\n"
		     "      -- i j k step = %d %d %d %d\n"
		     "      -- ags[seg].npts = %d\n"
		     "      -- Aborting.\n",
		     MAX_AGS, seg, i, j, k, ags_step[seg/6], ags1[seg].npts);
	      return 1;
	    }
	  }
	}
      }
    }

    /* malloc the space for ags[seg].da, ags[seg].db */
    s1 = sizeof(double) * ags1[seg].npts;
    if (!(ags1[seg].da = (double*) malloc(s1))) {
      printf("Ooops, ags[%d].da (%d) malloc failed...\n", seg, s1);
      exit(-1);
    }
    s2 = sizeof(float) * ags1[seg].npts * (ags1[seg].npts-1) / 2;
    if (!(ags1[seg].db = (float*) malloc(s2))) {
      printf("Ooops, ags[%d].db (%d) malloc failed...\n", seg, s2);
      exit(-1);
    }
    ss += s1 + s2;

    /* now calculate sums of basis_signal-squared and
       basis_signal_1*basis_signal_2 for single segments */
    pair = 0;
    for (point1=0; point1<ags1[seg].npts; point1++) {
      /* calculate sum of basis_signal-squared */
      pos1 = ags1[seg].grid_pos[point1];
      sum = 0.0;
      for (i=0; i<GRID_SEGS; i++) {
	if (!sel[seg][i]) continue;
	k = basis[pos1].hi_time[i];
	sig1 = basis[pos1].signal[i];
	for (j=basis[pos1].lo_time[i]; j<=k; j++) {
	  sum += (double) (sig1[j] * sig1[j]);
	}
      }
      ags1[seg].da[point1] = sum;
      
      for (point2=0; point2<point1; point2++) {
	/* calculate sum of basis_signal_1*basis_signal_2 */
	sum = 0.0;
	pos2 = ags1[seg].grid_pos[point2];
	for (i=0; i<GRID_SEGS; i++) {
	  if (!sel[seg][i]) continue;
	  k = basis[pos1].hi_time[i];
	  sig1 = basis[pos1].signal[i];
	  sig2 = basis[pos2].signal[i];
	  for (j=basis[pos1].lo_time[i]; j<=k; j++) {
	    sum += (double) (sig1[j] * sig2[j]);
	  }
	}
	ags1[seg].db[pair++] = sum;
      }
    }
    if (!quiet) printf("AGS init: segment seg %i has %i points, %i pairs, %i kB.\n",
		       seg, ags1[seg].npts, pair, (s1+s2)>>10);
  }
  if (!quiet) printf("ags[]da,db total = %dMB\n", ss>>20);

  return 0;
} /* grid_init */

/* ================================================================ */

double coarse_grid_1(const Event_Signal *asig,  /* observed signals */
		     int seg, Interaction *ints,
		     double *chisq0, double min_e_fraction)
{
  /*
    performs the coarse grid search algorithm for single-segment events
    returns chi-squared
    ints[] array is set to contain the resulting best-fit double-interaction pair
  */

  double *da, dd[MAX_AGS], sum, a, e1, e2, e1a, e2a, e1b, e2b, c2, chisq;
  const float  *sig1, *sig2, *db;
  int    i, j, k, n, pos1, pos2, pos1b, pos2b, pair;


  /* calculate chisq with no interactions */
  *chisq0 = 0.0;
  for (i=0; i<MEAS_SEGS; i++) {
    if (!sel[seg][i]) continue;
    sig1 = asig->signal[i];
    for (j=0; j<TIME_STEPS; j++) {
      *chisq0 += (double) (sig1[j] * sig1[j]);
    }
  }
  /* if (!quiet) printf("chisq0: %7.4f\n", *chisq0); */

  /* precalculate sums of (observed signal)*(basis signal) for AGS */
  for (n=0; n<ags1[seg].npts; n++) {
    pos1 = ags1[seg].grid_pos[n];
    sum = 0.0;
    for (i=0; i<GRID_SEGS; i++) {
      if (!sel[seg][i]) continue;
      k = basis[pos1].hi_time[i];
      sig1 = basis[pos1].signal[i];
      sig2 = asig->signal[i];
      for (j=basis[pos1].lo_time[i]; j<=k; j++) {
	sum += (double) (sig1[j] * sig2[j]);
      }
    }
    dd[n] = sum;
  }

  /* now do the grid search over the course grid for the hit segment */
  chisq = 100.0;

  pos1  = pos2  = -1;
  pos1b = pos2b = -1;
  e1a   = e2a   = 0.0;
  e1b   = e2b   = 0.0;
  pair  = 0;
  da    = ags1[seg].da;
  db    = ags1[seg].db;
  for (n=1; n<ags1[seg].npts; n++) {
    for (k=0; k<n; k++) {

      /*
	do 2-parameter linear least-squares for e1 and e2
	------------------------------
	e1 = fractional energy of 1st interaction
	e2 = fractional energy of 2nd interaction
	c2 = chi-squared
	------------------------------
	if we were using matinv, this would be:
	beta[1] = dd[n];
	beta[2] = dd[k];
	alpha[1][1] = ags1[seg].da[n];
	alpha[2][1] = ags1[seg].db[pair];
	alpha[1][2] = ags1[seg].db[pair];
	alpha[2][2] = ags1[seg].da[k];
	matinv(alpha, 2, 2);
	if (alpha[1][1]==0.0 || alpha[2][2]==0.0) ERROR;
	e1 = beta[1]*alpha[1][1] + beta[2]*alpha[1][2];
	e2 = beta[1]*alpha[2][1] + beta[2]*alpha[2][2];
      */
      a  = (da[n]*da[k] - db[pair]*db[pair]);
      e2 = (da[n]*dd[k] - db[pair]*dd[n])/a;
      e1 = (dd[n] - db[pair]*e2)/da[n];

      c2 = *chisq0 - dd[n]*e1 - dd[k]*e2;
      if (c2<chisq) {
	/* require that neither interaction is < 10% of total */
	if (e2<min_e_fraction) {
	  e2 = min_e_fraction;
	  e1 = (dd[n] - db[pair]*e2)/da[n];
	  c2 = *chisq0 - dd[n]*e1 - (2.0*dd[k] - db[pair]*e1 - da[k]*e2)*e2;
	} else if (e1<min_e_fraction) {
	  e1 = min_e_fraction;
	  e2 = (dd[k] - db[pair]*e1)/da[k];
	  c2 = *chisq0 - dd[k]*e2 - (2.0*dd[n] - db[pair]*e2 - da[n]*e1)*e1;
	}
	if (c2<chisq) {
	  /* this grid point pair has lower chi-squared than previous best
	     so save the parameters */
	  chisq = c2;
	  pos1b = pos1;
	  pos2b = pos2;
	  e1b = e1a;
	  e2b = e2a;
	  pos1 = ags1[seg].grid_pos[n];
	  pos2 = ags1[seg].grid_pos[k];
	  e1a  = e1;
	  e2a  = e2;
	  /* printf("grid: chisq, e1, e2, n, k, pair: %f %f %f %d %d %d\n",
	     chisq, e1, e2, n, k, pair); */
	}
      }
      pair++;
    }
  }

  ints[0].pos = pos1;
  ints[1].pos = pos2;
  ints[2].pos = pos1b;
  ints[3].pos = pos2b;
  ints[0].seg = ints[1].seg = ints[2].seg = ints[3].seg = seg;
  if (pos1<0 || pos2<0) {  /* (this should never happen) */
    printf("HRMM. No good grid positions found!\n");
    return chisq;
  }
  ints[0].r = basis[pos1].ir;
  ints[0].p = basis[pos1].ip;
  ints[0].z = basis[pos1].iz;
  ints[0].e = e1a;
  ints[1].r = basis[pos2].ir;
  ints[1].p = basis[pos2].ip;
  ints[1].z = basis[pos2].iz;
  ints[1].e = e2a;

  if (pos1b<0 || pos2b<0) return chisq;
  ints[2].r = basis[pos1b].ir;
  ints[2].p = basis[pos1b].ip;
  ints[2].z = basis[pos1b].iz;
  ints[2].e = e1b;
  ints[3].r = basis[pos2b].ir;
  ints[3].p = basis[pos2b].ip;
  ints[3].z = basis[pos2b].iz;
  ints[3].e = e2b;

  return chisq;
} /* coarse_grid_1 */

/* ================================================================ */

double refine_grid_1(const Event_Signal *asig,  /* observed signals */
		     double chisq, double chisq0, double min_e_fraction,
		     Interaction *ints)
{
  /*
    performs the adaptive part of grid search algorithms
    i.e. refines the grid search to a 1x1x1 grid
    this version does it for two interactions in one hit segment

    chisq = previous best chisq
            at present value of course rid search local_pars[]
    chisq0 = sum(obs. signal squared) = chisq with no interactions

    returns new refined value of chi-squared
    ints[] array contains the starting coarse-grid double-interaction pair
           and on exit is set to contain the new best-fit fine-grid pair
  */

  double sum1, sum2, chisq2, a, e1, e2, c2;
  double da1, db1, dd1, da2[27], dd2[27];
  const float  *sig1, *sig2;
  int    i, j, k, pos, pos1, pos2, pos11, pos22[27];
  int    ix1, iy1, iz1, ix2, iy2, iz2, i1, i2, i3, iter;
  int    npts2, pt2, seg1, seg2;


  chisq2 = chisq;
  seg1 = ints[0].seg;
  seg2 = ints[1].seg;

  /* now do adaptive part of AGS; i.e. check neighboring sites */
  for (iter = 0; iter < 6; iter++ ) {

    /* there are 27 possible neighboring grid points
       for each of the two interactions */
    /* here we precalculate the sums around the 2nd interaction */
    ix2 = ints[1].r - 0.5;
    iy2 = ints[1].p - 0.5;
    iz2 = ints[1].z - 0.5;
    npts2 = 0;
    for (i1=ix2; (i1<ix2+3 && i1<=maxir[seg2]); i1++) {
      if (i1<0) continue;
      for (i2=iy2; (i2<iy2+3 && i2<=maxip[seg2]); i2++) {
	if (i2<0) continue;
	for (i3=iz2; (i3<iz2+3 && i3<=maxiz[seg2]); i3++) {
	  if (i3<0) continue;
	  pos = grid_pos_lu[seg2][i1][i2][i3];
	  if (pos < 0) continue;
	  /* calculate sum of (basis signal)**2
	     calculate sum of (observed signal)*(basis signal) */
	  sum1 = 0.0;
	  sum2 = 0.0;
	  for (i=0; i<GRID_SEGS; i++) {
	    if (!sel[seg1][i]) continue;
	    k = basis[pos].hi_time[i];
	    sig1 = basis[pos].signal[i];
	    sig2 = asig->signal[i];
	    for (j=basis[pos].lo_time[i]; j<=k; j++) {
	      sum1 += (double) (sig1[j] * sig1[j]);
	      sum2 += (double) (sig1[j] * sig2[j]);
	    }
	  }
	  da2[npts2] = sum1;
	  dd2[npts2] = sum2;
	  pos22[npts2++] = pos;
	}
      }
    }

    /* now we calculate the sums around the 1st interaction
       and look at all combinations with neighbors to the 2nd interaction */
    pos1 = -1;
    pos2 = -1;
    ix1 = ints[0].r - 0.5;
    iy1 = ints[0].p - 0.5;
    iz1 = ints[0].z - 0.5;
    for (i1=ix1; (i1<ix1+3 && i1<=maxir[seg1]); i1++) {
      if (i1<0) continue;
      for (i2=iy1; (i2<iy1+3 && i2<=maxip[seg1]); i2++) {
	if (i2<0) continue;
	for (i3=iz1; (i3<iz1+3 && i3<=maxiz[seg1]); i3++) {
	  if (i3<0) continue;
	  pos11 = grid_pos_lu[seg1][i1][i2][i3];
	  if (pos11 < 0) continue;
	  /* calculate sum of (basis signal)**2
	     calculate sum of (observed signal)*(basis signal) */
	  da1 = 0.0;
	  dd1 = 0.0;
	  for (i=0; i<GRID_SEGS; i++) {
	    if (!sel[seg1][i]) continue;
	    k = basis[pos11].hi_time[i];
	    sig1 = basis[pos11].signal[i];
	    sig2 = asig->signal[i];
	    for (j=basis[pos11].lo_time[i]; j<=k; j++) {
	      da1 += (double) (sig1[j] * sig1[j]);
	      dd1 += (double) (sig1[j] * sig2[j]);
	    }
	  }
	  for (pt2=0; pt2<npts2; pt2++) {
	    /* calculate sum of (basis signal 1)*(basis signal 2) */
	    db1 = 0.0;
	    for (i=0; i<GRID_SEGS; i++) {
	      if (!sel[seg1][i]) continue;
	      k = basis[pos11].hi_time[i];
	      sig1 = basis[pos11].signal[i];
	      sig2 = basis[pos22[pt2]].signal[i];
	      for (j=basis[pos11].lo_time[i]; j<=k; j++) {
		db1 += (double) (sig1[j] * sig2[j]);
	      }
	    }

	    /* do the 2-parameter linear least-squares for e1 and e2 */
	    a  = da1*da2[pt2] - db1*db1;
	    e2 = (da1*dd2[pt2] - db1*dd1)/a;
	    e1 = (dd1 - db1*e2)/da1;
	    c2 = chisq0 - dd1*e1 - dd2[pt2]*e2;
	    if (c2<chisq2-0.00000001) {
	      /* require that neither interaction is < 10% of total */
	      if (e2<min_e_fraction) {
		e2 = min_e_fraction;
		e1 = (dd1 - db1*e2)/da1;
		c2 = chisq0 - dd1*e1 -(2.0*dd2[pt2] - db1*e1 - da2[pt2]*e2)*e2;
	      } else if (e1<min_e_fraction) {
		e1 = min_e_fraction;
		e2 = (dd2[pt2] - db1*e1)/da2[pt2];
		c2 = chisq0 - dd2[pt2]*e2 - (2.0*dd1 - db1*e2 - da1*e1)*e1;
	      }
	      if (c2<chisq2-0.00000001) {
		/* this grid point pair has lower chi-squared than previous best
		   so save the parameters */
		chisq2 = c2;
		pos1 = pos11;
		pos2 = pos22[pt2];
		ints[0].e = e1;
		ints[1].e = e2;
		/* printf("grid: chisq, e1, e2, n, k, pt2: %f %f %f %d %d %d",
		   chisq2, e1, e2, n, k, pt2 */
	      }
	    }
	  }
	}
      }
    }

    /* if no improved pair was found then quit this part of the algorithm */
    if (pos1<0 || pos2<0) break;
    /* else save the new parameters */
    ints[0].pos = pos1;
    ints[1].pos = pos2;
    ints[0].r = basis[pos1].ir;
    ints[0].p = basis[pos1].ip;
    ints[0].z = basis[pos1].iz;
    ints[1].r = basis[pos2].ir;
    ints[1].p = basis[pos2].ip;
    ints[1].z = basis[pos2].iz;
    if (!quiet) {
      printf("** grid2: chisq, pars: %7.4f", chisq2);
      for (i=0; i<2; i++) printf(", %6.3f %6.3f %6.3f %6.3f",
				 ints[i].r, ints[i].p, ints[i].z, ints[i].e);
      printf("\n");
    }
  }

  return chisq2;
} /* refine_grid_1 */
