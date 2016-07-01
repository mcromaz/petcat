/*** time alignment routines */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "petcat.h"

int time_sef(int *tr, int len) {		/* from DecompIF.c */
	/* return index at which short energy filter has max value */

	int s1, s2, s3, s4;
	int i, max, max_i;
	#ifdef SAMPLE25
	static int PU_RISE = 20;	/* rise time, flat top parameters for energy filter */
	static int PU_FLAT = 10;
	#else
	static int PU_RISE = 40;	/* rise time, flat top parameters for energy filter */
	static int PU_FLAT = 20;
	#endif
	int MIN_T, MAX_T;
	static int dd[1024];

  #ifdef SAMPLE25
	MIN_T = 20.0;
	MAX_T = len - 20.0;
  #else
	MIN_T = 40.0;
	MAX_T = len - 40.0;
	#endif

	/* sums for derivative / pileup check */
	s1 = s2 = 0.0;
  s3 = 0;

  /* use PU_* trapezoid to calculate derivative of signal */
  for (i = MIN_T; i < MIN_T + PU_RISE; i++) {
    s1 += tr[i+PU_RISE+PU_FLAT] - tr[i]; /* first value of derivative */
  }
  dd[MIN_T] = s1;
  s2 += s1;
  s3++;

  for (i = MIN_T; i < MAX_T - (2*PU_RISE+PU_FLAT); i++) {
    /* remaining values of derivative */
    s1 += (tr[i+2*PU_RISE+PU_FLAT] - tr[i+PU_RISE+PU_FLAT] -
     tr[i+PU_RISE] + tr[i]);
    dd[i+1] = s1;
    s2 += s1;
    s3++;
  }

  /* check for max */
  s2 /= s3;
	max = -1.0;
	max_i = -1;
  for (i = MIN_T; i <= MAX_T - (2*PU_RISE+PU_FLAT); i++) {
    s4 = abs(s2 - dd[i]);		/* variation of derivitive from avg */
    if (max < s4) {
	    max = s4;
			max_i = i;
	  }
  }
	return max_i;
}

float t_cfd(int *buf, int buf_len)
{
  /* implement a Ge-style CFD to get the time from the signal trace
     buf = signal trace (input)
     returns CFD time as a float, in time steps, or -1.0 on error */

  int deriv[MAX_TRACE_LEN];
  int cfd[MAX_TRACE_LEN];
  int i, imax = 0, max_deriv = 0;

  assert(buf != 0 && buf_len > (5 - 2 * CFD_INT_LEN + CFD_DELAY) && buf_len < MAX_TRACE_LEN);

  deriv[0] = 0;
  for (i = 0; i < CFD_INT_LEN; i++) {
    deriv[0] += buf[i+CFD_INT_LEN] - buf[i];
  }
  for (i = 1; i < buf_len - 5 - 2*CFD_INT_LEN; i++) {
    deriv[i] = deriv[i-1] +
       buf[i+2*CFD_INT_LEN] - 2 * buf[i+CFD_INT_LEN] + buf[i];
    if (max_deriv < deriv[i]) {
      /* imax = time step where current pulse has its maximum value */
      max_deriv = deriv[i];
      imax = i;
    }
  }
  for (i = 0; i < buf_len - 5 - 2*CFD_INT_LEN - CFD_DELAY; i++) {
    cfd[i] = deriv[i] - deriv[i+CFD_DELAY]/CFD_FRACTION;
  }
  for (i = imax + CFD_DELAY; i > 0; i--) {
    if (cfd[i] <= 0 && cfd[i+1] > 0) {
      /* interpolate zero crossing and return time in steps, as a float */
      return ((float) i) - ((float) cfd[i]) / ((float) (cfd[i+1] - cfd[i]));
      break;
    }
  }
  return -1.0;
} /* t_cfd */

int align_cfd_1(int *traces, int tr_len, float *delay1, float *delay0)
{
  /* align the traces in an observed event,
     using the times from t_cfd() as a t0 reference.
     This routine uses interpolation to align over fractional time steps.
     It is used in eb2esig to get events for gdecomp */

  int j, k;
  int *tr;
  int n = 0, it;
  float s0 = 0.0f;
  float t, tt, tmin = 1000.0;

  assert(delay0 != 0 && delay1 != 0);

  for (j = 0; j < NSIGS; j++) {
    tr = traces + j * tr_len;
    if (net(tr, tr_len)) {
      t = t_cfd(tr, tr_len);
      if (t < 0.0) return -1;
      if (j < NSIGS - 1) t -= delay1[j];  // DCR: we are 99% sure about the - sign
      if (tmin > t) tmin = t;
      if (j < NSIGS - 1) {  // a segment signal, not the CC
	s0 += delay0[j];  // to calculate a mean of the delay0 values for hit segments
	/* note that the average value for delay0 has already been subtracted
	   from the individual values in the calling program, namely eb2esig */
	n++;
      }
    }
  }
  //t += 4.0f;
  if (n > 0 ) {
    t += s0 / ((float) n);  // correct for mean delay0
  } else {
    fprintf(stderr, "Warning! No net in event\n");
  }
  if (t < 0.0) return 0;

  for (j=0; j<NSIGS; j++) {
    tt = t;
    if (j < NSIGS-1) tt += delay1[j];  // a segment signal, not the CC, so correct for delay1
    if (tt <= 0.0f) {
      fprintf(stderr, "Warning! tt < 0 (%f) for seg %d in align_cfd_1\n", tt, j);
      continue;
    }
    /* do the interpolating and shifting of the traces */
    it = tt;             // integer shift
    tt -= (float) it;    // fractional shift for interpolation
    //for (k=0; k<50; k++) {  // do the shift & interpolate
    for (k=0; k<TIME_STEPS; k++) {  // do the shift & interpolate
      traces[j * tr_len + k] = (1.0f - tt) * traces[j * tr_len + k + it] + tt * traces[j * tr_len + k + it + 1];
    }
  }

  return (int) (100.0f * t);  /* DCR modified to propagate time offset to gdecomp */
} /* align_cfd_1 */
