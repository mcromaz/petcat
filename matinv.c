#include <math.h>
/* ======================================================================= */
int matinv(double *array, int norder, int dim)
{
  double amax, save, d1;
  int i, j, k, ik[1000], jk[1000];

  for (k = 0; k < norder; ++k) {
    /* find largest element array(i,j) in rest of matrix */
    amax = 0.f;
    while (1) {
      for (i = k; i < norder; ++i) {
	for (j = k; j < norder; ++j) {
	  if (fabs(amax) - (d1 = array[i + j*dim], fabs(d1)) <= 0.f) {
	    amax = array[i + j*dim];
	    ik[k] = i;
	    jk[k] = j;
	  }
	}
      }
      if (amax == 0.f) return 0;
      /* interchange rows and columns to put amax in array(k,k) */
      i = ik[k];
      if (i < k) continue;
      if (i > k) {
	for (j = 0; j < norder; ++j) {
	  save = array[k + j*dim];
	  array[k + j*dim] = array[i + j*dim];
	  array[i + j*dim] = -save;
	}
      }
      j = jk[k];
      if (j >= k) break;
    }
    if (j > k) {
      for (i = 0; i < norder; ++i) {
	save = array[i + k*dim];
	array[i + k*dim] = array[i + j*dim];
	array[i + j*dim] = -save;
      }
    }
    /* accumulate elements of inverse matrix */
    for (i = 0; i < norder; ++i) {
      if (i != k) array[i + k*dim] = -array[i + k*dim] / amax;
    }
    for (i = 0; i < norder; ++i) {
      for (j = 0; j < norder; ++j) {
	if (i != k && j != k)
	  array[i + j*dim] += array[i + k*dim] * array[k + j*dim];
      }
    }
    for (j = 0; j < norder; ++j) {
      if (j != k) array[k + j*dim] /= amax;
    }
    array[k + k*dim] = 1.f / amax;
  }
  /* restore ordering of matrix */
  for (k = norder-1; k >= 0; --k) {
    j = ik[k];
    if (j > k) {
      for (i = 0; i < norder; ++i) {
	save = array[i + k*dim];
	array[i + k*dim] = -array[i + j*dim];
	array[i + j*dim] = save;
      }
    }
    i = jk[k];
    if (i > k) {
      for (j = 0; j < norder; ++j) {
	save = array[k + j*dim];
	array[k + j*dim] = -array[i + j*dim];
	array[i + j*dim] = save;
      }
    }
  }
  return 0;
} /* matinv */

/* ======================================================================= */
int matinv_(double *array, int *norder, int *dim)
{
  return matinv(array, *norder, *dim);
}
int matinv__(double *array, int *norder, int *dim)
{
  return matinv(array, *norder, *dim);
}
