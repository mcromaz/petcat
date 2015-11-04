/*** oneInt.c: routines in support of single interaction point search */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "petcat.h"

int neigh[36][4] = {
  {1, 5, 6, -1}, {2, 0, 7, -1}, {3, 1, 8, -1}, {4, 2, 9, -1}, {5, 3, 10, -1}, {0, 4, 11, -1},
  {7, 11, 12, 0}, {8, 6, 13, 1}, {9, 7, 14, 2}, {10, 8, 15, 3}, {11, 9, 16, 4}, {6, 10, 17, 5},
  {13, 17, 18, 6}, {14, 12, 19, 7}, {15, 13, 20, 8}, {16, 14, 21, 9}, {17, 15, 22, 10}, {12, 16, 23, 11},
  {19, 23, 24, 12}, {20, 18, 25, 13}, {21, 19, 26, 14}, {22, 20, 27, 15}, {23, 21, 28, 16}, {18, 22, 29, 17},
  {25, 29, 30, 18}, {26, 24, 31, 19}, {27, 25, 32, 20}, {28, 26, 33, 21}, {29, 27, 34, 22}, {24, 28, 35, 23},
  {31, 35, -1, 24}, {32, 30, -1, 25}, {33, 31, -1, 26}, {34, 32, -1, 27}, {35, 33, -1, 28}, {30, 34, -1, 29}
};

void numGridPts(int stride) {
  /* diagnostic routine to determine grid pts/seg on fine/course grids */

  int s, r, phi, z;
  int numGridPts;

  assert(stride >= 1 && stride <= 4);

  for (s = 0; s < SSEG; s++ ) {
    numGridPts = 0;
    for (r = 0; r < SRAD; r += stride) {
      for (phi = 0; phi < SPHI; phi += stride) {
        for (z = 0; z < SZZZ; z += stride) {
          if (grid_pos_lu[s][r][phi][z] >= 0)  {
            numGridPts++;
          }
        }
      }      
    }
    if (s % 6 == 0) { printf("\nseg %2d:%2d   ", s, s + 5); }
    printf("%5d", numGridPts);
  }
  printf("\n");
  return; 
}

Basis_Point *oneIntSearchHelper(const Event_Signal *asig, int netSeg, int *neighSegArr, float *chisq) {

  Basis_Point *b;
  int r, phi, z, idx;
  int minDiffPos;
  float diff, minDiffVal;
  int i, n, neighSeg;

  minDiffVal = 10000.;
  minDiffPos = -1;
  for (r = 0; r < SRAD; r++) {
    for (phi = 0; phi < SPHI; phi++) {
      for (z = 0; z < SZZZ; z++) {
        idx = grid_pos_lu[netSeg][r][phi][z];
        if (idx != -1) { 
          b = &basis[idx];
          /* calc a chisq for cc and each of the neighbors */
          diff = 0.;
          for (i = 0 ; i < 50; i++) {	/* net segment */
            diff += fabs(asig->signal[netSeg][i] - b->signal[netSeg][i]);
          }
          for (n = 0; n < 4; n++) {  /* and neighbours */
            neighSeg = neighSegArr[n];
            if (neighSeg != -1) {
              for (i = 0; i < 50; i++) {
                diff += fabs(asig->signal[neighSeg][i] - b->signal[neighSeg][i]);
              }
            }
          }
        }
        if (diff < minDiffVal) {
          minDiffPos = idx;
          minDiffVal = diff;
        }
      }
    }
  }
  if (minDiffPos == -1) { return 0;}
  *chisq = (double) minDiffVal;
  return &basis[minDiffPos];
}

int badTrPred(const Event_Signal *asig, int nseg, int *seg) {

  int i, j;
  float sum;

  for (i = 0; i < nseg; i++) {	/* scaleToOne prob, 0 net tr */
    sum = 0.;
    for (j = 42; j < 50; j++) {
      sum += asig->signal[seg[i]][j];
    }
    if (sum == 0.) { return 1; }
  }
  return 0;
}

void scaleToOne(const Event_Signal *asig, Event_Signal *xsig, int netSeg) {

  float sum, scaleFactor;
  int i, n, ns;

  memcpy(xsig, asig, sizeof(Event_Signal));

  sum = 0.;
  for (i = 42; i < 50; i++) {
    sum += asig->signal[netSeg][i];
  }
  scaleFactor = sum / 8.;
  assert(scaleFactor != 0);

  for (i = 0; i < 50; i++) {
    xsig->signal[netSeg][i] /= scaleFactor;
  }

  for (n = 0; n < 4; n++) {
    ns = neigh[netSeg][n];
    if (ns != -1) {
      for (i = 0; i < 50; i++) {
        xsig->signal[ns][i] /= scaleFactor; 
      }
    }
  }

}

void subOverlapImage(const Event_Signal *asig, Event_Signal *xsig, int overlap, int net, int othernet, Basis_Point *otherpt) {

  float e1, e2, efrac;
  int i;

  memcpy(xsig, asig, sizeof(Event_Signal));

  e1 = asig->seg_energy[net]; 
  e2 = asig->seg_energy[othernet]; 
  efrac = e2 / (e1 + e2);
  scaleToOne(asig, xsig, net);
  for (i = 0; i < 50; i++) {
    xsig->signal[overlap][i] -= (1. - efrac) * otherpt->signal[overlap][i];
  }
  return;
}

int oneIntSearchIter(const Event_Signal *asig, Event_Signal *bsig, int numNet, int *netSeg, simpleInteraction *ints, double *chisq) {

  Event_Signal xsig;
  Basis_Point *b1, *b2, *tmp;
  float e1, e2, diff1, diff2;
  int neigh1[4], neigh2[4];
  int seg1, seg2, k1, k2, overlap, i, j;

  assert(asig != NULL && bsig != NULL);		/* pre */
  assert(numNet == 1 || numNet == 2);

  b1 = NULL;
  b2 = NULL;

  if (numNet == 1) {
    for (i = 0; i < 4; i++) {
      neigh1[i] = neigh[netSeg[0]][i];
    }
    b1 = oneIntSearchHelper(asig, netSeg[0], neigh1, &diff1);
  }

  if (numNet == 2) {
    assert(netSeg[0] != netSeg[1]);
    e1 = asig->seg_energy[netSeg[0]]; 
    e2 = asig->seg_energy[netSeg[1]]; 
    if (e1 > e2) {
      seg1 = netSeg[1];	// seg1 is the lowE pt
      seg2 = netSeg[0];
    } else {
      seg1 = netSeg[0];
      seg2 = netSeg[1];
    }

    /* generate neighbour arrays */
    for (i = 0; i < 4; i++) {
      k1 = neigh[seg1][i]; 
      neigh1[i] = (k1 == seg2) ? -1 : k1;
      k2 = neigh[seg2][i]; 
      neigh2[i] = (k2 == seg1) ? -1 : k2;
    }
    /* determine if there is an overlapping seg - there can only be one  */
    overlap = -1;		/* overlapping image charge segment */
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        if (neigh1[i] == neigh2[j] && neigh1[i] != -1) { 
          overlap = neigh1[i];
          goto there;
        }
      }
    }
there:
    
    /* generate neighbours for 1, excluding net seg from 2, and do search on 1 */
/*
    printf("dbg: net1=%d, neigh=%d,%d,%d,%d\n", seg1, neigh1[0], neigh1[1], neigh1[2], neigh1[3]);
    printf("dbg: net2=%d, neigh=%d,%d,%d,%d\n", seg2, neigh2[0], neigh2[1], neigh2[2], neigh2[3]);
*/
    scaleToOne(asig, &xsig, seg1);
    b1 = oneIntSearchHelper(&xsig, seg1, neigh1, &diff1);

    if (overlap == -1) {		/* segments can be searched independantly */
      scaleToOne(asig, &xsig, seg2);
      b2 = oneIntSearchHelper(&xsig, seg2, neigh2, &diff2);
    } else {
      /* subtract out fraction of overlap seg */
      subOverlapImage(asig, &xsig, overlap, seg2, seg1, b1);
      b2 = oneIntSearchHelper(&xsig, seg2, neigh2, &diff2);
    }
  }  

  if (numNet == 2 && seg2 == netSeg[0]) { /* swap pts to maintain orig interaction point ordering */
    tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  if (b1 != NULL) {
    ints[0].x = b1->x; 
    ints[0].y = b1->y; 
    ints[0].z = b1->z; 
    ints[0].e = asig->seg_energy[seg1]; 
    ints[0].diff = diff1; 
  }
  if (b2 != NULL) {
    ints[1].x = b2->x; 
    ints[1].y = b2->y; 
    ints[1].z = b2->z; 
    ints[1].e = asig->seg_energy[seg2]; 
    ints[1].diff = diff2; 
  }
  *chisq = (numNet == 1) ? diff1 : (diff1 + diff2) / 2.;		/* hack! */

  return 0;
}

int oneIntSearch2(const Event_Signal *asig, Event_Signal *bsig, int numNet, int *netSeg, simpleInteraction *ints, double *chisq) {

  Event_Signal xsig;
  Basis_Point *b[4];	/* up to four basis points */	
  int neighSegArr[4];
  float diff[4], adiff;
  int i, j, k, r, seg;

  assert(asig != NULL && bsig != NULL);		/* pre */
  assert(numNet >= 1 && numNet <= 4);

  if (numNet == 1) {
    for (i = 0; i < 4; i++) {
      neighSegArr[i] = neigh[netSeg[0]][i];
    }
    b[0] = oneIntSearchHelper(asig, netSeg[0], neighSegArr, &diff[0]);
    if (b[0] == 0) {return 1;}
  }
  else {
    for (i = 0; i < numNet; i++) {
      //assert(numNet == 2);
      /* generate neighbor array for for netSeg[i] */
      seg = netSeg[i];
      for (j = 0; j < 4; j++ ) {
        r = neigh[seg][j]; 
        /* if r is in the net list, set to -1 */
        for (k = 0; k < numNet; k++) {
          if (netSeg[k] == r) {
            r = -1;
            break;
          }
        }
        neighSegArr[j] = r;
      }
      /* norm asig and search */
      scaleToOne(asig, &xsig, netSeg[i]);
      b[i] = oneIntSearchHelper(&xsig, seg, neighSegArr, &diff[i]);
      if (b[i] == 0) { return 1;}   /* search failed */
    }
  }

  /* now, numNet basis points ident */
  for (i = 0; i < numNet; i++) {
    ints[i].x = b[i]->x; 
    ints[i].y = b[i]->y; 
    ints[i].z = b[i]->z; 
    ints[i].e = asig->seg_energy[netSeg[i]]; 
    ints[i].diff = diff[i]; 
  }

  adiff = 0.;
  for (i = 0; i < numNet; i++) {
    adiff += diff[i];
  }
  *chisq = adiff / numNet;	/* hack! */
  
  return 0; 
}

int oneIntSearch1(const Event_Signal *asig, Event_Signal *bsig, int numNet, int *netSeg, simpleInteraction *ints, double *chisq) {

  Event_Signal xsig;
  static int cnt = 0;
  Basis_Point *b1, *b2;
  float diff1, diff2;
  int neighSegArr[4];
  int i, k;

  assert(asig != NULL && bsig != NULL);		/* pre */
  assert(numNet == 1 || numNet == 2);

  b1 = NULL;
  b2 = NULL;

  if (numNet == 1) {
    for (i = 0; i < 4; i++) {
      neighSegArr[i] = neigh[netSeg[0]][i];
    }
    b1 = oneIntSearchHelper(asig, netSeg[0], neighSegArr, &diff1);
  }
  if (numNet == 2) {
    /* generate neighbour array for 1, excluding net seg from 2 and do search */
/*
    asigDiag(asig, netSeg[0]);
    scaleToOne(asig, &xsig, netSeg[0]);
    asigDiag(&xsig, netSeg[0]);
    if (cnt++ > 10 ) { exit(1); }
*/
    for (i = 0; i < 4; i++) {
      k = neigh[netSeg[0]][i]; 
      neighSegArr[i] = (k == netSeg[1]) ? -1 : k;
    }
    //printf("-dbg: netSeg[0]=%d, neigh=%d,%d,%d,%d\n", netSeg[0], neighSegArr[0], neighSegArr[1], neighSegArr[2], neighSegArr[3]);
    scaleToOne(asig, &xsig, netSeg[0]);
    b1 = oneIntSearchHelper(&xsig, netSeg[0], neighSegArr, &diff1);
    /* generate neighbour array for 2, excluding net seg from 1 and do search */
    for (i = 0; i < 4; i++) {
      k = neigh[netSeg[1]][i]; 
      neighSegArr[i] = (k == netSeg[0]) ? -1 : k;
    }
    //printf("-dbg: netSeg[1]=%d, neigh=%d,%d,%d,%d\n", netSeg[1], neighSegArr[0], neighSegArr[1], neighSegArr[2], neighSegArr[3]);
    scaleToOne(asig, &xsig, netSeg[1]);
    b2 = oneIntSearchHelper(&xsig, netSeg[1], neighSegArr, &diff2);
  }

  if (b1 != NULL) {
    ints[0].x = b1->x; 
    ints[0].y = b1->y; 
    ints[0].z = b1->z; 
    ints[0].e = asig->seg_energy[netSeg[0]]; 
    ints[0].diff = diff1; 
  }
  if (b2 != NULL) {
    ints[1].x = b2->x; 
    ints[1].y = b2->y; 
    ints[1].z = b2->z; 
    ints[1].e = asig->seg_energy[netSeg[1]]; 
    ints[1].diff = diff2; 
  }
/*
  printf("ints[0]: x=%5.2f, y=%5.2f, z=%5.2f, e=%5.2f\n", b1->x, b1->y, b1->z, asig->seg_energy[netSeg[0]]);
  printf("ints[1]: x=%5.2f, y=%5.2f, z=%5.2f, e=%5.2f\n", b2->x, b2->y, b2->z, asig->seg_energy[netSeg[1]]);
  exit(1);
*/
  *chisq = (numNet == 1) ? diff1 : (diff1 + diff2) / 2.;		/* hack! */
  return 0;
}


FILE *fou;
int open = 0;

void asigDiag(const Event_Signal *asig, int seg) {

  int *buf;
  int i, n, neighSeg, num;

  buf= calloc(4096, sizeof(int));
  assert(buf != 0);

  if (open == 0) {
    fou = fopen("out.spn", "w");
    open = 1;
  }

  printf("asigDiag: seg=%d, seg_ener=%f, total_energy=%f\n", seg, asig->seg_energy[seg], asig->total_energy);
  for (i = 0; i < 50; i++) {
    buf[i] = (int) (10000. * asig->signal[seg][i]);
  }
  for (n = 0; n < 4; n++) {
    neighSeg = neigh[seg][n];
    if (neighSeg != -1) {
      for (i = 0; i < 50; i++) {
        buf[i + (n + 1) * 100] = (int) (10000. * asig->signal[neighSeg][i]);
      }
    }
  }
  num = fwrite(buf, sizeof(int), 4096, fou);
  assert(num == 4096);
  free(buf);
}

int normDiag(const Event_Signal *asig, int nseg, int *seg) {

  return 0;
}
