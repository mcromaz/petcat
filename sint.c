/*** sint: single interaction search */

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "petcat.h"
#include "vegcat.h"

short tmp[4096];  /* tmp buffer for tr manipulation */

m3eb *to_m3eb(Mario *m) {
  m3eb *a;
  int i, j;
  a = calloc(1, sizeof(m3eb));
  a->tr_len = 300;
  a->cc_ener = m->ccEnergy;
  for (i = 0; i < 36; i++) {a->seg_ener[i] = m->segEnergy[i]; }
  for (i = 0; i < 37; i++) {
    for (j = 0; j < a->tr_len; j++) {
      a->tr[i][j] = (float) m->wf[i][j];
      //if (i == 36) { printf("%d, %d\n", j, m->wf[i][j]);}
    }
  }
  return a;
}

void pp_tr(m3eb *x, char *filename) {
  FILE *fou;
  int i, j;
  fou = fopen(filename, "w");
  fprintf(fou, "seg,ti,v\n");
  for (i = 0; i <37; i++) {
    for (j = 0; j < 300; j++) {
      fprintf(fou, "%d,%d,%f\n", i + 1, j + 1, x->tr[i][j]);
    }
  }
  fclose(fou);
  return;
}

int ident_net(m3eb *a, int *netlist, int pickoff) {

  float threshold = 100., le_avg, te_avg;
  int num_net = 0, i, j;
  for (i = 0; i < 36; i++) {
    for (j = 0, le_avg = 0.; j < 30; j++) { le_avg += a->tr[i][j]; }
    le_avg /= 30.;
    for (j = pickoff, te_avg = 0.; j < pickoff + 30; j++) { te_avg += a->tr[i][j]; }
    te_avg /= 30.;
    if (te_avg - le_avg > threshold) { netlist[num_net++] = i; }
  }
  return num_net;
}

float smooth[1024];
float diff[1024];
int t_fcfd(float *tr, int len, float cfd_arm_thresh, int debug) {
  float frac = 0.3;
  int offset = 20;
  int hf = 2; // half filter length
  int start, crossing;
  float sum, ravg = 0.;
  int i, j;
  FILE *fou;

  for (i = 0; i < 30; i++) {ravg += tr[i];}
  ravg /= 30.;
  ravg = 0.;//hack

  memset(smooth, 0, 1024 * sizeof(float)); // extra care for diag
  memset(diff, 0, 1024 * sizeof(float));

  for (i = hf; i < len - hf; i++) {
    sum = 0.;
    for (j = -hf; j <= hf; j++) { sum += tr[i + j]; }
    smooth[i] = sum / (2 * hf + 1) - ravg;
  }
  for (i = hf; i < len - hf - offset; i++) {
    diff[i] = smooth[i] - frac * smooth[i + offset];
  }

  start = -1, crossing = -1;
  for (i = hf; i < len - hf - offset; i++) {
    if (diff[i] < cfd_arm_thresh) { start = i; break; }
  }
  if (start != -1) {
    for (i = start; i < len - hf - offset; i++) {
      if (diff[i] < 0.) { crossing = i - 1; break; }
    }
  }

  if (debug == 1) {
    fou = fopen("t_fcfd.csv", "w");
    fprintf(fou, "ti,tr,smooth,diff\n");
    for (i = 0; i < len; i++) {
      fprintf(fou, "%d,%f,%f,%f\n", i + 1, tr[i], smooth[i], diff[i]);
    }
    fclose(fou);
  }
  return crossing;
}

int t_fled(float *tr, int len, float frac) {
  int hf = 2; // half filter length
  float sum, ravg, max, thresh;
  int i, j, max_i, t;
  memset(smooth, 0, 1024 * sizeof(float)); // extra care for diag
  for (i = 0; i < 30; i++) {ravg += tr[i];}
  ravg /= 30.;
  for (i = hf; i < len - hf; i++) {
    sum = 0.;
    for (j = -hf; j <= hf; j++) { sum += tr[i + j]; }
    smooth[i] = sum / (2 * hf + 1) - ravg;
  }
  max_i = -1, max = 0.;
  for (i = hf; i < len - hf; i++) {
    if (tr[i] > max) {
      max = tr[i], max_i = i;
    }
  }
  if (max_i == -1) { return -1; }
  thresh = frac * max;
  for (i = hf; i < len - hf; i++) {
    if (tr[i] > thresh) {t = i - 1; break;}
  }
  return t;
}

float tr_f[50];
void basis_cfd(int seg) {
  FILE *fou;
  Basis_Point *bPt;
  int cnt = 0, r, p, z, offset, t_seg, t_cc, t1, t2, i;
  fou = fopen("basis_cfd.csv", "w");
  fprintf(fou, "ir,ip,iz,x,y,z,r,t_seg,t_cc,t1,t2\n");
  for (r = 0; r < SRAD; r++ ) {
    for (p = 0; p < SPHI; p++) {
      for (z = 0; z < SZZZ; z++) {
        offset = grid_pos_lu[seg][r][p][z];
        if (offset >= 0) {
          cnt++;
          bPt = basis + offset;
          for (i = 0; i < 50; i++) {tr_f[i] = - (float) bPt->signal[seg][i];}
          t_seg = t_fcfd(tr_f, 50, 0.1, 0);
          for (i = 0, t1 = -1; i < 50; i++) {tr_f[i] = (float) bPt->signal[36][i];}
          t_cc = t_fcfd(tr_f, 50, 0.1, 0);
          for (i = 0, t1 = -1; i < 50; i++) {
            if (tr_f[i] < -0.2) {t1 = i - 1; break;}
          }
          for (i = 0, t2 = -1; i < 50; i++) {
            if (tr_f[i] < -0.2) {t2 = i - 1; break;}
          }
          fprintf(fou, "%d,%d,%d,%f,%f,%f,%f,%d,%d,%d,%d\n", bPt->ir, bPt->ip, bPt->iz,bPt->x, bPt->y, bPt->z,
          sqrt(bPt->x * bPt->x + bPt->y * bPt->y),t_seg,t_cc, t1, t2);
        }
      }
    }
  }
  printf("\nseg %d, %d basis points\n", seg, cnt);
  printf("exiting!!\n");
  fflush(stdout);
  exit(1);
}

float buf[1024];
double mvavg7max(float *x, int len) {

    double max;
    int i, j, max_i;

    max = -10000.;
    for (i = 3; i < len - 3; i++) {
      buf[i] = 0.;
      for (j = i - 3; j < i + 4; j++) {
        buf[i] += x[i + j];
      }
      if (buf[i] > max) {
        max = buf[i];
        max_i = i;
      }
    }
    return max / 7.;
}

float tr_gain[37];

m3eb *norm1(m3eb *x, int offset) {

  m3eb *y;
  double avg_bl, max, ampl_cc;
  int s, i, max_i;

  y = malloc(sizeof(m3eb));
  y->cc_ener = x->cc_ener;
  y->tr_len = x->tr_len;
  for (i = 0; i < 36; i++) { y->seg_ener[i] = x->seg_ener[i]; }
  // offset corr
  for (s = 0; s < 37; s++) {
    avg_bl = 0.;
    for (i = offset - 30; i < offset; i++) { avg_bl += x->tr[s][i]; }
    avg_bl /= 30.;
    for (i = 0; i < x->tr_len; i++) {
      y->tr[s][i] = x->tr[s][i] - avg_bl;
    }
  }

  ampl_cc = mvavg7max(&(y->tr[36][offset]), y->tr_len - offset);

  for (s = 0; s < 37; s++) {
    for (i = 0; i < y->tr_len; i++) {
      y->tr[s][i] = (y->tr[s][i] / ampl_cc); //  * (tr_gain[s] / tr_gain[36]);
    }
  }
  return y;
}

void sint_init(char *basisName, char *trGainName) {

  FILE *f_tr_gain;
  int i;

  if (read_basis(basisName) != 0) {
    fprintf(stderr, "could not open basis file %s\n", basisName);
    exit(1);
  }

  f_tr_gain = fopen(trGainName, "r");
  for (i = 0; i < 37; i++) {
    fscanf(f_tr_gain, "%f", tr_gain + i);
  }
  //basis_cfd(15);
  //exit(1);
  return;
}

int netlist[36];
extern int neigh[36][4];

float x2meas(m3eb *a, int bidx, int netseg, int offset, float *meas) {

  Basis_Point *b;
  int i, s, debug = 1;
  float diff, x2, x2_net, x2_neigh[4];

  b = basis + bidx;
  meas[0] = (float) bidx;

  x2_net = 0;
  for (i = 0; i < 50; i++) {
    diff = (a->tr[netseg][i + offset] - b->signal[netseg][i]);
    x2_net += diff * diff;
  }
  meas[1] = x2_net;

  for (s = 0; s < 4; s++) {  // iterate over neighbors
    x2_neigh[s] = 0.;
    for (i = 0; i < 50; i++) {
      diff = (a->tr[neigh[netseg][s]][i + offset] - b->signal[neigh[netseg][s]][i]);
      x2_neigh[s] += diff * diff;
    }
  }

  for (i = 0; i < 4; i++) { meas[i + 2] = x2_neigh[i]; }

  x2 = x2_net;
  for (i = 0; i < 4; i++) { x2 += x2_neigh[i]; }

  meas[6] = x2;
  return x2;
}

float meas_tbl[10000][7];

struct crys_intpts *sint(Mario *m, sdiag *sd) {

  struct crys_intpts *x;
  m3eb *a, *na;
  Basis_Point *b;
  FILE *fdbg;
  float min, meas;
  int t, num_net, seg, r, p, z, offset, cnt, debug = 1, bidx, i;

  min = 1000000.;
  a = to_m3eb(m);
  t = t_fled(&(a->tr[36][0]), a->tr_len, 0.2);
  num_net = ident_net(a, netlist, 200);
  na = norm1(a, 134);
  pp_tr(na, "tmpu.csv");

  if (num_net != 1) { return 0; }
  seg = netlist[0];
  cnt = 0;
  for (r = 0; r < SRAD; r++ ) {
    for (p = 0; p < SPHI; p++) {
      for (z = 0; z < SZZZ; z++) {
        bidx = grid_pos_lu[seg][r][p][z];
        if (bidx >= 0) {
          (void) x2meas(na, bidx, seg, 134, &meas_tbl[cnt++][0]);
        }
      }  // z
    }  // phi
  }  // r

  if (debug == 1) {
    fdbg = fopen("meas.csv", "w");
    fprintf(fdbg, "x,y,z,n,l,r,u,d,x2\n");

    for (i = 0; i < cnt; i++) {
      b = basis + (int) meas_tbl[i]
      [0]; // [0] is the basis offset
      fprintf(fdbg, "%5.2f,%5.2f,%5.2f,", b->x, b->y, b->z);
      fprintf(fdbg, "%f,%f,%f,%f,%f,%f\n", meas_tbl[i][1], meas_tbl[i][2],
      meas_tbl[i][3], meas_tbl[i][4], meas_tbl[i][5],meas_tbl[i][6]);
    }
    fclose(fdbg);
    exit(1);
  }


  printf("t = %d, num_net = %d, seg = %d, seg_e = %f, cc_e = %f\n", t, num_net, netlist[0],
a->seg_ener[netlist[0]], a->cc_ener);
  x = calloc(1, sizeof(struct crys_intpts));
  x->num = num_net;
  x->cfd = (float) t;
  x->intpts[0].seg = netlist[0];
  x->intpts[0].seg_ener = a->seg_ener[netlist[0]];
  x->intpts[0].e = 42. ;
  return x;
}
