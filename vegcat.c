#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <math.h>
#include "petcat.h"
#include "vegcat.h"
#include "neigh.h"

int main(int argc, char **argv) {

  FILE *fin, *fou, *ftr, *fl;
  char *inputFile, *inputFileList, s[80];
  int cnt, cnt1, cnt2, cnt3, numTok, numRuns;
  struct gebData ghdr;
  Mario *rawEvts, *evt;
  Event_Signal e;
  struct decomp_state *a;
  struct crys_intpts *x;
  sdiag *sd;
  preprocCnt pcnt;
  postprocCnt postCnt;
  int holenum, xtalnum;
  Basis_Point *b;
  char ch;
  int verboseFlag = 0, fileListFlag = 0, sintFlag = 0, testFlag = 0;
  int stat, numhdr = 0, numevt = 0, maxevt = 100, num, seg, i, j, k, l, m, n;
  char seglabel[] = {'n', 'l', 'r', 'u', 'd'};
  struct option opts[] = {{"verbose", no_argument, 0, 'v'},
                          {"sint", no_argument, 0, 's'},
                          {"filename", required_argument, 0, 'f'},
                          {"numevts", required_argument, 0, 'n'},
                          {"filelist", required_argument, 0, 'l'},
                          {"debug", no_argument, 0, 'd'},
                          { 0, 0, 0, 0}};
  struct {
    char *filename;
    int run;
    int seg;
    Mario *rawEvts;
    int numEvts;
  } runList[32];

  struct config {
    int holenum;
    int xtalnum;
    char *basisName;
    char *detMapName;
    char *filterName;
    char *trGainName;
    char *xTalkParsName;
  } cfg = {109, 0, "../coinc/q4a8_basis_xt.dat", "../coinc/detmap_Q4pos4_CC09.txt", "../coinc/filter.txt",
        "../coinc/tr_gain_Q4pos4_CC09.txt", "../coinc/q4a8_xtalk_pars_in.txt"};

  while ((ch = getopt_long(argc, argv, "vsf:l:n:t", opts, 0)) != -1) {
    switch(ch) {
    case 'v': verboseFlag = 1;
              fprintf(stdout, "I'm verbose ..\n");
              break;
    case 's': sintFlag = 1;
              break;
    case 'f': inputFile = optarg;
              break;
    case 'l': inputFileList = optarg;
              fileListFlag = 1;
              break;
    case 'n': maxevt = atoi(optarg);
              break;
    case 't': testFlag = 1;
              break;
    default: fprintf(stderr, "usage: vegcat [-vf:l:n:]\n");
             exit(1);
    }
  }
  printf("vegcat\n");

  if (testFlag == 1) {
    fou = fopen("test1.csv", "w");
    assert(fou != 0);
    fprintf(fou, "x0,y0,z0,dx,dy,dz\n");
    fprintf(stdout, "basis test\n");
    a = dl_decomp_init(cfg.basisName, 1); // 1 to suppress diag info
    if (a == 0) { fprintf(stderr, "decomp init failed!\n"); exit(1); }
    //(void) read_basis(cfg.basisName);
    cnt = 0, cnt1 = 0, cnt2 = 0;
    for (i = 0; i < MAX_SRAD; i++) {
      for (j = 0; j < MAX_SPHI; j++) {
        for (k = 0; k < MAX_SZZZ; k++) {
          if (grid_pos_lu[3][i][j][k] > 0) {
            b = basis + grid_pos_lu[3][i][j][k];
            // convert b to event signal
            memset(&e, 0, sizeof(Event_Signal));
            for (m = 0; m < 37; m++) {
              for (n = 0; n < 50; n++) {
                e.signal[m][n] = b->signal[m][n];
              }
            }
            e.total_energy = 1000.;
            e.seg_energy[3] = 1000.;
            e.core_e[0] = e.core_e[1] = e.core_e[2] = e.core_e[3] = 1000.;
            x = dl_decomp(a, &e, &postCnt);
            if (x->num == 1) {
              fprintf(fou, "%f,%f,%f,%f,%f,%f\n", b->x, b->y, b->z,
                fabs(b->x - x->intpts[0].x), fabs(b->y - x->intpts[0].y),
                fabs(b->z - x->intpts[0].z));
              cnt1++;
            }
            if (x->num == 2) {
              cnt2++;
              fprintf(stdout, "%f,%f,%f\n", b->x, b->y, b->z);
            }
            if (x->num == 3) { cnt3++; }
            //printf("cnt = %d, x->num = %d\n", cnt, x->num);
            cnt++;
          }
        }
      }
    }
    fprintf(stdout, "cnt = %d, cnt1 = %d, cnt2 = %d, cnt3 = %d\n", cnt, cnt1, cnt2, cnt3);
    exit(1);
  }

  for (i = 0; i < 32; i++) {
    runList[i].rawEvts = malloc(maxevt * sizeof(Mario));
    runList[i].filename = malloc(80);
    runList[i].numEvts = 0;
  }

  if (fileListFlag == 1) {
    fl = fopen(inputFileList, "r");
    if (fl == 0) { fprintf(stderr, "could not open file %s\n", inputFileList); exit(1);}
    cnt = 0;
    while (cnt < 32 && fgets(s, 80, fl) != 0) {
      numTok = sscanf(s, "%s %d %d", runList[cnt].filename, &runList[cnt].run, &runList[cnt].seg);
      if (numTok != 3) { fprintf(stderr, "wrong fmt - line %d, %s\n", cnt + 1, inputFileList); exit(1);}
      cnt++;
    }
  } else {
    strncpy(runList[0].filename, inputFile, 80);
    runList[0].run = 1, runList[0].seg = 15; // defaults
    cnt = 1;
  }
  numRuns = cnt;

  assert(( fou = fopen("out.csv", "w")) != 0);

  if (sintFlag == 1) {
    sint_init(cfg.basisName, cfg.trGainName);
  }
  else {        /* std decomp init */
    stat = startPreProcess(100, cfg.detMapName, cfg.filterName, cfg.trGainName,
             cfg.xTalkParsName);
    if (stat < 0) { fprintf(stderr, "startPreProcess failed!\n"); exit(1); }
    a = dl_decomp_init(cfg.basisName, 1); // 1 to suppress diag info
    if (a == 0) { fprintf(stderr, "decomp init failed!\n"); exit(1); }
  }

  for (i = 0; i < numRuns; i++) {
    fin = fopen(runList[i].filename, "r");
    if (fin == 0) { fprintf(stderr, "could not open file %s .. skipping\n", inputFile); continue;}
    while ((runList[i].numEvts < maxevt) && (fread(&ghdr, sizeof(struct gebData), 1, fin) == 1)) {
      if (ghdr.type == 100) {
        num = fread(runList[i].rawEvts + runList[i].numEvts++, sizeof(Mario), 1, fin);
        assert(num == 1);
      }
      else {
        fseek(fin, ghdr.length, SEEK_CUR);
      }
    }
    fclose(fin);
  }

  ftr = fopen("tr.csv", "w");
  if (ftr == 0) { fprintf(stderr, "could not open file tr.csv\n"); exit(1);}
  fprintf(ftr, "run,evt,seg,ch,val\n");
  for (i = 0; i < numRuns; i++) {
    for (j = 0; j < runList[i].numEvts; j++) {
      evt = runList[i].rawEvts + j;  // j'th evt in run i
      for (k = 0; k < 5; k++) {
        seg = (k == 0) ? runList[i].seg : neigh[runList[i].seg][k - 1];
        for (l = 0; l < 300; l++) {
          fprintf(ftr, "%d,%d,%c,%d,%d\n", runList[i].run, j + 1, seglabel[k],
                    l + 1, evt->wf[seg][l]);
        }
      }
    }
  }
  fclose(ftr);

  sd = calloc(1, sizeof(sdiag));  /* single interaction diagnositcs */

  fprintf(fou, "run, evt, tled, int, seg, x, y, z, e\n");
  for (i = 0; i < numRuns; i++) {
    rawEvts = runList[i].rawEvts;
    for (j = 0; j < runList[i].numEvts; j++) {
      if (sintFlag == 0) {
        stat = preProcessMario(rawEvts + j, &e, &pcnt);
        x = dl_decomp(a, &e, &postCnt);
      } else {
        x = sint(rawEvts + j, sd);
        exit(1); // hack
      }
      x->crystal_id = cfg.holenum * 4 + cfg.xtalnum; /* HLC -- Set crystal_id properly
					      so later rotations, etc.
					      make sense. */
      for (k = 0; k < x->num; k++) {
      fprintf(fou, "%d, %d, %5.1f, %d, %d, %5.2f, %5.2f, %5.2f, %7.2f\n", runList[i].run, j + 1, x->cfd, k + 1,
        x->intpts[k].seg, x->intpts[k].x, x->intpts[k].y, x->intpts[k].z, x->intpts[k].e);
      }
    }
    fprintf(stdout, "%s, %d evts\n", runList[i].filename, runList[i].numEvts);
  }
  return 0;
}
