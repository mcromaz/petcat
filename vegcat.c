#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include "petcat.h"
#include "vegcat.h"

int main(int argc, char **argv) {

  FILE *fin, *fou;
  char *inputFile;
  struct gebData ghdr;
  Mario *rawevts;
  Event_Signal e;
  struct decomp_state *a;
  struct crys_intpts *x;
  preprocCnt pcnt;
  postprocCnt postCnt;
  int holenum, xtalnum;
  char ch;
  int verboseFlag, stat, numhdr = 0, numevt = 0, maxevt = 100, num, i, j;
  struct option opts[] = {{"verbose", no_argument, 0, 'v'},
                          {"filename", required_argument, 0, 'f'},
                          {"numevts", required_argument, 0, 'n'},
                          { 0, 0, 0, 0}};

  while ((ch = getopt_long(argc, argv, "vf:n:", opts, 0)) != -1) {
    switch(ch) {
    case 'v': verboseFlag = 1;
              fprintf(stdout, "I'm verbose ..\n");
              break;
    case 'f': inputFile = optarg;
              break;
    case 'n': maxevt = atoi(optarg);
              break;
    default: fprintf(stderr, "usage: vegcat [-v]\n");
             exit(1);
    }
  }
  printf("vegcat\n");
  printf("sizeof(Mario) = %d\n", sizeof(Mario));

  rawevts = malloc(maxevt * sizeof(Mario));

  fin = fopen(inputFile, "r");
  if (fin == 0) { fprintf(stderr, "could not open file %s\n", inputFile); exit(1);}
  fou = fopen("out.csv", "w");

  stat = startPreProcess(100, cfg.detMapName, cfg.filterName, cfg.trGainName,
           cfg.xTalkParsName);
  if (stat < 0) { fprintf(stderr, "startPreProcess failed!\n"); exit(1); }

  a = dl_decomp_init(cfg.basisName, 1); // 1 to suppress diag info
  if (a == 0) { fprintf(stderr, "decomp init failed!\n"); exit(1); }

  while ((numevt < maxevt) && (fread(&ghdr, sizeof(struct gebData), 1, fin) == 1)) {
    numhdr++;
    if (ghdr.type == 100) {
      num = fread(&rawevts[numevt], sizeof(Mario), 1, fin);
      assert(num == 1);
      numevt++;
    }
    else {
      fseek(fin, ghdr.length, SEEK_CUR);
    }
  }
  printf("here!!\n");

  fprintf(fou, "evt, int, seg, x, y, z, e\n");

  for (i = 0; i < numevt; i++) {
    stat = preProcessMario(rawevts + i, &e, &pcnt);
    e.time = ghdr.timestamp;
    x = dl_decomp(a, &e, &postCnt);
    x->crystal_id = cfg.holenum * 4 + cfg.xtalnum; /* HLC -- Set crystal_id properly
					      so later rotations, etc.
					      make sense. */
    for (j = 0; j < x->num; j++) {
      fprintf(fou, "%d, %d, %d, %5.2f, %5.2f, %5.2f, %7.2f\n", i + 1, j + 1,
        x->intpts[j].seg, x->intpts[j].x, x->intpts[j].y, x->intpts[j].z, x->intpts[j].e);
    }
  }
  fprintf(stdout, "numhdr = %d, numevt = %d\n", numhdr, numevt);
}
