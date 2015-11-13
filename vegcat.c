#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include "petcat.h"
#include "vegcat.h"

int main(int argc, char **argv) {

  FILE *fin, *fou;
  char *inputFile;
  struct gebData ghdr1, ghdr2;
  Mario mario;
  Event_Signal e;
  struct decomp_state *a;
  struct crys_intpts *x;
  preprocCnt pcnt;
  postprocCnt postCnt;
  int holenum, xtalnum;
  char ch;
  int verboseFlag, stat, numhdr = 0, numevt = 0, num, i;
  struct option opts[] = {{"verbose", no_argument, 0, 'v'},
                          {"filename", required_argument, 0, 'f'},
                          { 0, 0, 0, 0}};

  while ((ch = getopt_long(argc, argv, "vf:", opts, 0)) != -1) {
    switch(ch) {
    case 'v': verboseFlag = 1;
              fprintf(stdout, "I'm verbose ..\n");
              break;
    case 'f': inputFile = optarg;
              break;
    default: fprintf(stderr, "usage: vegcat [-v]\n");
             exit(1);
    }
  }
  printf("vegcat\n");

  fin = fopen(inputFile, "r");
  if (fin == 0) { fprintf(stderr, "could not open file %s\n", inputFile); exit(1);}
  fou = fopen("out.csv", "w");

  stat = startPreProcess(100, cfg.detMapName, cfg.filterName, cfg.trGainName,
           cfg.xTalkParsName);
  if (stat < 0) { fprintf(stderr, "startPreProcess failed!\n"); exit(1); }

  a = dl_decomp_init(cfg.basisName, 1); // 1 to suppress diag info
  if (a == 0) { fprintf(stderr, "decomp init failed!\n"); exit(1); }

  ghdr2.type = 1;
  ghdr2.length = sizeof(struct crys_intpts);

  fprintf(fou, "evt, int, seg, x, y, z, e\n");

  while (fread(&ghdr1, sizeof(struct gebData), 1, fin) == 1) {
    numhdr++;
    if (ghdr1.type == 100) {
      numevt++;
      num = fread(&mario, sizeof(Mario), 1, fin);
      assert(num == 1);
      stat = preProcessMario(&mario, &e, &pcnt);
      e.time = ghdr1.timestamp;
      x = dl_decomp(a, &e, &postCnt);
      x->crystal_id = cfg.holenum * 4 + cfg.xtalnum; /* HLC -- Set crystal_id properly
					      so later rotations, etc.
					      make sense. */
      for (i = 0; i < x->num; i++) {
        fprintf(fou, "%d, %d, %d, %5.2f, %5.2f, %5.2f, %7.2f\n", numevt, i + 1,
          x->intpts[i].seg, x->intpts[i].x, x->intpts[i].y, x->intpts[i].z, x->intpts[i].e);
      }
      /*
      ghdr2.timestamp = x->timestamp;
      num = fwrite(&ghdr2, sizeof(struct gebData), 1, fou);
      assert(num == 1);
      num = fwrite(x, sizeof(struct crys_intpts), 1, fou);
      assert(num == 1);
      */
    }
    else {
      fseek(fin, ghdr1.length, SEEK_CUR);
    }
  }
  fprintf(stdout, "numhdr = %d, numevt = %d\n", numhdr, numevt);
}
