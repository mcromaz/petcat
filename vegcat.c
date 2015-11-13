#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "petcat.h"
#include "vegcat.h"

int main(int argc, char **argv) {

  FILE *fin, *fou, *fl;
  char *inputFile, *inputFileList, s[80];
  int cnt, numTok, numRuns;
  struct gebData ghdr;
  Mario *rawEvts;
  Event_Signal e;
  struct decomp_state *a;
  struct crys_intpts *x;
  preprocCnt pcnt;
  postprocCnt postCnt;
  int holenum, xtalnum;
  char ch;
  int verboseFlag = 0, fileListFlag = 0;
  int stat, numhdr = 0, numevt = 0, maxevt = 100, num, i, j, k;
  struct option opts[] = {{"verbose", no_argument, 0, 'v'},
                          {"filename", required_argument, 0, 'f'},
                          {"numevts", required_argument, 0, 'n'},
                          {"filelist", required_argument, 0, 'l'},
                          { 0, 0, 0, 0}};
  struct {
    char *filename;
    int run;
    Mario *rawEvts;
    int numEvts;
  } runList[32];

  while ((ch = getopt_long(argc, argv, "vf:l:n:", opts, 0)) != -1) {
    switch(ch) {
    case 'v': verboseFlag = 1;
              fprintf(stdout, "I'm verbose ..\n");
              break;
    case 'f': inputFile = optarg;
              break;
    case 'l': inputFileList = optarg;
              fileListFlag = 1;
              break;
    case 'n': maxevt = atoi(optarg);
              break;
    default: fprintf(stderr, "usage: vegcat [-vf:l:n:]\n");
             exit(1);
    }
  }
  printf("vegcat\n");

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
      numTok = sscanf(s, "%s %d", runList[cnt].filename, &runList[cnt].run);
      if (numTok != 2) { fprintf(stderr, "wrong fmt - line %d, %s\n", cnt + 1, inputFileList); exit(1);}
      cnt++;
    }
  } else {
    strncpy(runList[0].filename, inputFile, 80);
    runList[0].run = 1; // default
    cnt = 1;
  }
  numRuns = cnt;

  assert(( fou = fopen("out.csv", "w")) != 0);

  stat = startPreProcess(100, cfg.detMapName, cfg.filterName, cfg.trGainName,
           cfg.xTalkParsName);
  if (stat < 0) { fprintf(stderr, "startPreProcess failed!\n"); exit(1); }

  a = dl_decomp_init(cfg.basisName, 1); // 1 to suppress diag info
  if (a == 0) { fprintf(stderr, "decomp init failed!\n"); exit(1); }

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

  fprintf(fou, "run, evt, int, seg, x, y, z, e\n");

  for (i = 0; i < numRuns; i++) {
    rawEvts = runList[i].rawEvts;
    for (j = 0; j < runList[i].numEvts; j++) {
      stat = preProcessMario(rawEvts + j, &e, &pcnt);
      x = dl_decomp(a, &e, &postCnt);
      x->crystal_id = cfg.holenum * 4 + cfg.xtalnum; /* HLC -- Set crystal_id properly
					      so later rotations, etc.
					      make sense. */
      for (k = 0; k < x->num; k++) {
      fprintf(fou, "%d, %d, %d, %d, %5.2f, %5.2f, %5.2f, %7.2f\n", runList[i].run, j + 1, k + 1,
        x->intpts[k].seg, x->intpts[k].x, x->intpts[k].y, x->intpts[k].z, x->intpts[k].e);
      }
    }
    fprintf(stdout, "%s, %d evts\n", runList[i].filename, runList[i].numEvts);
  }
  return 0;
}
