#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "petcat.h"
#include "vegcat.h"
//#include "neigh.h"

int main(int argc, char **argv) {
  char *outdir="Mar14";

  FILE *fin, *fou, *ftr, *fl, *fave;
  char *inputFile, *inputFileList, s[200];
  int cnt, numTok, numRuns;
  struct gebData ghdr;
  Mario *rawEvts, *evt;
  Event_Signal e;
  struct decomp_state *a;
  struct crys_intpts *x;
  preprocCnt pcnt;
  postprocCnt postCnt;
  int holenum, xtalnum;
  char ch;
  int verboseFlag = 0, fileListFlag = 0;
  //int stat, numhdr = 0, numevt = 0, maxevt = 100, num, seg, i, j, k, l;
  int stat, numhdr = 0, numevt = 0, maxevt = 200000, num, seg, i, j, k, l;
  char seglabel[] = {'n', 'l', 'r', 'u', 'd'};
  const int maxruns = 1;// used to be 200
  struct option opts[] = {{"verbose", no_argument, 0, 'v'},
                          {"filename", required_argument, 0, 'f'},
                          {"numevts", required_argument, 0, 'n'},
                          {"filelist", required_argument, 0, 'l'},
                          {"segid", required_argument, 0, 's'},
                          { 0, 0, 0, 0}};
  struct {
    char *filename;
    int run;
    int seg;
    Mario *rawEvts;
    int numEvts;
  } runList[maxruns];

  inputFileList = malloc(200 * sizeof(char));

  while ((ch = getopt_long(argc, argv, "vf:l:n:s:", opts, 0)) != -1) {
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
    case 's': seg = atoi(optarg);
              break;
    default: fprintf(stderr, "usage: vegcat [-vf:l:n:]\n");
             exit(1);
    }
  }
  //seg = atoi(argv[1]); //Segment num
  printf("vegcat id %i\n", seg);


  for (i = 0; i < maxruns; i++) {
    runList[i].rawEvts = malloc(maxevt * sizeof(Mario));
    runList[i].filename = malloc(80);
    runList[i].numEvts = 0;
  }

  if (fileListFlag == 1) {
    sprintf(s, "%s.txt", inputFileList);
    fl = fopen(s, "r");
    if (fl == 0) { fprintf(stderr, "could not open file %s\n", inputFileList); exit(1);}
    cnt = 0;
    //while (cnt < maxruns && fgets(s, 80, fl) != 0) {
    while ( fgets(s, 80, fl) != 0) {
      numTok = sscanf(s, "%s %d %d", runList[cnt].filename, &runList[cnt].run, &runList[cnt].seg);
      if (numTok != 3) { fprintf(stderr, "wrong fmt - line %d, %s\n", cnt + 1, inputFileList); exit(1);}
      //if (runList[cnt].seg == seg) break;
      if (runList[cnt].run == seg) break; //this alos compatible
      //cnt++;
    }
  } else {
    strncpy(inputFileList, inputFile, 80);// inputFileList is also used for filename of output
    strncpy(runList[0].filename, inputFile, 80);
    runList[0].run = 1, runList[0].seg = 15; // defaults
    cnt = 1;
  }
  //numRuns = cnt;
  
  numRuns = 1;
  //assert(( fou = fopen("out.csv", "w")) != 0);
  sprintf(s, "./%s/%s%i_out.csv", outdir, inputFileList, seg);
  assert((fou = fopen(s, "w")) != 0);

  stat = startPreProcess(100, cfg.detMapName, cfg.filterName, cfg.trGainName,
           cfg.xTalkParsName);
  if (stat < 0) { fprintf(stderr, "startPreProcess failed!\n"); exit(1); }

  //a = dl_decomp_init(cfg.basisName, 0); // 1 to suppress diag info
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

  sprintf(s, "./%s/%s%i_avepos.csv",outdir, inputFileList, seg);
  fave = fopen(s, "w");
  //fprintf(fave,"filename, evtnum, avex, avey, avez, aveE\n");
  fprintf(fave,"runnum, segnum, evtnum, avex, avey, avez, aveE\n");
  fprintf(fou, "run, evt, int, seg, x, y, z, efrac, e, chisq\n");
  for (i = 0; i < numRuns; i++) {
    rawEvts = runList[i].rawEvts;
    double ave[5] = {0.0};//RT
    for (j = 0; j < runList[i].numEvts; j++) {
      stat = preProcessMario(rawEvts + j, &e, &pcnt);
      x = dl_decomp(a, &e, &postCnt);
      x->crystal_id = cfg.holenum * 4 + cfg.xtalnum; /* HLC -- Set crystal_id properly
					      so later rotations, etc.
					      make sense. */
      for (k = 0; k < x->num; k++) {
      fprintf(fou, "%d, %d, %d, %d, %5.2f, %5.2f, %5.2f, %7.2f, %7.2f, %7.2f, %7.2f, %7.2f, %i\n", runList[i].run, j + 1, k + 1,
	      x->intpts[k].seg, x->intpts[k].x, x->intpts[k].y, x->intpts[k].z, x->intpts[k].e, x->intpts[k].seg_ener, x->norm_chisq, x->t0, x->cfd, x->trig_time); // seg_ener added
      ave[0] += x->intpts[k].x/((double) (x->num * runList[i].numEvts));
      ave[1] += x->intpts[k].y/((double) (x->num * runList[i].numEvts));
      ave[2] += x->intpts[k].z/((double) (x->num * runList[i].numEvts));
      ave[3] += x->intpts[k].e/((double) (x->num * runList[i].numEvts));
      ave[4] += x->intpts[k].seg_ener/((double) (x->num * runList[i].numEvts));
      }
    }
    fprintf(stdout, "%s, %d evts\n", runList[i].filename, runList[i].numEvts);
    fprintf(stdout, "avex:%f, avey:%f, avez:%f, avee1:%f avee2:%f\n\n", ave[0], ave[1], ave[2], ave[3], ave[4]);
    //fprintf(fave, "%s, %d evts, ", runList[i].filename, runList[i].numEvts);
    fprintf(fave, "%i, %i, %d, ", runList[i].run, runList[i].seg, runList[i].numEvts);
    fprintf(fave, "%f, %f, %f, %f\n", ave[0], ave[1], ave[2], ave[3]);//ave[4] to be added later
  }
  fclose(fave);
  return 0;
}
