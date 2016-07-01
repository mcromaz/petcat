/*** petcat: offline decomp, preprocessor from GRETINA online src */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <assert.h>
#include <signal.h>
#include <getopt.h>
#include "petcat.h"
#include "mode3io.h"
#include "diagcnt.h"
#include "wrap.h"
#include "lh.h"
extern int preProcessMario(Mario *mario, Event_Signal *event, preprocCnt *diagcnt);
int h[4096];

int gotsignal = 0;
void breakhandler(int dummy) {
  printf("Got break signal, Aborting cleanly at next read...\n");
  gotsignal = 1;
}

int main(int argc, char **argv) {
  int m, outsideX;
  int stat;
  unsigned short int *buf;
  Event_Signal e;
  int ie;
  int inclevtlen, len = 0, num, i;
  long long int lenLong;
  int numxtalevts = 0;
  FILE *fin, *fou, *fspn, *fcfg;
  struct decomp_state *a;
  struct crys_intpts *x;
  mode3stream *evtstr;
  mode3Cnt *cnt;
  preprocCnt *pcnt;
  postprocCnt *postCnt;
  int holenum, xtalnum, singlePoint;
  int read3cnt;
  int new;	/* buffer new?, true!  */
  char inputFile[80], basisFile[80], detMapFile[80],
    filterFile[80], trGainFile[80], xTalkParsFile[80], s[80];
  struct gebData gh;  /* global header output with mode 2 data */
  Mario *mario;
  int verboseFlag = 0;
  char ch;

  struct option opts[] = {{"verbose", no_argument, 0, 'v'},
                          { 0, 0, 0, 0}};

  while ((ch = getopt_long(argc, argv, "v", opts, 0)) != -1) {
    switch(ch) {
    case 'v': verboseFlag = 1;
              fprintf(stdout, "I'm verbose ..\n");
              break;
    default: fprintf(stderr, "usage: petcat [-v]\n");
             exit(1);
    }
  }

  logMsg(verboseFlag, "the number is %d", 42);

  outsideX = 0;

  gotsignal = 0;
  signal(SIGINT, breakhandler);

  for (i = 0; i < 4096; i++) { h[i] = 0; }

  cnt = Calloc(sizeof(mode3Cnt), 1);
  pcnt = Calloc(sizeof(preprocCnt), 1);
  postCnt = Calloc(sizeof(postprocCnt), 1);
  mario = Calloc(sizeof(Mario), 1);
  pcnt->ener_cc = lh_init("ener_cc", 16384);

  initializePostCnt(postCnt);
  logMsg(verboseFlag, "initialized postCnt\n");

  fcfg = fopen("petcat.cfg", "r");
  assert(fcfg != 0);
  fscanf(fcfg, "%s %d", s, &holenum);
  fscanf(fcfg, "%s %d", s, &xtalnum);
  fscanf(fcfg, "%s %s", s, inputFile);
  fscanf(fcfg, "%s %s", s, basisFile);
  fscanf(fcfg, "%s %s", s, detMapFile);
  fscanf(fcfg, "%s %s", s, filterFile);
  fscanf(fcfg, "%s %s", s, trGainFile);
  fscanf(fcfg, "%s %s", s, xTalkParsFile);
  fclose(fcfg);

  logMsg(verboseFlag, "holenum = %d, xtalnum = %d\n", holenum, xtalnum);
  logMsg(verboseFlag,"inputFile = %s, basisFile = %s\n", inputFile, basisFile);

  /* Check ending of filename to decide if it's compressed or not.
     Compressed files are read in using the pipe (popen) method.
     HLC - May 29, 2013 */
  size_t flen = strlen(inputFile);
  size_t end1 = strlen(".gz");
  size_t end2 = strlen(".gzip");
  int compressed = ((strncmp(inputFile + flen - end1, ".gz", end1) == 0) ||
		     (strncmp(inputFile + flen - end2, ".gzip", end2) == 0));
  if (!compressed) {
    fin = fopen(inputFile, "r");
  } else {
    char inputFile2[80];
    sprintf(inputFile2, "zcat %s", inputFile);
    printf("Input file is compressed, opening in pipe: %s", inputFile2);
    fin = popen(inputFile2, "r");
  }
  assert(fin != 0);
  fou = fopen("out.dat", "w");
  assert(fou != 0);

  buf = Calloc(OBUF_LEN, sizeof(unsigned short int));

  inclevtlen = 286; //get_evt_len(fin);
  fprintf(stdout, "inclevtlen = %d (in short ints) for %s \n", inclevtlen, inputFile);

  stat = startPreProcess(inclevtlen, detMapFile, filterFile, trGainFile, xTalkParsFile);
  if (stat < 0) {
    fprintf(stderr, "startPreProcess failed!\n");
    exit(1);
  }

  a = dl_decomp_init(basisFile, 1); // 1 to suppress diag info
  assert(a != 0);
  //numGridPts(3);
  //exit(1);

  gh.type = 1;
  gh.length = sizeof(struct crys_intpts);

  e.total_energy = 1;
  while (((holenum > 100 && (lenLong = readMarioFormat(fin, mario, cnt)) > 0 && !gotsignal)) ||
	 ((len = read3(fin, holenum, xtalnum, buf, cnt)) > 0 && !gotsignal && holenum < 100)) {
    //fprintf(stdout, "-- read3(), len = %d\n", len);
    new = 1;
    if (holenum < 100) {
      while ( (holenum < 100 && (stat = preProcess(buf, len, inclevtlen, &e, new, pcnt)) >= 0 && !gotsignal) ) {
	new = 0;
	x = dl_decomp(a, &e, postCnt);
	log2intpts(x);
	x->crystal_id = holenum*4 + xtalnum; /* HLC -- Set crystal_id properly
						so later rotations, etc.
						make sense. */
	//printf("x->pad = %d\n", x->pad);
	gh.timestamp = x->timestamp;
	num = fwrite(&gh, sizeof(struct gebData), 1, fou);
	assert(num == 1);

	for (m=0; m<x->num; m++) {
	  if (x->intpts[m].x > 40.) { printf("%0.2f\n", x->intpts[m].x); outsideX++;}
	}
	num = fwrite(x, sizeof(struct crys_intpts), 1, fou);
	assert(num == 1);
	ie = (int) e.total_energy;
	if (ie >= 0 && ie < 4096) {
	  h[ie]++;
	}
      }
      bzero(buf, OBUF_LEN * sizeof(unsigned short int));
    } else if (holenum > 100 && (stat = preProcessMario(mario, &e, pcnt)) >= 0 && !gotsignal)  {
      e.time = lenLong;
      new = 0;
      x = dl_decomp(a, &e, postCnt);
      log2intpts(x);
      x->crystal_id = holenum*4 + xtalnum; /* HLC -- Set crystal_id properly
					      so later rotations, etc.
					      make sense. */
      //printf("x->pad = %d\n", x->pad);
      gh.timestamp = x->timestamp;
      num = fwrite(&gh, sizeof(struct gebData), 1, fou);
      assert(num == 1);

      for (m=0; m<x->num; m++) {
	if (x->intpts[m].x > 40.) { printf("%0.2f\n", x->intpts[m].x); outsideX++;}
      }
      num = fwrite(x, sizeof(struct crys_intpts), 1, fou);
      assert(num == 1);
      ie = (int) e.total_energy;
      if (ie >= 0 && ie < 4096) {
	h[ie]++;
      }
    }
  }
  ppMode3Cnt(cnt);
  ppPreprocCnt(pcnt);
  // ppPostprocCnt(postCnt);
  lh2spe(pcnt->ener_cc);

  fspn = fopen("segener.spn", "w");
  fwrite(pcnt->seg, sizeof(int), 36 * 4096, fspn);
  fclose(fspn);

  printf("# outside X radius: %d\n", outsideX);
}
