#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include "petcat.h"
#include "vegcat.h"
#include "neigh.h"
#define MAX_EVT 100

int counter = 0, interval = 600;
char outdir[] = "basis_wevt/seg15";
char outfilename[100];

struct evtList {
  char *filename;
  int seg;
  Mario *rawEvts;
  int numEvts;
} runList[2] = {{"./basis_noise/seg15.dat", 15, 0, 0},
                {"./basis_noise/seg15.dat", 15, 0, 0}};

struct evtList pList[] = {{"none", 15, 0, 1},
			  {"none", 15, 0, 1},
			  {"none", 15, 0, 1},
			  {"none", 15, 0, 1},
			  {"none", 15, 0, 1},
			  {"none", 15, 0, 1}};  // fake it

float scratch[37][300];

//for calculating the timing of CFD Jan12 RT
double cf = -0.5;//should be negative value
int cfdelay = 10;
double thoffset[2] = {0.0};
int threshold = -50;
int runn[2][MAX_EVT];
double x[2][MAX_EVT], y[2][MAX_EVT], z[2][MAX_EVT];

Mario *wsum(Mario *evt0, Mario *evt1, float weight) {

  Mario *wevt;
  float s0, s1;
  int i, j;

  wevt = calloc(1, sizeof(Mario));

  for (i = 0; i < 37; i++) {
    for (j = 0; j < 300; j++) {
      s0 = (float)evt0->wf[i][j];
      s1 = (float) evt1->wf[i][j];
      wevt->wf[i][j] = (short int) (weight * s0 + (1. - weight) * s1);
    }
  }

  for (i = 0; i < 36; i++) {
    wevt->segEnergy[i] = weight * evt0->segEnergy[i] + (1.0 - weight) * evt1->segEnergy[i];
  }
  wevt->ccEnergy = weight * evt0->ccEnergy + (1.0 - weight) * evt1->ccEnergy;
  return wevt;
}

int main(int argc, char **argv) {

  FILE *fin, *fou[6], *ftr, *fgeoin, *fgeoout;
  struct gebData ghdr, defaulthdr = {100, sizeof(Mario), 0ll};
  Mario *evt;
  char evtlabel[] = {'a', 'b', 's'}, seglabel[] = {'n', 'l', 'r', 'u', 'd', 'c'};
  struct option opts[] = {{"verbose", no_argument, 0, 'v'},
                          {"sum-only", no_argument, 0, 's'},
                          {"append", no_argument, 0, 'a'},
                          {"ratio", required_argument, 0, 'r'}, 
                          {"cfd", no_argument, 0, 'c'}, //Use CFD
			  { 0, 0, 0, 0}};
  int verboseFlag = 0, sumOnlyFlag = 0, appendFlag = 0;
  int i, j, k, num, seg, baselineFlag = 1, CFDFlag = 0;
  double ratio = 0.5, ratioarray[6]={.1,.2,.4,.6,.8,.9};   // default
  char ch;

  for (i = 0; i < 2; i++) {
    runList[i].rawEvts = malloc(MAX_EVT * sizeof(Mario));
  }

  while ((ch = getopt_long(argc, argv, "vsar:c", opts, 0)) != -1) {
    switch(ch) {
    case 'v': verboseFlag = 1;
              fprintf(stdout, "I'm verbose ..\n");
              break;
    case 's': sumOnlyFlag = 1;
              break;
    case 'a': appendFlag = 1;
              break;
    case 'r': ratio = atof(optarg);
              break;
    case 'c': CFDFlag = 1;
              printf("CFDFalg is now set to 1\n");
              break;
    default: fprintf(stderr, "usage: cevt [-vsar:c]\n");
             exit(1);
    }
  }
  for (i = 0; i < 2; i++) {
    //counter = 0; //Not necessary to choose same geometrical points.
    fin = fopen(runList[i].filename, "r");
    sprintf(outfilename, "inputs/geo_basis_noise%i.txt", runList[i].seg);
    fgeoin = fopen(outfilename, "r");
    fgets(outfilename,20,fgeoin);//Dummy

    if (fin == 0) { fprintf(stderr, "could not open file %s\n", runList[i].filename); exit(1);}
    while ((runList[i].numEvts < MAX_EVT) && (fread(&ghdr, sizeof(struct gebData), 1, fin) == 1)) {
      if (ghdr.type == 100) {
	num = fread(runList[i].rawEvts + runList[i].numEvts, sizeof(Mario), 1, fin);
	assert(num == 1);
	fgets(outfilename,100,fgeoin);
	sscanf(outfilename,"%d, %lf, %lf, %lf\n",&runn[i][runList[i].numEvts],&x[i][runList[i].numEvts],&y[i][runList[i].numEvts],&z[i][runList[i].numEvts]);
	//printf("%d, %f, %f, %f\n",runn[i][runList[i].numEvts],x[i][runList[i].numEvts],y[i][runList[i].numEvts],z[i][runList[i].numEvts]);
	if (counter % 60 == 0) {
	  runList[i].numEvts++;
	}
 	counter++;
      }
      else {
	fseek(fin, ghdr.length, SEEK_CUR);
      }
    }
    fclose(fin);
  }
  fclose(fgeoin);



  sprintf(outfilename, "%s_geo.txt", outdir);
  fgeoout = fopen(outfilename, "w");
  fprintf(fgeoout, "runn, x, y, z \n");

  char outname[20];
  for(i = 0; i < 6; i++){
    sprintf(outname, "%s_%i.dat", outdir,(int)(10 * ratioarray[i]));
    fou[i] = (appendFlag == 1) ? fopen(outname, "a") : fopen(outname, "w");
    if (fou[i] == 0) { fprintf(stderr, "could not open file %s\n",outname); exit(1);}
  }

  for (i = 0; i < runList[0].numEvts; i++){
    for (j = 0; j < runList[0].numEvts; j++){
      for(k=0; k<6; k++){
	pList[k].rawEvts = wsum(runList[0].rawEvts + i, runList[1].rawEvts + j, ratioarray[k]);
	assert( 1 == fwrite(&defaulthdr, sizeof(struct gebData), 1, fou[k]));
	assert( 1 == fwrite(pList[k].rawEvts, sizeof(Mario), 1, fou[k]));
      }
      fprintf(fgeoout, "%07d, %f, %f, %f \n", 1000000+1000*i+j, x[0][i], y[0][i], z[0][i]);
      fprintf(fgeoout, "%07d, %f, %f, %f \n", 1000000+1000*i+j, x[1][j], y[1][j], z[1][j]);
    }
  }
  for(k=0; k<6; k++){
    fclose(fou[k]);
  }

  fclose(fgeoout);



  /*printf("list 1: %d evts, list 2: %d evts\n", runList[0].numEvts, runList[1].numEvts);

  if (appendFlag == 0) {
    CFDFlag == 1 ? sprintf(outname, "tr_cfd.csv") : sprintf(outname, "tr.csv");
    ftr = fopen(outname, "w");
    if (ftr == 0) { fprintf(stderr, "could not open file %s\n", outname); exit(1);}
    fprintf(ftr, "evt, seg, ch, val\n");
    for (i = 0; i < 3; i++) {
      if (sumOnlyFlag == 1 && i != 2) { continue; }
          evt = pList[i].rawEvts;  // take 1st evt
	  for (j = 0; j < 5; j++) {  // 'n', 'l', 'r', 'u', 'd'(, 'c')
	    seg = (j == 0) ? pList[i].seg : neigh[pList[i].seg][j - 1];
          for (k = 0; k < 300; k++) {
            fprintf(ftr, "%c,%c,%d,%d\n", evtlabel[i], seglabel[j], k + 1, evt->wf[seg][k]);
	  }
      }
      j = 5;
      seg = 36; // Center Contact
      for (k = 0; k < 300; k++) {
	fprintf(ftr, "%c,%c,%d,%d\n", evtlabel[i], seglabel[j], k + 1, evt->wf[seg][k]);
      }
    }
    fclose(ftr);
    }*/
  return 0;
}
