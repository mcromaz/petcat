#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include "petcat.h"
#include "vegcat.h"
#include "neigh.h"
#define MAX_EVT 100

struct evtList {
  char *filename;
  int seg;
  Mario *rawEvts;
  int numEvts;
} runList[2] = {{"../coinc-data/WFOUT-Run0111Segment15.dat", 15, 0, 0},
                {"../coinc-data/WFOUT-Run0113Segment15.dat", 15, 0, 0}};

struct evtList pList[] = {{"none", 15, 0, 1},
                        {"none", 15, 0, 1},
                        {"none", 15, 0, 1}};  // fake it

float scratch[37][300];
int verboseFlag = 0, sumOnlyFlag = 0, appendFlag = 0, fileListFlag = 0;
int i, j, k, num, seg, baselineFlag = 1, CFDFlag = 0;

//for calculating the timing of CFD Jan12 RT
double cf = -0.5;//should be negative value
int cfdelay = 10;
double thoffset[2] = {0.0};
int threshold = -50;
FILE* frunlistout;


float cft(Mario *evt){
  double sumwf[300], zerocross;
  float s0, s1;
  int i = 36, j, k;
  int thflg = 0;

  for (j = 0; j<cfdelay; j++) {
    s1 = cf * (double)evt->wf[i][j];
    sumwf[j] = s1;
  }
  for (j = cfdelay; j < 300; j++) {
    s0 = (double)evt->wf[i][j - cfdelay];
    s1 = cf * (double)evt->wf[i][j];
    sumwf[j] = s0 + s1;
  }

  for(j= 10; j < 300 - 10; j++) {
    //if(thflg == 1 && sumwf[j] < 0. && sumwf[j+1] > 0.){
    if(thflg == 1) {
      s0 = 0.;
      s1 = 0.;
      for(k = 0; k < 5; k++){
	s0 += sumwf[j - k]/5.;
	s1 += sumwf[j + 1 + k]/5.;
      }
      if(s0 < 0. && s1 > 0. && s1 - s0 > sumwf[j+1] - sumwf[j])  break;
    }else if(thflg == 0 && sumwf[j] < threshold){
      thflg = 1;
    }
  }
  //zerocross =  (double) j + sumwf[j] / (sumwf[j]-sumwf[j+1]); 
  zerocross =  (double) j -2. + 5. * s0 / (s0 - s1); 

  if(thflg == 1){
    return zerocross;
  } else {
    return -1000.;
  }
}


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
  wevt->ccEnergy = weight * evt0->ccEnergy + (1 - weight) * evt1->ccEnergy;

  return wevt;
}

Mario *avg_thoffset(struct evtList *x, int baselineFlag, double refoffset, char *suffix) {
  Mario *avgEvt, *evt;
  float avg;
  int i, j, k, n = 0;
  int Deltat;
  int refwf[300];
  int thesegid = 36;
  FILE *foffset;
  char outname[100];

  avgEvt = calloc(1, sizeof(Mario));
  memset(scratch, 0, 37 * 300 * sizeof(float));

  sprintf(outname, "./cevtout/offset/dummy/%s_cfd%s.csv", x->filename, suffix);
  //foffset = fopen("output/offset.csv", "w");
  foffset = fopen(outname, "w");
  fprintf(foffset, "\n RefOffset = %f \n evtnum, Offset of the evt, Deltat\n", refoffset);
  
  for (i = 0; i < x->numEvts; i++) {
    evt = x->rawEvts + i;
    memset(refwf, 0, 300 * sizeof(int));
    for (k = 0; k < 300; k++) {
      refwf[k] = (float) evt->wf[thesegid][k];
    }

    if (baselineFlag != 0) {
      avg = 0.;
      for (j = 0; j < 30; j++) { // According to the article, 40 points are used for the baseline estimation.. RT Jan11
	avg += refwf[j];
      }
      avg /= 30.;
      for (j = 0; j < 300; j++) {
	refwf[j] -= avg;
      }
    }
    
    for (j = 0; j < 300; j++) {
      avgEvt->wf[thesegid][j] = (short) refwf[j];
    }
    Deltat = (int) (refoffset - cft(avgEvt));
    assert(Deltat*Deltat < 90000);

    //Now fill scrach[][] again with applying time offset.
    if(Deltat * Deltat >100.){ // strange event should be discarded at this point
      //x->numEvts =  x->numEvts - 1; <- Bug
      n++; // Add dummy int
      printf("Evt. %i of %s\t Ref offset: %f, offset: %f, Delta T: %i -> Discarded. \n",i , x->filename , refoffset,cft(avgEvt),Deltat); // Print discarded events
      continue;
    }
    //printf("%f, %f, %i\n", refoffset, cft(avgEvt), Deltat); // debug
    fprintf(foffset, "%i, %f, %i\n", i, cft(avgEvt), Deltat); 
	

    if(Deltat < 0) {
      for (j = 0; j < 37; j++) {
	for (k = 0; k < 300 + Deltat; k++) {
	  scratch[j][k] += (float) evt->wf[j][k - Deltat];
	}
	if (baselineFlag != 0) { // Because the triger timing is shifted for each event, besline subtraction should be done evt by evt. RT Jan13
	  avg = 0.;
	  for (k = 0; k < 30; k++) {
	    avg +=  (float) evt->wf[j][k - Deltat];
	  }
	  avg /= 30.;
	  for (k = 0; k < 300 + Deltat; k++) {
	    scratch[j][k] -= avg;
	  }
	}
      }
    } else {
      for (j = 0; j < 37; j++) {
	for (k = Deltat; k < 300; k++) {
	  scratch[j][k] += (float) evt->wf[j][k - Deltat];
	}
	if (baselineFlag != 0) {
	  avg = 0.;
	  for (k = Deltat; k < 30 + Deltat; k++) {
	    avg +=  (float) evt->wf[j][k - Deltat];
	  }
	  avg /= 30.;
	  for (k = Deltat; k < 300; k++) {
	    scratch[j][k] -= avg;
	  }
	}
      }
    }
    
    avgEvt->ccEnergy += evt->ccEnergy;
    for (j = 0; j < 36; j++) {
       avgEvt->segEnergy[j] += evt->segEnergy[j];
    }
  }

  //Divide by "accepted num evts"
  x->numEvts = x->numEvts - n;
  for (j = 0; j < 37; j++) {
    for (k = 0; k < 300; k++) {
      scratch[j][k] /= (float) x->numEvts;
      avgEvt->wf[j][k] = (short) scratch[j][k];
    }
  }

  avgEvt->ccEnergy /= (float) x->numEvts;
  for (i = 0; i < 36; i++) { avgEvt->segEnergy[i] /= (float) x->numEvts; }

  fclose(foffset);

  return avgEvt;
}

Mario *avg(struct evtList *x, int baselineFlag) {// Old fcn

  //struct evtList *y; // not used
  Mario *avgEvt, *evt;
  float avg;
  int i, j, k;

  avgEvt = calloc(1, sizeof(Mario));
  memset(scratch, 0, 37 * 300 * sizeof(float));

  for (i = 0; i < x->numEvts; i++) {
    evt = x->rawEvts + i;
    for (j = 0; j < 37; j++) {
      for (k = 0; k < 300; k++) {
        scratch[j][k] += (float) evt->wf[j][k];
      }
    }
  }

  if (baselineFlag != 0) {
    for (i = 0; i < 37; i++) {
      avg = 0.;
      for (j = 0; j < 30; j++) { // According to the article, 40 points are used for the baseline estimation.. RT Jan11
        avg += scratch[i][j];
      }
      avg /= 30.;
      for (j = 0; j < 300; j++) {
        scratch[i][j] -= avg;
      }
    }
  }

  for (j = 0; j < 37; j++) {
    for (k = 0; k < 300; k++) {
      scratch[j][k] /= (float) x->numEvts;
      avgEvt->wf[j][k] = (short) scratch[j][k];
    }
  }

  for (i = 0; i < x->numEvts; i++) {
    evt = x->rawEvts + i; //Old one was "evt = x->rawEvts;" which doesn't incriment. Jan11
    avgEvt->ccEnergy += evt->ccEnergy;
    for (j = 0; j < 36; j++) {
      avgEvt->segEnergy[j] += evt->segEnergy[j];
    }
  }
  avgEvt->ccEnergy /= (float) x->numEvts;
  for (i = 0; i < 36; i++) { avgEvt->segEnergy[i] /= (float) x->numEvts; }

  return avgEvt;
}

float cc_avg_cft(struct evtList *x, int baselineFlag) { // Get the constant fraction timing of Center contact channel. 36th
  Mario *avgEvt, *evt;
  float avg;
  int i, j, k;

  avgEvt = calloc(1, sizeof(Mario));
  memset(scratch, 0, 37 * 300 * sizeof(float));

  for (i = 0; i < x->numEvts; i++) {
    evt = x->rawEvts + i;
    for (k = 0; k < 300; k++) {
      scratch[36][k] += (float) evt->wf[36][k];
    }
  }

  if (baselineFlag != 0) {
    avg = 0.;
    for (j = 0; j < 30; j++) { // According to the article, 40 points are used for the baseline estimation.. RT Jan11
      avg += scratch[36][j];
    }
    avg /= 30.;
    for (j = 0; j < 300; j++) {
      scratch[36][j] -= avg;
    }
  }

  for (k = 0; k < 300; k++) {
    scratch[36][k] /= (float) x->numEvts;
    avgEvt->wf[36][k] = (short) scratch[36][k];
  }

  return cft(avgEvt);

}


int cevt(double ratio, char *datname1, int runn1, int segn1, char *datname2, int runn2, int segn2){

  FILE *fin, *fou, *ftr;
  struct gebData ghdr, defaulthdr = {100, sizeof(Mario), 0ll};
  Mario *evt;
  char evtlabel[] = {'a', 'b', 's'}, seglabel[] = {'n', 'l', 'r', 'u', 'd', 'c'};
  char suffix[50];
  sprintf(suffix, "_Run%iSeg%i_Run%iSeg%i", runn1, segn1, runn2, segn2);
  //int i, j, k, num, seg, baselineFlag = 1, CFDFlag = 0;
  /*  char *filename;
  int seg;
  Mario *rawEvts;
  int numEvts;
  runList[0]. = (struct evtList){datname1, segn1, 0, 0};
  runList[1] = (struct evtList){datname2, segn2, 0, 0};
  */

  for (i = 0; i < 2; i++) {
    runList[i].rawEvts = malloc(MAX_EVT * sizeof(Mario));
  }
  runList[0].filename = datname1;
  runList[0].seg = segn1;
  runList[0].numEvts = 0;
  runList[1].filename = datname2;
  runList[1].seg = segn2;
  runList[1].numEvts = 0;


  pList[0].seg = segn1;
  pList[1].seg = segn2;
  pList[2].seg = segn1;

  //printf("%i\n",runList[i].numEvts);

  for (i = 0; i < 2; i++) {
    fin = fopen(runList[i].filename, "r");
    if (fin == 0) { fprintf(stderr, "could not open file %s\n", runList[i].filename); exit(1);}
    while ((runList[i].numEvts < MAX_EVT) && (fread(&ghdr, sizeof(struct gebData), 1, fin) == 1)) {
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

  if(CFDFlag == 1) {
    thoffset[0] = cc_avg_cft(runList + 0, baselineFlag); // Deduce the reference time

    pList[0].rawEvts = avg_thoffset(runList + 0, baselineFlag, thoffset[0], suffix);
    pList[1].rawEvts = avg_thoffset(runList + 1, baselineFlag, thoffset[0], suffix); //use theoffset[0] not [1]
  } else {
    pList[0].rawEvts = avg(runList + 0, baselineFlag);
    pList[1].rawEvts = avg(runList + 1, baselineFlag);
  }

  pList[2].rawEvts = wsum(pList[0].rawEvts, pList[1].rawEvts, ratio);
  
  char outname[100];
  CFDFlag == 1 ? sprintf(outname, "./cevtout/cevt%s_cfd.dat",suffix) : sprintf(outname, "./cevtout/cevt%s.dat",suffix);
  fou = (appendFlag == 1) ? fopen(outname, "a") : fopen(outname, "w");
  if (fou == 0) { fprintf(stderr, "could not open file %s\n",outname); exit(1);}

  for (i = 0; i < 3; i++) {
    if (sumOnlyFlag == 1 && i != 2) { continue; }
    assert( 1 == fwrite(&defaulthdr, sizeof(struct gebData), 1, fou));
    assert( 1 == fwrite(pList[i].rawEvts, sizeof(Mario), 1, fou));
  }
  fclose(fou);

  fprintf(frunlistout, "%s %d %d\n", outname, 1000*runn1 + runn2, segn1);
  printf("\n%s %d %d\n", outname, 1000*runn1 + runn2, segn1);
	
  printf("list 1: %d evts, list 2: %d evts\n", runList[0].numEvts, runList[1].numEvts);

  if (appendFlag == 0) {
    CFDFlag == 1 ? sprintf(outname, "./cevtout/tr/tr_cfd%s.csv",suffix) : sprintf(outname, "./cevtout/tr/tr%s.csv",suffix);
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
    close(ftr);
  }
  return 0;
}

int main(int argc, char **argv) {

  struct option opts[] = {{"verbose", no_argument, 0, 'v'},
                          {"sum-only", no_argument, 0, 's'},
                          {"append", no_argument, 0, 'a'},
                          {"ratio", required_argument, 0, 'r'}, 
                          {"cfd", no_argument, 0, 'c'}, //Use CFD
                          {"filelist", required_argument, 0, 'l'},//Filelist
			  { 0, 0, 0, 0}};
  /*  int verboseFlag = 0, sumOnlyFlag = 0, appendFlag = 0, fileListFlag = 0;
      int i, j, k, num, seg, baselineFlag = 1, CFDFlag = 0;*/
  double ratio = 0.5;   // default
  char ch;
  char *inputFileList, s[80];;
  char *datname[20];
  int runn[20], segn[20];
  int cnt, numTok, numRuns;
  FILE *fl;
  int runA, runB;

  inputFileList = malloc(80 * sizeof(char));
  
  while ((ch = getopt_long(argc, argv, "vsar:cl:", opts, 0)) != -1) {
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
    case 'l': inputFileList = optarg;
              fileListFlag = 1;
              break;
    default: fprintf(stderr, "usage: cevt [-vsar:cl:]\n");
             exit(1);
    }
  }

  /*
  for (i = 0; i < 20; i++) {
    runListTemp[i].rawEvts = malloc(maxevt * sizeof(Mario));
    runListTemp[i].filename = malloc(80);
    runListTemp[i].numEvts = 0;
  }
  */
  
  char outname[100];
  if (fileListFlag == 1) {
    CFDFlag ? sprintf(outname,"cevtrunlist_cfd.txt") : sprintf(outname,"cevtrunlist.txt");
    frunlistout = fopen(outname, "w");
    fl = fopen(inputFileList, "r");
    if (fl == 0) { fprintf(stderr, "could not open file %s\n", inputFileList); exit(1);}
    cnt = 0;
    while (cnt < 20 && fgets(s, 80, fl) != 0) {
      datname[cnt] = malloc(80 * sizeof(char));
      numTok = sscanf(s, "%s %d %d", datname[cnt], &runn[cnt], &segn[cnt]);
      if (numTok != 3) { fprintf(stderr, "wrong fmt - line %d, %s\n", cnt + 1, inputFileList); exit(1);}
      cnt++;
    }
    /*else {
      strncpy(runListTemp[0].filename, inputFile, 80);
      runListTemp[0].run = 1, runListTemp[0].seg = 15; // defaults
      cnt = 1;
      }*/
    numRuns = cnt;
    //numRuns = 2;
    for (runA = 0; runA < numRuns - 1; runA++){
      for (runB = runA + 1 ; runB < numRuns; runB++){
	//printf("A:%i %s %d %d\n", runA, datname[runA], runn[runA], segn[runA]);
      	//printf("B:%i %s %d %d\n", runB, datname[runB], runn[runB], segn[runB]);
      	cevt(0.5, datname[runA], runn[runA], segn[runA], datname[runB], runn[runB], segn[runB]);
      }
    }
  } else {
    exit(1);
  }
  return 0;
}
