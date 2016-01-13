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

//for calculating the timing of CFD Jan12 RT
double cf = -0.5;//should be negative value
int cfdelay = 10;
double thoffset[2] = {0.0};
int threshold = -50;

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
    if(thflg == 1 && sumwf[j] < 0. && sumwf[j+1] > 0.){
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

Mario *avg_thoffset(struct evtList *x, int baselineFlag, double refoffset) {
  Mario *avgEvt, *evt;
  float avg;
  int i, j, k;
  int Deltat;
  int refwf[300];
  int thesegid = 36;

  avgEvt = calloc(1, sizeof(Mario));
  memset(scratch, 0, 37 * 300 * sizeof(float));

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
      x->numEvts =  x->numEvts - 1;
      printf("Evt. %i of %s\t Ref offset: %f, offset: %f, Delta T: %i -> Discarded. \n",i , x->filename , refoffset,cft(avgEvt),Deltat); // Print discarded events
      continue;
    }
    //printf("%f, %f, %i\n", refoffset, cft(avgEvt), Deltat); // debug

    if(Deltat < 0) {
      for (j = 0; j < 37; j++) {
	for (k = 0; k < 300 + Deltat; k++) {
	  scratch[j][k] += (float) evt->wf[j][k - Deltat];
	}
      }
    } else {
      for (j = 0; j < 37; j++) {
	for (k = Deltat; k < 300; k++) {
	  scratch[j][k] += (float) evt->wf[j][k - Deltat];
	}
      }
    }
  
    avgEvt->ccEnergy += evt->ccEnergy;
    for (j = 0; j < 36; j++) {
       avgEvt->segEnergy[j] += evt->segEnergy[j];
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

  //Divide by "accepted num evts"
  for (j = 0; j < 37; j++) {
    for (k = 0; k < 300; k++) {
      scratch[j][k] /= (float) x->numEvts;
      avgEvt->wf[j][k] = (short) scratch[j][k];
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

int main(int argc, char **argv) {

  FILE *fin, *fou, *ftr;
  struct gebData ghdr, defaulthdr = {100, sizeof(Mario), 0ll};
  Mario *evt;
  char evtlabel[] = {'a', 'b', 's'}, seglabel[] = {'n', 'l', 'r', 'u', 'd'};
  struct option opts[] = {{"verbose", no_argument, 0, 'v'},
                          {"sum-only", no_argument, 0, 's'},
                          {"append", no_argument, 0, 'a'},
                          {"ratio", required_argument, 0, 'r'},
                          { 0, 0, 0, 0}};
  int verboseFlag = 0, sumOnlyFlag = 0, appendFlag = 0;
  int i, j, k, num, seg, baselineFlag = 1;
  double ratio = 0.5;   // default
  char ch;

  for (i = 0; i < 2; i++) {
    runList[i].rawEvts = malloc(MAX_EVT * sizeof(Mario));
  }

  while ((ch = getopt_long(argc, argv, "vsar:", opts, 0)) != -1) {
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
    default: fprintf(stderr, "usage: cevt [-vsar:]\n");
             exit(1);
    }
  }

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

  thoffset[0] = cc_avg_cft(runList + 0, baselineFlag);
  thoffset[1] = cc_avg_cft(runList + 1, baselineFlag); // It was necessary to be executed to adjust baseline

  pList[0].rawEvts = avg_thoffset(runList + 0, baselineFlag, thoffset[0]);
  pList[1].rawEvts = avg_thoffset(runList + 1, baselineFlag, thoffset[0]); //use theoffset[0] not [1]

  pList[2].rawEvts = wsum(pList[0].rawEvts, pList[1].rawEvts, ratio);

  fou = (appendFlag == 1) ? fopen("cevt.dat", "a") : fopen("cevt.dat", "w");
  if (fou == 0) { fprintf(stderr, "could not open file cevt.dat\n"); exit(1);}
  for (i = 0; i < 3; i++) {
    if (sumOnlyFlag == 1 && i != 2) { continue; }
    assert( 1 == fwrite(&defaulthdr, sizeof(struct gebData), 1, fou));
    assert( 1 == fwrite(pList[i].rawEvts, sizeof(Mario), 1, fou));
  }
  fclose(fou);

  printf("list 1: %d evts, list 2: %d evts\n", runList[0].numEvts, runList[1].numEvts);

  if (appendFlag == 0) {
    ftr = fopen("tr.csv", "w");
    if (ftr == 0) { fprintf(stderr, "could not open file tr.csv\n"); exit(1);}
    fprintf(ftr, "evt, seg, ch, val\n");
    for (i = 0; i < 3; i++) {
      if (sumOnlyFlag == 1 && i != 2) { continue; }
      for (j = 0; j < 5; j++) {  // 'n', 'l', 'r', 'u', 'd'
          //seg = (j == 0) ? runList[i].seg : neigh[runList[i].seg][j - 1];
          //evt = runList[i].rawEvts + 7;  // take first evt
          seg = (j == 0) ? pList[i].seg : neigh[pList[i].seg][j - 1];
          evt = pList[i].rawEvts;  // take 1st evt
          for (k = 0; k < 300; k++) {
            fprintf(ftr, "%c,%c,%d,%d\n", evtlabel[i], seglabel[j], k + 1, evt->wf[seg][k]);
        }
      }
    }
    fclose(ftr);
  }
  return 0;
}
