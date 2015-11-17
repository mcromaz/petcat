#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "petcat.h"
#include "vegcat.h"
#include "neigh.h"
#define MAX_EVT 100

int main(int argc, char **argv) {

  FILE *fin, *ftr;
  struct {
    char *filename;
    int seg;
    Mario *rawEvts;
    int numEvts;
  } runList[2] = {{"../coinc-data/WFOUT-Run0114Segment15.dat", 15, 0, 0},
                  {"../coinc-data/WFOUT-Run0114Segment9.dat", 9, 0, 0}};
  struct gebData ghdr;
  Mario *evt;
  char evtlabel[] = {'a', 'b', 's'}, seglabel[] = {'n', 'l', 'r', 'u', 'd'};
  int i, j, k, num, seg;

  for (i = 0; i < 2; i++) {
    runList[i].rawEvts = malloc(MAX_EVT * sizeof(Mario));
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

  printf("list 1: %d evts, list 2: %d evts\n", runList[0].numEvts, runList[1].numEvts);

  ftr = fopen("tr.csv", "w");
  if (fin == 0) { fprintf(stderr, "could not open file tr.csv\n"); exit(1);}

  fprintf(ftr, "evt, seg, ch, val\n");
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 5; j++) {  // 'n', 'l', 'r', 'u', 'd'
        seg = (j == 0) ? runList[i].seg : neigh[runList[i].seg][j - 1];
        evt = runList[i].rawEvts;  // take first evt
        for (k = 0; k < 300; k++) {
          fprintf(ftr, "%c, %c, %d, %d\n", evtlabel[i], seglabel[j], k + 1, evt->wf[seg][k]);
      }
    }
  }
}
