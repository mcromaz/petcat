/*** diagcnt.c: counters for preprocessor and decomp routines */

#include <stdio.h>
#include <assert.h>
#include "diagcnt.h"

int cntbits(long long int x) {

  int i, cnt;
  for (i = 0, cnt = 0; i < 63; i++) {
    if ((x & (1ll << i)) != 0) {
      cnt++;
    }
  }
  assert(cnt >= 0 && cnt < 64);
  if (cnt == 0) printf("++ cntbits returned 0\n");
  return cnt;
}

void ppPreprocCnt(preprocCnt *x) {
  int i;
  fprintf(stdout, "--preprocessor counters\n");
  fprintf(stdout, "numMode3Evts = %d\n", x->numMode3Evts);
  fprintf(stdout, "numPartial = %d\n", x->numPartial);
  fprintf(stdout, "numFull = %d\n", x->numFull);
  fprintf(stdout, "pp58 = %d\n", x->pp58);
  for (i = 0; i < 50; i++) {
    fprintf(stdout, "[%5d] ", (x->full)[i]);
    if (((i + 1) % 10) == 0) fprintf(stdout, "\n");
  }
  fflush(stdout);
}

void ppPostprocCnt(postprocCnt *x) {
  int i;
  fprintf(stdout, "--postprocessor counters\n");
  fprintf(stdout, "oneInteractionPoint(to start):\n");
  for (i=0; i<36; i++) {
    fprintf(stdout, "[%5d] ", (x->oneInt)[i]);
    if (((i+1)%6) == 0) fprintf(stdout, "\n");
  }
  fprintf(stdout, "twoInteractionPoints(to start):\n");
  for (i=0; i<36; i++) {
    fprintf(stdout, "[%5d] ", (x->twoInt)[i]);
    if (((i+1)%6) == 0) fprintf(stdout, "\n");
  }
  fprintf(stdout, "threeInteractionPoints(to start):\n");
  for (i=0; i<36; i++) {
    fprintf(stdout, "[%5d] ", (x->threeInt)[i]);
    if (((i+1)%6) == 0) fprintf(stdout, "\n");
  }
  fprintf(stdout, "> threeInteractionPoints(to start):\n");
  for (i=0; i<36; i++) {
    fprintf(stdout, "[%5d] ", (x->manyInt)[i]);
    if (((i+1)%6) == 0) fprintf(stdout, "\n");
  }
  fprintf(stdout, "combinedDiffSeg = %d\n", x->combinedDiffSeg);
  fprintf(stdout, "combinedWithinSeg = %d\n", x->combinedWithinSeg);
  fprintf(stdout, "combined3to2:\n");
  for (i=0; i<36; i++) {
    fprintf(stdout, "[%5d] ", (x->combined3to2)[i]);
    if (((i+1)%6) == 0) fprintf(stdout, "\n");
  }
  fprintf(stdout, "combined3to1:\n");
  for (i=0; i<36; i++) {
    fprintf(stdout, "[%5d] ", (x->combined3to1)[i]);
    if (((i+1)%6) == 0) fprintf(stdout, "\n");
  }
  fprintf(stdout, "combined2to1:\n");
  for (i=0; i<36; i++) {
    fprintf(stdout, "[%5d] ", (x->combined2to1)[i]);
    if (((i+1)%6) == 0) fprintf(stdout, "\n");
  }
  fflush(stdout);
}

void initializePostCnt(postprocCnt *x) {
  int i;
  for (i=0; i<36; i++) {
    x->oneInt[i] = 0;
    x->twoInt[i] = 0;
    x->threeInt[i] = 0;
    x->manyInt[i] = 0;
    x->combined3to2[i] = 0;
    x->combined3to1[i] = 0;
    x->combined2to1[i] = 0;
  }
  x->combinedDiffSeg = 0;
  x->combinedWithinSeg = 0;
}
