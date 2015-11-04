/*** diagcnt.h: counters for preprocessor and decomp routines */

#ifndef _diagcnt_h
#define _diagcnt_h

#include "lh.h"

typedef struct {
  int numMode3Evts;
  int numPartial;
  int numFull;
  int pp58;
  int full[64];
  lh *ener_cc;
  int seg[36][4096];
} preprocCnt;

void ppPreprocCnt(preprocCnt *cnt);
int cntbits(long long int);

typedef struct {
  int combinedDiffSeg;
  int combinedWithinSeg;
  int combined3to2[36];
  int combined3to1[36];
  int combined2to1[36];
  int oneInt[36];
  int twoInt[36];
  int threeInt[36];
  int manyInt[36];
} postprocCnt;

void ppPostprocCnt(postprocCnt *cnt);

#endif /* _diagcnt_h */
