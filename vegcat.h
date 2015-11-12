#ifndef _vegcat_h
#define _vegcat_h

struct config {
  int holenum;
  int xtalnum;
  char *basisName;
  char *detMapName;
  char *filterName;
  char *trGainName;
  char *xTalkParsName;
} cfg = {109, 0, "q4a8_basis_xt.dat", "detmap_Q4pos4_CC09.txt", "filter.txt",
      "tr_gain_Q4pos4_CC09.txt", "q4a8_xtalk_pars_in.txt"};

struct gebData {
  int type;
  int length;
  long long int timestamp;
};

typedef struct {
  float ccEnergy;
  float segEnergy[36];
  float pad;
  short wf[37][300];
} Mario;

int preProcessMario(Mario *mario, Event_Signal *event, preprocCnt *diagcnt);

#endif
