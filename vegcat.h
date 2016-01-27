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
} cfg = {109, 0, "../coinc/q4a8_basis_xt.dat", "../coinc/detmap_Q4pos4_CC09_Ecal.txt", "../coinc/filter.txt",
	 "../coinc/tr_gain_Q4pos4_CC09.txt", "../coinc/q4a8_xtalk_pars_in.txt"}; //RT mod

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
