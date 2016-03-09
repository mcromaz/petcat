#ifndef _vegcat_h
#define _vegcat_h

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

typedef struct {
  float cc_ener;
  float tr_len;
  float seg_ener[36];
  float tr[37][1024];
} m3eb; //mode 3, event-built

int preProcessMario(Mario *mario, Event_Signal *event, preprocCnt *diagcnt);

/* single interaction */

typedef struct sdiag {
  int stat;
} sdiag;

void sint_init(char *basisName, char *trGainName);
struct crys_intpts *sint(Mario *m, sdiag *sd);

#endif
