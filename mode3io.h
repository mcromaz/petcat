
/*** mode3io.h */

#include <stdio.h>

#define SBUF_LEN (32 * 1024 * 1024)
#define OBUF_LEN (32 * 1024 * 1024)
#define MAX_EVT3_LEN 2048 /* placeholder only */

struct gebData {
  int type;
  int length;
  long long int timestamp;
};

typedef struct {
  FILE *fin;
  unsigned short int *sbuf;
  unsigned short int *curr;
  int id;
  int len;
} mode3stream;

typedef struct {
  int numHdr;
  int numMode3Hdr;
  int numMode3Evts;
  int numMode3AcceptedEvts;
  int numObuf;
  int numGoodButIgnored;
  int holextal[256];
} mode3Cnt;

typedef struct {
  float ccEnergy;
  float segEnergy[36];
  float pad;
  short wf[37][300];
} Mario;

/* prototypes */
mode3stream *init3stream(char *filename, int id);
long long int readMarioFormat(FILE *fin, Mario *mario, mode3Cnt *cnt);
int read3(FILE *fin, int holenumreqd, int xtalreqd, unsigned short int *obuf, mode3Cnt *cnt);
void gh_log3(struct gebData *hdr);
void ppMode3Cnt(mode3Cnt *cnt);
