/*** mode3io: routines to read mode3 data */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <netinet/in.h>
#include "petcat.h"
#include "mode3io.h"
#include "wrap.h"

unsigned short int sbuf[SBUF_LEN];

int get_evt_len(FILE *fin) {
  /*** gives the length (in short ints) of the first mode3 evt including 0xaaaaaaaa sentinnel */
  unsigned short int *sbuf;
  static int sbuf_len = 16384;
  int num, i;
	
  assert(fin != 0);
  sbuf = Calloc(sbuf_len, sizeof(unsigned short int));
  num = fread(sbuf, sizeof(unsigned short int), sbuf_len, fin);
  assert(num == sbuf_len);
  rewind(fin);
	
  for (i = 0; i < sbuf_len; i += 2) {
    if (sbuf[i] == 0xaaaa && sbuf[i + 1] == 0xaaaa) {
    return 2 * (ntohs(sbuf[i + 2]) & 0x7ff) + 2;  /* inclusive event length in short ints */
    }
  }
  return -2;	/* mode 3 evt not found in first part of file */
}

/*** fill: fill a buffer with non-fragmented mode 3 events */

int fill(unsigned short int *sbuf, FILE *fin) {

  static int fraglen = 0;
  static unsigned short int fragbuf[MAX_EVT3_LEN];
  int totlen;
  int i;

  assert(sbuf != 0);
  assert(fin != 0);
  assert(fraglen >= 0);

  if (fraglen == 0) {
    (void) memcpy(sbuf, fragbuf, fraglen);
  }

  totlen = fread(sbuf + fraglen, sizeof(unsigned short int), SBUF_LEN - fraglen, fin);
  if (totlen <= 0) { /* assume eof */
    return 0;
  }
  /* copy fragment at end of buf to frag */
  for (i = totlen - 1; i > 0; i--) {
   if (sbuf[i - 1] == 0xaaaa && sbuf[i] == 0xaaaa) {
     /* copy fragment */
     /* return len */
   }  
  }
  fprintf(stderr, "fill: no sentinel found in mode 3 buffer");
  return -1;
}

mode3stream *init3stream(char *filename, int id) {

  mode3stream *s;

  assert(filename != 0);
  assert(id >= 0);

  s = Calloc(sizeof(mode3stream), 1);
  s->fin = fopen(filename, "r");
  if (s->fin == 0) {
    fprintf(stderr, "init3stream: unable to open file %s\n", filename);
    return 0;
  }
  s->sbuf = Malloc(sizeof(unsigned short int) * SBUF_LEN);
  s->len = fill(s->sbuf, s->fin);
  s->id = id;
  
  return s;
}

unsigned short int *next(mode3stream m3s) {

  /* jump event length, is there an aaaa, if so return it */
  /* start from where you are an search, report error */
  return 0;
}

int read3(FILE *fin, int holenumreqd, int xtalreqd, unsigned short int *obuf, mode3Cnt *cnt) {

  struct gebData hdr;
  int numhdr = 0, num3evts = 0, numbadhdr = 0, numinobuf = 0, evts, num;
  unsigned short int *s, *o;
  int xtal, holenum, inclevtlen, i, j;	/* inclevtlen - length of a single evt in bytes including "aaaaaaaa" */
  char *c;
	
  assert(fin != 0);

  s = sbuf;
  o = obuf;
	
  while (fread(&hdr, sizeof(struct gebData), 1, fin) == 1) {
    cnt->numHdr++;
    //gh_log3(&hdr);
    if (hdr.type != 2) {	/* skip, n.b. mode 3 corresponds to type 2 */
      fseek(fin, hdr.length, SEEK_CUR);
      continue;
    }
	
    cnt->numMode3Hdr++;
    assert(hdr.length <= SBUF_LEN * 2);
    num = fread(sbuf, sizeof(char), hdr.length, fin);
    if (num != hdr.length) {	/* fragmented data or messed-up hdr.length */
      fprintf(stderr, "fragment data or erroneous header length (num = %d, hdr.length = %d)\n", num, hdr.length);
      break;
    }
		
    inclevtlen = 4 * (ntohs(sbuf[2]) & 0x7ff) + 4;  /* incl evt length in bytes */
    evts = hdr.length / inclevtlen;
    assert((hdr.length % inclevtlen) == 0);
    cnt->numMode3Evts += evts;

    for (i = 0, s = sbuf; i < evts; i++) {
      holenum = ((ntohs(s[3])) >> 8) & 0x1f;
      xtal = ((ntohs(s[3])) >> 6) & 0x3;
      assert(holenum >= 0 && holenum < 32);
      assert(xtal >= 0 && xtal < 4);
      assert(((holenum << 2) + xtal) < 256);
      (cnt->holextal)[(holenum << 2) + xtal]++;
      if ((holenum == holenumreqd) && (xtal == xtalreqd)) {
        cnt->numMode3AcceptedEvts++;
        numinobuf++;
        /*
        for (j = 0; j < inclevtlen / 2; j++) {
          o[j] = ntohs(s[j]);  
        }
        */
        swab((char *) s, (char *) o, inclevtlen);
        o += inclevtlen / 2;
      }
      cnt->numObuf++;
      
      if (o + inclevtlen / 2 >= obuf + OBUF_LEN) {
        cnt->numGoodButIgnored += evts - i;
        //fprintf(stdout, "read3() - len = %d\n", o - obuf);
        return (o - obuf); /* no room for next evt, just bail */
      }
      s += inclevtlen / 2;
    } 	/* end for */
  }	/* end while */

  return (o - obuf);	/* length of valid data in obuf, again in shorts */			
}

void ppMode3Cnt(mode3Cnt *x) {

  int i, n0, n1, n2, n3;
  fprintf(stdout, "numHdr = %d\n", x->numHdr);
  fprintf(stdout, "numMode3Hdr = %d\n", x->numMode3Hdr);
  fprintf(stdout, "numMode3Evts = %d\n", x->numMode3Evts);
  fprintf(stdout, "numMode3AcceptedEvts = %d\n", x->numMode3AcceptedEvts);
  fprintf(stdout, "numObuf = %d\n", x->numObuf);
  fprintf(stdout, "numGoodButIgnored = %d\n", x->numGoodButIgnored);

  for (i = 0; i < 32; i++) {
    n0 = (x->holextal)[i * 4];
    n1 = (x->holextal)[i * 4 + 1];
    n2 = (x->holextal)[i * 4 + 2];
    n3 = (x->holextal)[i * 4 + 3];
    if (n0 != 0 || n1 != 0 || n2 != 0 || n3 != 0) {
      fprintf(stdout, "hole %2d, [%8d] [%8d] [%8d] [%8d]\n", i, n0, n1, n2, n3);
    }
  
  }
}

long long int readMarioFormat(FILE *fin, Mario *mario, mode3Cnt *cnt) {

  struct gebData hdr;
  int num = 0;
  
  while (fread(&hdr, sizeof(struct gebData), 1, fin) == 1) {
    cnt->numHdr++;

    if (hdr.type != 100) {	/* skip */
      fseek(fin, hdr.length, SEEK_CUR);
      continue;
    }    
    cnt->numMode3Hdr++;
    num = fread(mario, 1, hdr.length, fin);
    break;
  }
  if (num > 0) {
    return hdr.timestamp;
  } else {
    return num;
  }

}
