/*** fal.c: access functions for a given digitizer data format */

#include <stdio.h>
#include "petcat.h"

struct fmt_defn get = {0xaaaaaaaa, ener_nd, time_nd, vec_nd, ch_nd, tr_nd, tr_len_nd, evt_len_nd};

int ener_od(unsigned short int *e) {

    int ener;

    ener = (int) ((e[7] & 0x3f) << 16) + e[4];
    if (ener > 2097151) {
      ener -= 4194304;
    }
    ener = (ener < 0) ? -ener : ener;

    return ener;
}

long long int time_od(unsigned short int *e) {

  return e[3] + (((unsigned long long int) e[2]) << 16) + (((unsigned long long int) e[5]) << 32);
 
}

int vec_od(unsigned short int *e) {
  return e[1] >> 3;
}

int ch_od(unsigned short int *e) {
  return e[1] & 0x7;
}

int tr_od(int *tr, int tr_len, int inv, unsigned short int *e) {

  int act_tr_len;
  unsigned short int *raw_tr;
  int i;

  act_tr_len = get.tr_len(e);
  /*printf("trace length = %d\n", trace_len); */
  if (tr_len > MAX_TR_LEN)  {
    return 1;
  }
  raw_tr = e + 12;
  for (i = 0; i < act_tr_len; i += 2) {
    tr[i] = -(raw_tr[i + 1] & 2047);
    if (raw_tr[i + 1] & 2048)
      tr[i] += 2048;
    if (inv != 0)                                /* invert */
      tr[i] = -tr[i];

    tr[i + 1] = -(raw_tr[i] & 2047);
    if (raw_tr[i] & 2048)
      tr[i + 1] += 2048;
    if (inv != 0)
      tr[i + 1] = -tr[i + 1];
  }
  return 0;
}

int tr_len_od(unsigned short int *e) {
  return e[0] * 2 - 14;
}

int evt_len_od(unsigned short int * e) {
  return  2 * e[0] + 2;
}

int ener_nd(unsigned short int *e) {

    int ener;

    ener = (int) ((e[7] & 0x1ff) << 16) + e[4];
    if (ener > 16777216) {
      ener -= 33554432;
    }
    ener = (ener < 0) ? -ener : ener;

    return ener;
}
int vec_nd(unsigned short int *e) {
  return (e[0] & 0xf800) >> 11;
}

int ch_nd(unsigned short int *e) {
  return e[1] & 0xf;
}

long long int time_nd(unsigned short int *e) {
  return e[3] + (((unsigned long long int) e[2]) << 16) + (((unsigned long long int) e[5]) << 32);
}

int tr_nd(int *tr, int tr_len, int inv, unsigned short int *e) {

  int act_tr_len, i; 
  short int *raw_tr;

  act_tr_len = get.tr_len(e);
  //printf("actual trace length = %d\n", act_tr_len);
  if (tr_len > MAX_TR_LEN)  {
    printf("bad tr len\n");
    return 1;
  }
  raw_tr = (short int *) (e + 14);
  for (i = 0; i < tr_len; i += 2) {
    tr[i] = raw_tr[i + 1];
    tr[i + 1] = raw_tr[i];
    //printf("tr[%d] = %d\n", i, tr[i]);
    //printf("tr[%d] = %d\n", i + 1, tr[i + 1]);
  }
  //printf("last i = %d, tr[%d] = %d, tr[%d] = %d\n", i, i, tr[i], i + 1, tr[i + 1]);
  //printf("last i = %d, raw_tr[%d] = %d, raw_tr[%d] = %d\n", i, i, raw_tr[i], i + 1, raw_tr[i + 1]);
  if (inv) {
    for (i = 0; i < tr_len; i++)
      tr[i] = -tr[i];
  }
  return 0;
}

int tr_len_nd(unsigned short int *e) {
  return ((e[0] * 2) & 0x3ff) - 16;
}

int evt_len_nd(unsigned short int * e) {
  return  2 * (e[0] & 0x3ff) + 2;
}
