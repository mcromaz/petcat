/*** misc util routines */

#include <assert.h>
#include "petcat.h"

int segevt_id(uint2 *e) {		/* from segevt.c */

  assert(e != 0);				/* pre */
  assert(e[0] == 0xaaaa && e[1] == 0xaaaa);
  return (get.vec(e + 2) << 4) + get.ch(e + 2);
}

int segevt_ener(uint2 *e, int shift) { /* from segevt.c */

  int rawE, corrE;

  assert(e != 0);				/* pre */
  assert(e[0] == 0xaaaa && e[1] == 0xaaaa);

  int i=0;
  
  return get.ener(e + 2) >> shift;
}

void adjoff(int *tr, int tr_len) { /* from eb_utils.c */

  int s = 0, i;

  assert(tr != 0);		/* pre */
  assert(tr_len > 25);

  for (i = 0; i < 25; i++) {
    s += tr[i];
  } 
  if (s >= 0)
    s =  (s + 7) / 25;
  else
    s =  (s - 7) / 25;
  //printf("tr_len = %d, s = %d\n", tr_len, s);
  for (i = 0; i < tr_len; i++) {
    tr[i] -= s;
  }
}

int net(int *tr, int tr_len) { 		/* from eb_utils.c */

  int avg = 0, i;

  assert(tr != 0);	/* pre */
  assert(tr_len > 0);

  for (i = tr_len - 10; i < tr_len; i++) {
    avg += tr[i];
    // printf("%d %d\n", i, tr[i]);
  }

  //printf("avg %d \n\n\n\n", avg);
  
  return (avg > 3*80) ? 1 : 0;
}
