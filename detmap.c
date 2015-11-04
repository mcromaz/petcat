/*** detmap.c: routines for using the detector map file */

#include <assert.h>
#include <stdio.h>
#include "eb.h"
#include "wrap.h"

pair *read_mapfile(char *filename, float *a0, float *a1) {

  FILE *fin;
  pair *map;
  int vec, chan, xtal, slot, num, i, j;
  float u, v;
  long long int slot_pattern = 0ll;
  int cc_vec;

  assert(filename != 0);	/* pre */
  fin = fopen(filename, "r");
  if (fin == 0)
    return 0;

  map = Malloc(4096 * sizeof(pair));
  for (i = 0; i < 4096; i++) {
    map[i].a = -1;
    map[i].b = -1;
  }
  while ((num = fscanf(fin, "%x %d %d %d %f %f\n", &vec, &chan, &xtal, &slot, &u, &v)) >= 0) { 
    assert(num == 6);
    assert(vec > 0 && vec <= 256 && chan >= 0 && chan < 10);
    assert(xtal >= 0 && xtal < EBUF_SIZE);
    assert(slot >= 0 && slot < 37);
    map[(vec << 4) + chan].a = slot;
    map[(vec << 4) + chan].b = xtal;
    if (a0 != 0) {
      a0[slot] = u;
    }
    if (a1 != 0) {
      a1[slot] = v;
    }
    slot_pattern |= (1ll << slot);
    if (slot == 36) {
      cc_vec = vec;
    }
      
  }
  assert(slot_pattern == 0x1fffffffff);		/* all 37 slots read */

  /* add entries for slots 37, 38, 39 (unused cc energies) in map file */
  j = 0;
  for (i = 0; i < 4; i++) {
    vec = i + 0x3;
    if (vec != cc_vec) {                        /* do not include cc in detmap file, it goes in slot 36 */
      map[(vec << 4) + 9].a = 37 + j;
      j++;
    }
  }

  return map;
}
