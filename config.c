/*** config.c - routines which support config of preproc/decomp */

/*
  read_mapfile() from detmap.c
  read_filterfile() from filter.c
  read_param() from eb_util.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include "petcat.h"
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

struct filter *read_filterfile(char *filename) {

  FILE *fin;
  struct filter *fltr;
  char s[80];
  int num;

  assert(filename != 0);	/* pre */

  fin = fopen(filename, "r");
  if (fin == 0)
    return 0;

  fltr = Malloc(sizeof(struct filter));

  num = fscanf(fin, "%s %f %f\n", s, &(fltr->emin), &(fltr->emax));
  assert(num == 3);
  assert(strcmp(s, "egate") == 0);
  assert(fltr->emin >= 0 && fltr->emax < 16384 && fltr->emin <= fltr->emax);

  num = fscanf(fin, "%s %d %d\n", s, &(fltr->xmin), &(fltr->xmax));
  assert(num == 3);
  assert(strcmp(s, "xgate") == 0);
  assert(fltr->xmin >= 0 && fltr->xmax < 16384 && fltr->xmin <= fltr->xmax);

  num = fscanf(fin, "%s %d %d\n", s, &(fltr->segmin), &(fltr->segmax)); 
  assert(num == 3);
  assert(strcmp(s, "seg") == 0);
  assert(fltr->segmin >= 0 && fltr->segmax < 16384 && fltr->segmin <= fltr->segmax);

  fclose(fin);
  return fltr;
}

int read_param(const char *filename, const char *label, float *x, int len) {

  FILE *fin;
  char s[MAX_LINE_LEN], *tok0, *tok, *u, *s0;
  int cnt = 0;
  float a;

  assert(filename != 0);		/* pre */
  assert(label != 0);
  assert(x != 0);
  assert(len >= 0);

  fin = fopen(filename, "r");
  if (fin == 0)
    return -1;

  while (fgets(s, MAX_LINE_LEN, fin) != 0) {
    tok0 = strtok(s, " ");
    if ((u = index(tok0, ':')) != 0) {
      *u =  '\0';
      if (strcmp(label, tok0) == 0) {
        while (1) {
          s0 = fgets(s, MAX_LINE_LEN, fin);
          if (s0 == 0 || s[0] == '\n') {
            fclose(fin);
            return cnt;
          }
          s0 = s;
          while ((tok = strtok(s0, " \t\n")) != 0) {
            s0 = 0;
            a = (float) atof(tok);
            x[cnt++] = a;
            if (cnt >= len) {
              fclose(fin);
              return cnt;
            }
          }
        }
      }
    }
  } 
  fclose(fin);
  return -1;
}
