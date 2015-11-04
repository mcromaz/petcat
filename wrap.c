#include <stdlib.h>
#include <assert.h>
#include "wrap.h"

void *Malloc(size_t size) {
  void *p = malloc(size);
  assert(p != 0);
  return p;
}

void *Calloc(size_t nobj, size_t size) {
  void *p = calloc(nobj, size);
  assert(p != 0);
  return p;
}
