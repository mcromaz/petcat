/*** log3: routines supporting diagnostics for mode3 data */

#include <stdio.h>
#include <assert.h>
#include "mode3io.h"

/* log global header information */

static FILE *flog;

void gh_log3(struct gebData *hdr) {

  static int init = 0;

  if (!init) {
    flog = fopen("log", "w");
    assert(flog != 0);
    init = 1;
  }
  else {
    fprintf(flog, "hdr %d %d %lld\n", hdr->type, hdr->length, hdr->timestamp);
  }
}
