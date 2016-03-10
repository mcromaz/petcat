/* program to decompose event signals for GRETINA / GRETA.

   Author:  D.C. Radford    Aug 2004
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "petcat.h"
#include "wrap.h"

/*
Include this if you need to log for debugging, but it causes a circular build
dependency, since it comes from gretClust.
#include <logSend.h>
*/
int GRID_PTS;   /* actual number of grid points in the basis */
int GRID_SEGS;  /* actual number of signals calculated for each basis point */
int TIME_STEPS;  /*  actual time steps calculated/measured for each segment */
int SRAD;                 /* actual range of r in basis grid */
int SPHI;                 /* actual range of phi in basis grid */
int SZZZ;                 /* actual range of z in basis grid */

long int lrintf(float x);
#define DECOMP_PATHLEN 100

static int dump_mat = 0;
char mat_file_name[DECOMP_PATHLEN];
FILE *mat_file, *pos_file;

/* Penalty factor for 2 interactions compared to 1 interaction (net=1), was 1 */
#define PENALTY1  ((best == 6) ? 0.97f : 1.0f)
/* Penalty factor for 3 interactions compared to 2 interactions (net=1), was 0.98 */
#define PENALTY2  0.97f
/* Penalty factors for 2 interactions compared to 1 interaction (net=2), was 0.98 for both */
#define PENALTY3  0.97f
#define PENALTY4  ((best == jfit-1) ? 1.0f : PENALTY3)
/*
#define PENALTY1  0.0f//((best == 6) ? 0.97f : 1.0f)
#define PENALTY2  0.0f//0.97f
#define PENALTY3  0.0f//0.97f
#define PENALTY4  0.0f//((best == jfit-1) ? 1.0f : PENALTY3)
*/

Basis_Point   *basis;                 /* basis-signal data */
int           grid_pos_lu[SSEG][MAX_SRAD][MAX_SPHI][MAX_SZZZ];    /* basis-grid position look-up table */
int           maxir[SSEG], maxip[SSEG], maxiz[SSEG];  /* max. values of ir, ip, iz for each segment */
Adaptive_Grid ags1[SSEG];             /* Adaptive Grid Search coarse grid, 1 segment */

int           quiet;                  /* set to 1 to prevent output of diagnostic info */
int           average_sigs;           /* set to 1 to prevent output averaged observed and 
                                         fitted signals for single-segment (net=1) events */
int           *bad_segs;              /* list of segment numbers that should be disabled */

int open_mat_file() {
   char namebuf[DECOMP_PATHLEN + 5];

   if (mat_file_name[0] == 0) return -1;
   snprintf(namebuf, DECOMP_PATHLEN + 4, "%s.mat", mat_file_name);
   if (!(mat_file=fopen(namebuf, "w+"))) {
      printf("Cannot open %s\n", namebuf);
      return -1;
   }
   snprintf(namebuf, DECOMP_PATHLEN + 4, "%s.txt", mat_file_name);
   if (!(pos_file=fopen(namebuf, "w+"))) {
      printf("Cannot open %s\n", namebuf);
      return -1;
   }
   dump_mat = 1;
   return 0;
}
int close_mat_file() {
   dump_mat = 0;
   fclose(mat_file);
   fclose(pos_file);
   return 0;
}

int write_mat_file(struct decomp_state *di) {
   int retval = 0;
   char *s;
   if (dump_mat) {
      retval = fwrite(di->mat, 8192, 1, mat_file);
      if (retval != 1) { 
         printf("mat binary file write returned %d\n", retval);
         return -1;
      }
      s = dl_crys_intpts_2s(&di->pos);
      retval = fprintf(pos_file, "%s", s);
      free(s);
      if (retval < 0) {
         printf("mat text file write returned %d\n", retval);
         return -1;
      }
   }
   return 0;
}

int dl_decomp_init_global(char *basis_file, int q) {

  assert(basis_file != 0);
  
  if (read_basis(basis_file))
     return 1;
  if (grid_init())
     return 2;

  quiet = q;

  return 0;
}

struct decomp_state *dl_decomp_init_thread() {

  struct decomp_state *inst;

  inst = Calloc(1, sizeof(struct decomp_state));
  inst->bsig = Calloc(1, sizeof(Event_Signal));
  inst->ints = Calloc(2 * MAX_SEGS, sizeof(Interaction));
  inst->err = Calloc(1, sizeof(struct decomp_errcnt));
  inst->coal_dist = COAL_DIST_DEFAULT;
/*
  sendLogMsg("decompLib", "decomp_init_thread");
  sendToLog(0);
*/


  return inst;
}

/* this form of dl_decomp_init is retained for compatibility */
struct decomp_state *dl_decomp_init(char *basis_file_name, int set_quiet) {

   if (!dl_decomp_init_global(basis_file_name, set_quiet)) {
      return dl_decomp_init_thread();
   }
   return 0;
}

/* --------------------- */

int postprocess_events(Interaction *ints, int nints, float total_e,
		       int ppflag, float coal_dist,
		       double *x, double *y, double *z, double *e, int *s, postprocCnt *postCnt)
{
  /* prost-process the events
     convert cyl to cart coords
     combine closely-spaced interactions in neighboring segments
     optionally combine closely-spaced interactions in the same segment
     Input: ints, nints, total_e
     Output: x, y, z, e
     Control: ppflag = 2 to combine interactions in the same segment
     returns final number of interactions found
  */

  double efrac, elo, dinp, dinp12, dinp13, dinp23, e1, e2, dx, dy, dz;
  int    i, j, k, found, i12, i13, i23, seg[MAX_PARS];

  /* convert cyl to cart coords */
  for (i=0; i<nints; i++) {
    cyl_to_cart(ints[i].seg, &ints[i].r, &x[i], &y[i], &z[i], &e[i]);
    seg[i] = ints[i].seg;
    s[i] = seg[i];			/* mc */
    e[i] *= total_e;
  }
  found = nints;    /* will be final number of interactions found */

  int seghit[36] = {0};

  for (i=0; i<found-1; i++) {
    seghit[seg[i]]++;
  }
  for (i=0; i<36; i++) {
    if (seghit[i] == 1) { postCnt->oneInt[i]++; }
    else if (seghit[i] == 2) { postCnt->twoInt[i]++; }
    else if (seghit[i] == 3) { postCnt->threeInt[i]++; }
    else if (seghit[i] > 3) { postCnt->manyInt[i]++; }
  }

  /* check for closely-spaced interactions in neighboring segments
     that should be combined. Needs to be done in Cartesian coords */
  for (i=0; i<found-1; i++) {
    for (j=i+1; j<found; j++) {
      if (seg[i] == seg[j]) continue;
      dx = x[i] - x[j];
      dy = y[i] - y[j];
      dz = z[i] - z[j];
      if ((dx*dx + dy*dy + dz*dz) < coal_dist*coal_dist) {
	/* FIXME - We could add a fudge factor to coal_dist if desired */
	e1 = e[i];
	e2 = e[j];
	x[i] = (e1*x[i] + e2*x[j]) / (e1+e2);
	y[i] = (e1*y[i] + e2*y[j]) / (e1+e2);
	z[i] = (e1*z[i] + e2*z[j]) / (e1+e2);
	e[i] = e1+e2;
	postCnt->combinedDiffSeg++;
	found--;
	for (k=j; k<found; k++) {
	  x[k] = x[k+1];
	  y[k] = y[k+1];
	  z[k] = z[k+1];
	  e[k] = e[k+1];
	  seg[k] = seg[k+1];
	}
      } 
    }
  }

  /* This should be postprocessing*/
 /* if (found == 2 && seg[0] == seg[1] && ((e[0]<80) || (e[1]<80) && (fabs(e[0])>1 && fabs(e[1])>1))){ 
	e1 = e[0];
	e2 = e[1];
	x[0] = (e1*x[0] + e2*x[1]) / (e1+e2);
	y[0] = (e1*y[0] + e2*y[1]) / (e1+e2);
	z[0] = (e1*z[0] + e2*z[1]) / (e1+e2);
	e[0] = e1+e2;
	found--;
 } */

  if (ppflag < 2) return found;  /* return without doing the same-segment
				    distance- or energy-based coalescence...
				    this should be the default? */

  /* now combine closely-spaced events in the same segment */
  /* FIXME - this should really be done in grid coordinates, if at all
     But that's not easy since the code above needs to be done in
     cart coords, and can change the number of interactions */

  /* this was in the original lbl code -
     I'm not sure how they got it, or if it's really a good idea */
  /* FIXME - It should at least include some dependence on number of hit segments */
  if (total_e < 140.0) {
    elo = 0.35;
  } else if (total_e < 400.0) {
    elo = 0.28;
  } else if (total_e < 900.0) {
    elo = 0.20;
  } else {
    elo = 0.15;
  }

  if (found == 3 && seg[0] == seg[1] && seg[0] == seg[2]) {
    /* 3 interactions in a one-seg event */
    /* difference in position */
    dinp12 = sqrt((x[0]-x[1])*(x[0]-x[1]) + (y[0]-y[1])*(y[0]-y[1]) + (z[0]-z[1])*(z[0]-z[1]));
    dinp13 = sqrt((x[0]-x[2])*(x[0]-x[2]) + (y[0]-y[2])*(y[0]-y[2]) + (z[0]-z[2])*(z[0]-z[2]));
    dinp23 = sqrt((x[1]-x[2])*(x[1]-x[2]) + (y[1]-y[2])*(y[1]-y[2]) + (z[1]-z[2])*(z[1]-z[2]));

    if (e[0] < elo*0.7     || e[1] < elo*0.7     || e[2] < elo*0.7    ||
	dinp12 < coal_dist || dinp13 < coal_dist || dinp23 < coal_dist) {
      /* combine 2 of the 3 interactions */
      i12 = i13 = i23 = 0;
      if (dinp12 < coal_dist && dinp12 < dinp13 && dinp12 < dinp23) {
	/* separation less than 2 mm; combine the interactions */
	i12 = 1;
      } else if (dinp13 < coal_dist && dinp13 < dinp12 && dinp13 < dinp23) {
	/* separation less than 2 mm; combine the interactions */
	i13 = 1;
      } else if (dinp23 < coal_dist && dinp23 < dinp12 && dinp23 < dinp13) {
	/* separation less than 2 mm; combine the interactions */
	i23 = 1;
      } else if (e[0] < elo*0.7 && e[0] < e[1] && e[0] < e[2]) {
	/* first interaction is too weak; combine it with the nearest other one */
	if (dinp12 < dinp13) {
	  i12 = 1;
	} else {
	  i13 = 1;
	}
      } else if (e[1] < elo*0.7 && e[1] < e[0] && e[1] < e[2]) {
	/* second interaction is too weak; combine it with the nearest other one */
	if (dinp12 < dinp23) {
	  i12 = 1;
	} else {
	  i23 = 1;
	}
      } else if (e[2] < elo*0.7 && e[2] < e[0] && e[2] < e[1]) {
	/* third interaction is too weak; combine it with the nearest other one */
	if (dinp13 < dinp23) {
	  i13 = 1;
	} else {
	  i23 = 1;
	}
      }
      /* do the coalescence */
      if (i12) {
	i=0; j=1;
      } else if (i13) {
	i=0; j=2;
      } else {
	i=1; j=2;
      }
      e1 = e[i];
      e2 = e[j];
      x[i] = (e1*x[i] + e2*x[j]) / (e1+e2);
      y[i] = (e1*y[i] + e2*y[j]) / (e1+e2);
      z[i] = (e1*z[i] + e2*z[j]) / (e1+e2);
      e[i] = e1+e2;
      found--;
      postCnt->combinedWithinSeg++;
      postCnt->combined3to2[seg[0]]++;
      for (k=j; k<found; k++) {
	x[k] = x[k+1];
	y[k] = y[k+1];
	z[k] = z[k+1];
	e[k] = e[k+1];
	seg[k] = seg[k+1];
      }

      /* see if the remaining two interactions should also be combined */
      efrac = e1/(e1+e2);
      dinp  = sqrt((x[0]-x[1])*(x[0]-x[1]) + (y[0]-y[1])*(y[0]-y[1]) + (z[0]-z[1])*(z[0]-z[1]));
      if (dinp < coal_dist || efrac > (1.0 - elo) || efrac < elo) {
	/* separation less than 2 mm, or one interaction is too weak */
	/* combine the interactions */
	e1 = e[0];
	e2 = e[1];
	x[0] = (e1*x[0] + e2*x[1]) / (e1+e2);
	y[0] = (e1*y[0] + e2*y[1]) / (e1+e2);
	z[0] = (e1*z[0] + e2*z[1]) / (e1+e2);
	e[0] = e1+e2;
	found--;
	postCnt->combinedWithinSeg++;
	postCnt->combined3to1[seg[0]]++;
	postCnt->combined3to2[seg[0]]--; /* we've combined all 3, but I incremented 3to2 counter earlier */
      }
    }

  } else {    /* more than 1 segment hit, or only 2 interactions... */

    /* loop over the hit segments */
    for (i=0; i<found-1; i++) {
      if (seg[i] == seg[i+1]) {
	/* two interactions in this segment; see if they should be combined */
	j = i+1;
	e1 = e[i];
	e2 = e[j];
	/* difference in position */
	dinp  = sqrt((x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]));
	/* energy fraction for first interaction */
	efrac = e1/(e1+e2);
	if (dinp < coal_dist || efrac > (1.0 - elo) || efrac < elo) {
	  /* separation less than 2 mm, or one interaction is too weak */
	  /* do the coalescence */
	  x[i] = (e1*x[i] + e2*x[j]) / (e1+e2);
	  y[i] = (e1*y[i] + e2*y[j]) / (e1+e2);
	  z[i] = (e1*z[i] + e2*z[j]) / (e1+e2);
	  e[i] = e1+e2;
	  found--;
	  postCnt->combinedWithinSeg++;
	  postCnt->combined2to1[seg[i]]++;
	  for (k=j; k<found; k++) {
	    x[k] = x[k+1];
	    y[k] = y[k+1];
	    z[k] = z[k+1];
	    e[k] = e[k+1];
	    seg[k] = seg[k+1];
	  }
	}
      } /* seg[i] == seg[i+1] */
    } /* loop over found */
  }
  return found;
} /* postprocess_event */

/* --------------------- */

/* routine to actually do the decomposition (by calling decompose_1 or decompose_n) */

struct crys_intpts *dl_decomp(struct decomp_state *di, Event_Signal *asig, postprocCnt *postCnt) {

  /*** dl_decomp modified to return a crystal event regardless of its ability to decompose.
       If an event cannot be decomposed an empty crys_intpts struct is allocated with on its type and pad
       fileds set.  The pad field is 0 if event met standard decomp criteria or set to the following values if not:

       pad = 1   a null pointer was passed to dl_decomp()
           = 2   total energy below threshold
           = 3   no net charge segments in evt
           = 4   too many net charge segments
           = 5   chi^2 is bad following decomp (in this case crys_intpts is non-zero but post-processing step is not applied)
  */

  struct crys_intpts *ci;
  double xfinal[MAX_PARS], yfinal[MAX_PARS], zfinal[MAX_PARS], efinal[MAX_PARS];
	int sfinal[MAX_PARS];	/* keep track of segments hit */
  double esum, chisq, chisq2, t0;
  int    nseg, nints, found, i, ii, t, seg[TOT_SEGS];
  Basis_Point *b;
  simpleInteraction ints[MAX_SEGS];	/* for single interaction search */
  int stat;

  ci = (struct crys_intpts*) Calloc(1, sizeof(struct crys_intpts));
  /* ci->type = 0xabcd1234; */
  ci->type = 0xabcd5678;
  
  /* assert(asig != 0);	 */
  
  if (!asig) {
    printf("NULL asig\n");
    ci->pad = 1; /* Null pointer */
    return ci;
  }
  
  ci->timestamp = asig->time;
  for (i = 0; i < 4; i++) {
    ci->core_e[i] = asig->core_e[i];
  }
  ci->tot_e = asig->total_energy;
  ci->prestep = asig->prestep;
  ci->poststep = asig->poststep;

  if (asig->total_energy < 20.0) {
    printf("low total energy\n");
    ci->pad = 2; /* Total energy below threshold */
    return ci;
  }
  
  /* set bad segment signals to zero */
  for (i=0; bad_segs[i]>=0; i++) {
    ii = bad_segs[i];
    for (t=0; t<TIME_STEPS; t++) asig->signal[ii][t] = 0.0f;
    asig->seg_energy[ii] = 0.0f;
  }
  
  /* determine which segments have net energy */
  nseg = 0;
  esum = 0.0;
  for (i = 0; i < TOT_SEGS; i++) {
    seg[i] = 0;
    if (asig->seg_energy[i] > 30.0)  {
      esum += asig->seg_energy[i];
      seg[nseg++] = i;
    }
  }
  
  if (nseg == 0) {
    di->err->nonet++;
    printf("no net\n");
    ci->pad = 3; /* No net charge segments */
    return ci;
  }
  if (nseg > MAX_SEGS) {
    di->err->toomanynet++;
    printf("too many net - nseg = %d\n", nseg);
    ci->pad = 4; /* Too many net charge segments */
    for (i = 0; i < MAX_INTPTS; i++) {
      ci->intpts[i].seg = i;
      ci->intpts[i].seg_ener = asig->seg_energy[i];
    }
    return ci;
  }
  
  /* remove for high-rate runs with cc pileup */
  // if (fabs(esum - asig->total_energy) > 90.0) {
  //   di->err->sumener++;
  //   printf("bad sum\n");
  //   return 0;  
  // }
 
  (void) memset(di->bsig, 0, sizeof(Event_Signal));
  (void) memset(di->ints, 0, 2 * MAX_SEGS * sizeof(Interaction));

#ifdef ONEINT
  //printf("nseg=%d\n", nseg);
  if (nseg >= 1 && nseg <= 4) {
    stat = normDiag(asig, nseg, seg);
    if (badTrPred(asig, nseg, seg)) {
      fprintf(stderr, "bad net segment trace\n");
      ci->pad = 40;	/* oneInt search failed - bad net trace */
      return ci;
    }
    //stat = oneIntSearchIter(asig, di->bsig, nseg, seg, ints, &chisq);
    stat = oneIntSearch2(asig, di->bsig, nseg, seg, ints, &chisq);
    if (stat != 0) {
      ci->pad = 41;	/* oneInt search failed */
      return ci;
    }

    /* employ no post-processng, output ci directly */
    ci->num = nseg;
    ci->tot_e = asig->total_energy;
    for (i = 0, ci->baseline = 0.0; i < 20; i++) {  /* calc baseline */
      ci->baseline += (float) asig->signal[36][i];
    }
    ci->baseline /= 20.0;
    ci->chisq = chisq;
    ci->norm_chisq = 0;
    ci->timestamp = asig->time;
    assert(b != 0);
    for (i = 0; i < nseg; i++) {
      ci->intpts[i].x = ints[i].x;
      ci->intpts[i].y = ints[i].y;
      ci->intpts[i].z = ints[i].z;
      ci->intpts[i].e = asig->seg_energy[seg[i]];
      ci->intpts[i].seg = seg[i];
      ci->intpts[i].seg_ener = asig->seg_energy[seg[i]];
    }
    ci->pad = 42;	/* properly decomposed oneInt event */
    return ci; 
  }
  ci->pad = 43;		/* numNet out of range */ 
  return ci;
#else
  /* branch according to number of hit segments */
  if (nseg == 1) {
    //nints = decompose_1(asig, di->bsig, seg[0], di->ints, &t0, &chisq, 0, 1, 1, 1, 0, 0, 1, 1, 0.1);
    nints = decompose_1(asig, di->bsig, seg[0], di->ints, &t0, &chisq, 0, 1, 1, 1, 1, 1, 1, 1, 0.1);
  }
  else {
    printf("nseg >1\n");
    nints = decompose_n(asig, di->bsig, nseg, seg, 1, di->ints, &t0, &chisq);
  }
#endif

  /* no need to have such cut here */ 
  if (chisq > 99.9) {
    di->err->badchisq++;
    ci->pad = 5; /* Chi squared is bad */
    // return ci;
  }

  di->cnt++;

  /* post_process the results
     ppflag = 2; --  combine interactions in the same seg based on 
                     energy or distance
     ppflag = 1; -- do not combine interactions in the same seg 
                    based on energy or distance */
  found = postprocess_events(di->ints, nints, asig->total_energy,
			     2, di->coal_dist, xfinal, yfinal, zfinal, 
			     efinal, sfinal, postCnt);
  chisq2 = chisq / (float) (MEAS_SEGS * TIME_STEPS);
  chisq2 /= (4.0 / asig->total_energy) * (4.0 / asig->total_energy);
  
  /* crytal evt always returned, alloc crys_intpt at beginning of routine */
  // ci = (struct crys_intpts*) Calloc(1, sizeof(struct crys_intpts));
  // ci->type = 0xabcd1234;
  ci->num = found;
  ci->tot_e = asig->total_energy;

  /* calc baseline */
  for (i = 0, ci->baseline = 0.0; i < 20; i++) {
    ci->baseline += (float) asig->signal[36][i];
  }
  ci->baseline /= 20.0;
  
#ifdef TIMEOFFSET
  /* DCR: add in saved offset coming from time alignment 
     of signals in preproc */
  t0 += asig->time_offset;
  ci->intpts[MAX_INTPTS-1].z = asig->time_offset;  /* DCR UGLY hack; saves offset for comparing time spectra resolutions downstream; this line can be removed without changing proper functionality! */
#endif
  ci->t0 = t0;
  ci->chisq = chisq;
  ci->norm_chisq = chisq2;
#ifdef TIMESTAMP
  ci->timestamp = asig->time;
#endif
  for (i = 0; i < found; i++) {
    ci->intpts[i].x = xfinal[i];
    ci->intpts[i].y = yfinal[i];
    ci->intpts[i].z = zfinal[i];
    ci->intpts[i].e = efinal[i];
    ci->intpts[i].seg = sfinal[i];
    ci->intpts[i].seg_ener = asig->seg_energy[sfinal[i]];
  }
  if (dump_mat) {
    int j;
    for(i=0; i<4096; i++) {
      di->mat[i] = 0;
    }
    for (i=0; i<MEAS_SEGS; i++) {
      for (j=0; j<TIME_STEPS; j++) {
	di->mat[50*i + j] = lrintf(asig->signal[i][j] * 10000.0f);
	di->mat[2000 + 50*i + j] = lrintf(di->bsig->signal[i][j] * 10000.0f);
      }
    }
    di->pos = *ci;
  } 
  
  return ci;
} /* dl_decomp */

/* --------------------- */

/* converts an interaction-point structure into a string
   in the format of the old original gdecomp program */
char *dl_crys_intpts_2s(struct crys_intpts *x) {

  char *s;
  char int_fmt[] = "%11.2f %11.2f %11.2f %11.2f\n";
  int  hdr_line_len;
#ifdef TIMESTAMP
  char hdr_fmt_1[] = "%2d %2d %9.2f %6.2f %24.6f %9.4f %24lld\n";
  int  hdr_line_len_1 = 83;
#else
  char hdr_fmt_0[] = "%2d %2d %9.2f %6.2f %24.6f %9.4f %4d\n";
  int  hdr_line_len_0 = 63;
#endif
  int  int_line_len   = 48;
  int  out_line_len, num, i;

#ifdef TIMESTAMP
  hdr_line_len = hdr_line_len_1;
#else
  hdr_line_len = hdr_line_len_0;
#endif
  out_line_len = hdr_line_len + x->num * int_line_len + 1;
  s = (char*) Calloc(1, out_line_len);
#ifdef TIMESTAMP
  num = sprintf(s, hdr_fmt_1,
		x->num, -1, x->tot_e, x->t0, x->chisq, x->norm_chisq, x->timestamp);
#else
  num = sprintf(s, hdr_fmt_0,
		x->num, -1, x->tot_e, x->t0, x->chisq, x->norm_chisq, -1);
#endif

  /* assert(num == hdr_line_len);*/
  for (i = 0; i < x->num; i++) {
    num = sprintf(s + hdr_line_len + i * int_line_len, int_fmt,
		  x->intpts[i].x, x->intpts[i].y, x->intpts[i].z, x->intpts[i].e);
    /*  assert(num == int_line_len);*/
  }
  return s;
} /* dl_crys_intpts_2s */

/* --------------------- */

/* sets coalescence-length pamaeter, normally 2 mm by default */
void dl_set_coal_dist(struct decomp_state *inst, float d) {

  assert(d >= 0.0);
  inst->coal_dist = d;

  return;
} /* dl_set_coal_dist */

/* --------------------- */

/* extracts the net number of hit segments from the measured signal */
int num_net(Event_Signal *asig) {

  int nseg = 0, i;
  for (i = 0; i < 36; i++) {
    if (asig->seg_energy[i] > 30.0) {
      nseg++;
    }
  }
  return nseg;
} /* num_net */

/* --------------------- */

int decompose_1(const Event_Signal *asig,  /* observed signals */
		Event_Signal *bsig, /* fitted signals */
		int seg, Interaction *ints, double *t0, double *chisq_out, /* parameters */
		int grid2, int fit0, int fit1, int fit2, int fit3,         /* control */
		int fit4, int final_fit, int coalesce, double min_e_fraction)  /* control */
{
	
  /*
    performs the adaptive grid search algorithm for single-segment events,
    optionally followed by constrained least-squares
    
    returns number of best-fit interactions found
    ints[] array is set to contain the resulting best-fit double-interaction
    grid2, fit0, fit1, fit2, fit3, fit4, final_fit and coalesce control what parts
    of the decomposition algorithm get run
  */
  
  double chisq0, chisq2, chisq_cg;
  int    i, j, k, m, n, best = 0, best2, shift_t0 = 0, nints = 0;
  float  e1, e2;
  Event_Signal  asig_temp;
  Interaction   f_ints[11][4];     /* test/fitted interactions */
  Event_Signal  f_bsig[11];        /* calculated/fitted event data */
  double        f_t0[11], chisq[11];

  *t0 = 0.0;
  ints[2].seg = seg;
  ints[2].e   = 0.0;
  chisq2 = 100.0;
  for (i=0; i<7; i++) {
    chisq[i] = 100.0;
    f_t0[i] = 0.0;
  }
  if (fit0) grid2 = 0;
  
  memset(f_ints, 0, 7 * 4 * sizeof(Interaction));
  
  /* start off with a simple 1-interaction fit */
  if (fit0) {
    f_ints[6][0].seg = seg;
    f_ints[6][0].pos = -1;
    f_ints[6][0].r   = maxir[seg]/2;
    f_ints[6][0].p   = maxip[seg]/2;
    f_ints[6][0].z   = maxiz[seg]/2;
    f_ints[6][0].e   = 1.0;
    f_ints[6][1].seg = seg;
    f_ints[6][1].e   = 0.0;
    f_t0[6] = 3;   /* NOTE Changed, was f_t0[6] = 0 */
    chisq[6] = fitter(asig, &f_bsig[6], f_ints[6], &f_t0[6], 1, 0);
    best = 6;

    /* shift the remaining signal to bring t0 close to zero */
    asig_temp = *asig;
    shift_t0 = f_t0[best] - 0.0;
    if (shift_t0 > 0) {
      if (!quiet) printf("*----* shifting t0 by %d\n", shift_t0);
      for (i=0; i<MEAS_SEGS; i++) {
	for (j=0; j<TIME_STEPS-shift_t0; j++) {
	  asig_temp.signal[i][j] = asig_temp.signal[i][j+shift_t0];
	}
	for (j=TIME_STEPS-shift_t0; j<TIME_STEPS-1; j++) {
	  asig_temp.signal[i][j] = asig_temp.signal[i][TIME_STEPS-1];
	}
      }
    } else {
      shift_t0 = 0;
    }
    
    /* do the grid search over the course grid for the hit segment */
    chisq_cg = coarse_grid_1(&asig_temp, seg, f_ints[0], &chisq0, min_e_fraction);
    f_t0[0] = shift_t0;
    chisq[0] = eval_int_pos(asig, bsig, f_ints[0], f_t0[0], 2);
    if (chisq[0] < PENALTY1 * chisq[best]) best = 0;
    
  } else {

    /* do the grid search over the course grid for the hit segment */
    /* NOTE could guess that t0 is 2 or 3, rather than zero, at this point */
    chisq_cg = coarse_grid_1(asig, seg, f_ints[0], &chisq0, min_e_fraction);
    chisq[0] = eval_int_pos(asig, bsig, f_ints[0], f_t0[0], 2);
  }
  if (!quiet) {
    printf("** grid1: seg, chisq, pars: %2d, %7.4f, %7.4f", seg, chisq_cg, chisq[0]);
    for (i=0; i<2; i++) printf(", %6.3f %6.3f %6.3f %6.3f",
			       f_ints[0][i].r, f_ints[0][i].p,
			       f_ints[0][i].z, f_ints[0][i].e);
    printf("\n");
  }
  
  if (fit1) {
    /* do nonlinear least-squares fit
       starting from coarse grid point pair with lowest chisq */
    f_ints[1][0] = f_ints[0][0];
    f_ints[1][1] = f_ints[0][1];
    f_t0[1] = f_t0[0] - 1.0;
    if (f_t0[1] < 0.0) f_t0[1] = 0.0;
    chisq[1] = fitter(asig, &f_bsig[1], f_ints[1], &f_t0[1], 2, 0);
    if (chisq[1] < PENALTY1 * chisq[best]) best = 1;
  }

  if (fit2 && f_ints[0][2].pos >= 0 && f_ints[0][3].pos >= 0) {
    /* try nonlinear least-squares fit
       starting from coarse grid point pair with second best chisq */
    f_ints[2][0] = f_ints[0][2];
    f_ints[2][1] = f_ints[0][3];
    f_t0[2] = f_t0[0] - 1.0;
    if (f_t0[2] < 0.0) f_t0[2] = 0.0;
    chisq[2] = fitter(asig, &f_bsig[2], f_ints[2], &f_t0[2], 2, 0);
    if (chisq[2] < PENALTY1 * chisq[best]) best = 2;
  }

  if (grid2) {
    /* now do adaptive part of AGS,
       i.e. check neighboring sites on fine grid */
    chisq2 = refine_grid_1(asig, chisq_cg, chisq0, min_e_fraction, f_ints[0]);
    if (chisq2 < chisq_cg - 0.00000001) {
      chisq[0] = eval_int_pos(asig, bsig, f_ints[0], *t0, 2);
      if (!quiet) printf("** fine_grid1: seg, chisq: %2d, %7.4f, %7.4f\n", seg, chisq2, chisq[0]);
      if (chisq[0] < PENALTY1 * chisq[best]) best = 0;
      if (fit3) {
	/* try nonlinear least-squares fit
	   starting from new fine grid point pair */
	f_ints[3][0] = f_ints[0][0];
	f_ints[3][1] = f_ints[0][1];
	chisq[3] = fitter(asig, &f_bsig[3], f_ints[3], &f_t0[3], 2, 0);
	if (chisq[3] < PENALTY1 * chisq[best]) best = 3;
      }
    }
  }

  /*  if (fit4 && asig->seg_energy[seg] > 400.0 && chisq[best] > 0.028) {  NOTE hardwire limit */
    if (fit4 && asig->seg_energy[seg] > 150.0) {  /* NOTE hardwire limit */
    if (!quiet) printf("*** trying a fit with 3 interactions...\n");
    if (best == 6) {
      j = 0;
      if (chisq[1] < chisq[j]) j = 1;
      if (chisq[2] < chisq[j]) j = 2;
      if (chisq[3] < chisq[j]) j = 3;
    } else {
      j = best;
    }
    f_ints[4][0] = f_ints[j][0];
    f_ints[4][1] = f_ints[j][1];
    f_ints[4][2].seg = seg;
    f_ints[4][2].pos = -1;
    f_ints[4][2].r   = maxir[seg]/2;
    f_ints[4][2].p   = maxip[seg]/2;
    f_ints[4][2].z   = maxiz[seg]/2;
    f_ints[4][2].e   = asig->seg_energy[seg] / (3.0*asig->total_energy);
    f_t0[4] = f_t0[j];
    chisq[4] = fitter(asig, &f_bsig[4], f_ints[4], &f_t0[4], 3, 0);
    if (chisq[4] < PENALTY2 * chisq[best]) best = 4;
  }

  if (coalesce) {
    /* Try coalescing interactions, accepting/rejecting by chisq value */
    nints = 2;
    if (best == 4) {
      nints = 3;
    } else if (best == 6) {
      nints = 1;
    }
    best2 = best;
    if (nints == 3) {
      if (!quiet) printf("*** trying coalescence fits with 2 interactions...\n");
      for (j=7; j<10; j++) {
	if (j == 7) {
	  k=0; m=1; n=2;
	} else if (j == 8) {
	  k=1; m=0; n=2;
	} else {
	  k=2; m=1; n=0;
	}
	e1 = f_ints[best][m].e;
	e2 = f_ints[best][n].e;

	f_ints[j][0] = f_ints[best][k];
	f_ints[j][1] = f_ints[best][m];
	f_ints[j][1].r = (e1*f_ints[best][m].r + e2*f_ints[best][n].r) / (e1+e2);
	f_ints[j][1].p = (e1*f_ints[best][m].p + e2*f_ints[best][n].p) / (e1+e2);
	f_ints[j][1].z = (e1*f_ints[best][m].z + e2*f_ints[best][n].z) / (e1+e2);
	f_ints[j][1].e = e1+e2;
	f_t0[j] = f_t0[best];

	chisq[j] = fitter(asig, &f_bsig[j], f_ints[j], &f_t0[j], 2, 0);
	if (chisq[j] < chisq[best2]) best2 = j;
      }
      if (best2 != best && PENALTY2 * chisq[best2] <= chisq[best]) {
	best = best2;
	nints = 2;
      }
    }
	
    if (nints == 2) {
      if (!quiet) printf("*** trying coalescence fit with 1 interaction...\n");
      e1 = f_ints[best][0].e;
      e2 = f_ints[best][1].e;

      f_ints[10][0] = f_ints[best][0];
      f_ints[10][0].r = (e1*f_ints[best][0].r + e2*f_ints[best][1].r) / (e1+e2);
      f_ints[10][0].p = (e1*f_ints[best][0].p + e2*f_ints[best][1].p) / (e1+e2);
      f_ints[10][0].z = (e1*f_ints[best][0].z + e2*f_ints[best][1].z) / (e1+e2);
      f_ints[10][0].e = e1+e2;
      f_t0[10] = f_t0[best];

      chisq[10] = fitter(asig, &f_bsig[10], f_ints[10], &f_t0[10], 1, 0);

      if (PENALTY1 * chisq[10] <= chisq[best]) {
	best = 10;
	nints = 1;
      }
    }
  }

  /* decide which part of the algorithm gave the best chi-squared
     and, if necessary, finish up the nonlinear least-squares fit
     with more stringent convergence requirements */
  if (!quiet) printf(">>> g f1 f2 f3, best: %.4f %.4f %.4f %.4f %.4f %.4f, %.4f\n",
		     chisq[0], chisq[1], chisq[2], chisq[3],
		     chisq[4], chisq[6], chisq[best]);

  ints[0] = f_ints[best][0];
  ints[1] = f_ints[best][1];
  *t0 = f_t0[best];
  *chisq_out = chisq[best];
  nints = 2;
  if (best == 0) {
    if (!quiet) {  /* calculate bsig for matrix */
      chisq2 = eval_int_pos(asig, bsig, ints, *t0, 2);
      if (fabsf(chisq2 - chisq[best]) > 0.000001)
	printf("**** ACK!! best_chisq != eval_int_pos(); %f %f\n",
	       chisq[best], chisq2);
    }
  } else {
    if (best == 4) {
      nints = 3;
      ints[2] = f_ints[best][2];
    } else if (best == 6 || best == 10) {
      nints = 1;
    }
    if (final_fit) {
      *chisq_out = fitter(asig, bsig, ints, t0, nints, 1);
    } else {
      *bsig = f_bsig[best];
    }
  }

  return nints;
} /* decompose_1 */

/* ================================================================ */

int decompose_n(const Event_Signal *asig,  /* observed signals*/
		Event_Signal *bsig, /* fitted signals*/
		int nseg, int *seg, int coalesce, 
		Interaction *ints, double *t0, double *chisq_out)  /* parameters */
{
  /*
    performs the adaptive grid search algorithm for three-segment events,
    optionally followed by constrained least-squares

    returns number of best-fit interactions found
    ints[] array is set to contain the resulting best-fit multi-interaction
  */

  double chisq2;
  int    i, j, m, n, jseg, jfit, jint, segsave, best, shift_t0;
  float  e1, e2;
  Event_Signal  asig_temp;
  Interaction   f_ints[4*MAX_SEGS][3*MAX_SEGS];     /* test/fitted interactions */
  Event_Signal  f_bsig[4*MAX_SEGS];        /* calculated/fitted event data */
  double        f_t0[4*MAX_SEGS], chisq[4*MAX_SEGS];
  int           nints[4*MAX_SEGS];


  /* put segments in order of decreasing energy */
  for (i=1; i<nseg; i++) {
    for (j=i; j>0 && asig->seg_energy[seg[j-1]] < asig->seg_energy[seg[j]]; j--) {
      segsave  = seg[j-1];
      seg[j-1] = seg[j];
      seg[j]   = segsave;
    }
  }

  /* start off with a simple nseg-interaction fit, one in each hit segment */
  for (i=0; i<nseg; i++) {
    f_ints[0][i].seg = seg[i];
    f_ints[0][i].pos = -1;
    f_ints[0][i].r   = maxir[seg[i]]/2;
    f_ints[0][i].p   = maxip[seg[i]]/2;
    f_ints[0][i].z   = maxiz[seg[i]]/2;
    f_ints[0][i].e   = asig->seg_energy[seg[i]] / asig->total_energy;
  }
  nints[0] = nseg;
  f_t0[0] = 0;
  chisq[0] = fitter(asig, &f_bsig[0], f_ints[0], &f_t0[0], nseg, 0);
  best = 0;
  asig_temp = *asig;
  jfit = 1;

  /* loop throught the hit segments, trying to add interactions one at a time
     if the net energy gets too small, then quit  */
  for (jseg=0; jseg<nseg && asig->seg_energy[seg[jseg]] > 120.0; jseg++) {

    /* ----------------------------------------------
       now add a second interaction to the jseg-th segment;
       first subtract the calculated signal for the other segments
       from the observed signal */

    jint = nints[best] - nseg + jseg;
    for (i=0; i<jint; i++) {
      f_ints[jfit][i] = f_ints[best][i];
    }
    for (i=jint + 1; i<nints[best]; i++) {
      f_ints[jfit][i-1] = f_ints[best][i];
    }
    eval_int_pos(asig, bsig, f_ints[jfit], f_t0[best], nints[best]-1);
    for (i=0; i<MEAS_SEGS; i++) {
      for (j=0; j<TIME_STEPS; j++) {
	asig_temp.signal[i][j] -= bsig->signal[i][j];
      }
    }

    /* shift the remaining signal to bring t0 close to zero */
    shift_t0 = f_t0[best] - 1.0;
    if (shift_t0 > 0) {
      if (!quiet) printf("*----* shifting t0 by %d\n", shift_t0);
      for (i=0; i<MEAS_SEGS; i++) {
	for (j=0; j<TIME_STEPS-shift_t0; j++) {
	  asig_temp.signal[i][j] = asig_temp.signal[i][j+shift_t0];
	}
	for (j=TIME_STEPS-shift_t0; j<TIME_STEPS-1; j++) {
	  asig_temp.signal[i][j] = asig_temp.signal[i][TIME_STEPS-1];
	}
      }
    } else {
      shift_t0 = 0;
    }

    /* do 2-interation grid search in segment seg[jseg] */
    /* CHECKME min. energy fract */
    if (2 != decompose_1(&asig_temp, &f_bsig[jfit],
			 seg[jseg], &f_ints[jfit][jint], &f_t0[jfit], &chisq[jfit],
			 0, 0, 1, 0, 0, 0, 0, 0, 0.01)) {
      printf("*** ACK! decompose_1 gave number of interactions != 2...\n"
	     "         ... expect a crash very soon ...\n");
    } 
    if (!quiet) {
      printf("** grid1[%d]: chisq, pars: %7.4f", jseg, chisq[jfit]);
      for (i=jint; i<jint+2; i++) printf(", %6.3f %6.3f %6.3f %6.3f",
					 f_ints[jfit][i].r, f_ints[jfit][i].p,
					 f_ints[jfit][i].z, f_ints[jfit][i].e);
      printf("\n");
    }
    for (i=jint + 1; i<nints[best]; i++) {
      f_ints[jfit][i+1] = f_ints[best][i];
    }
    f_t0[jfit] = (f_t0[jfit] + f_t0[best] + (double) shift_t0)/2.0;
    nints[jfit] = nints[best] + 1;
    chisq2 = eval_int_pos(asig, &f_bsig[jfit], f_ints[jfit], f_t0[jfit], nints[jfit]);
    if (!quiet) {
      if (fabs(chisq[jfit]-chisq2) < 0.00001) {
	printf("***** chisq[%d] == chisq2: %f %f\n", jfit, chisq[jfit], chisq2);
      } else {
	printf("***** chisq[%d] != chisq2: %f %f\n", jfit, chisq[jfit], chisq2);
      }
    }
    chisq[jfit] = chisq2;
    if (chisq[jfit] < PENALTY3 * chisq[best]) best = jfit;

    /* do nonlinear least-squares fit with the extra interaction(s)*/
    for (i=0; i<nints[jfit]; i++) {
      f_ints[jfit+1][i] = f_ints[jfit][i];
    }
    /* f_t0[jfit+1] = f_t0[jfit]; */
    f_t0[jfit+1] = f_t0[best];
    nints[jfit+1] = nints[jfit];
    jfit++;
    chisq[jfit] = fitter(asig, &f_bsig[jfit], f_ints[jfit], &f_t0[jfit], nints[jfit], 0);
    if (chisq[jfit] < PENALTY4 * chisq[best]) best = jfit;

    if (!quiet) {
      printf("**>>> jint, jfit, best, chisq[best], nints[best] : %d %d %d %f %d\n",
	     jint, jfit, best, chisq[best], nints[best]);
    }
    jfit++;
  }

  if (coalesce) {
    /* Try coalescing interactions, accepting/rejecting by chisq value */
    /* loop throught the hit segments, looking for double interactions */
    for (jseg=0; jseg<nseg; jseg++) {
      m = n = -1;
      for (i = 0; i < nints[best]; i++) {
	f_ints[jfit][i] = f_ints[best][i];
	if (f_ints[best][i].seg == seg[jseg]) {
	  if (m < 0) {
	    m = i;
	  } else {
	    n = i;
	    break;
	  }
	}
      }
      if (n >= 0) {
	/* have found 2 interactions in this seg */
	if (!quiet)
	  printf("*** Try a coalescence fit for seg %d, ints %d, %d\n",
		 seg[jseg], m, n);
	e1 = f_ints[best][m].e;
	e2 = f_ints[best][n].e;
	f_ints[jfit][m].r = (e1*f_ints[best][m].r + e2*f_ints[best][n].r) / (e1+e2);
	f_ints[jfit][m].p = (e1*f_ints[best][m].p + e2*f_ints[best][n].p) / (e1+e2);
	f_ints[jfit][m].z = (e1*f_ints[best][m].z + e2*f_ints[best][n].z) / (e1+e2);
	f_ints[jfit][m].e = e1+e2;
	for (i = n; i < nints[best]-1; i++) {
	  f_ints[jfit][i] = f_ints[best][i+1];
	}
	f_t0[jfit] = f_t0[best];
	nints[jfit] = nints[best] - 1;

	chisq[jfit] = fitter(asig, &f_bsig[jfit], f_ints[jfit], &f_t0[jfit], nints[jfit], 0);
	if (PENALTY3 * chisq[jfit] <= chisq[best]) best = jfit;

	if (!quiet) {
	  printf("**>>> jfit, best, chisq[best], nints[best] : %d %d %f %d\n",
		 jfit, best, chisq[best], nints[best]);
	}
	jfit++;
      }
    }
  }

  /* do a final fit with tighter convergence criteria */
  for (i=0; i<nints[best]; i++) {
    f_ints[jfit][i] = f_ints[best][i];
  }
  f_t0[jfit] = f_t0[best];
  nints[jfit] = nints[best];
  chisq[jfit] = fitter(asig, &f_bsig[jfit], f_ints[jfit], &f_t0[jfit], nints[jfit], 1);
  if (chisq[jfit] < chisq[best]) best = jfit;
  
  *bsig = f_bsig[best];
  *t0 = f_t0[best];
  *chisq_out = chisq[best];
  for (i=0; i<nints[best]; i++) {
    ints[i] = f_ints[best][i];
  }

  return nints[best];
} /* decompose_n */
