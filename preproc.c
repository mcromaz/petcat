

/*** preproc.c : routine for GRETINA mode 3 preprocessing - taken from DecompIF.c (11-1) */

/* preProcess() differs from the online version in that:
	1) sbuf is given directly as an argument rather than derived from inptr linked list
	2) save2and3 functionality has been removed
	3)
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "petcat.h"
#include "diagcnt.h"
#include "lh.h"
#include "wrap.h"
#include "mode3io.h"

/* globals for config files */
static pair *map;
static float a0[37], a1[37], delay0[37], delay1[37];
static int gap = 10;
static struct filter *fltr;
static char *filterfile = "";
static char *trgainfile = "";
static float tr_gain[37];
static int evt_len, tr_len;
static int num_seg_evts, num_xtal_evts;
static int *cur_tr;
static unsigned short *ebuf;
static int notFull;
static int crystalId;

#define WITHINL 0

/* filenames stored in decompLib to make library build simpler */
/* need to rationalize these with *filterfile, *trgainfile */
#define DECOMP_PATHLEN 100
static char basisfile[DECOMP_PATHLEN];
static char xtalkparsfile[DECOMP_PATHLEN];
/* char detmapfile[]; */

/* globals for Heather's INL routines */
int iterations = 0;
double baselineFit[2];
double NstartFit[2];
float tau = 0.0002089;
float restingBase[1120] = {0.};
float sigmaBase[1120] = {0.};
int baseSamples[1120] = {0};
float Q[1120] = {0.};
int lastCalcBase[1120] = {0};
float inl[28][10][8192];
float enl[28][3][8192];
float e1 = 0.;
int waveform[2048];	/* global to hold cc waveform to be consistant with INL routines */

int preProcess(unsigned short int *sbuf, int inlentotal, int evt_len, Event_Signal *event, int new, preprocCnt *diagcnt) {

   struct inbuf *inlist=0, *gone;
   long long t = 0;
   int  x, v, e_0, e_1, t0, i, j, k;  /* remne e0, e1 to e_0, e_1 to not conflict w inl globals */
   int ener_cc, unshifted_ener_cc;
   int num_net, segs[36], net2[36][36];
   int locn, vec, ch, id, slot, *a, *b, *cur_tr_0, *saveref;
   float cal_ener_cc, fact, cc_avg;
   int ie, endofbuffer;
   int recordlen;
   unsigned long long mask;
   int segei;
   int prestep, poststep, ii;

   static int *ener_sp;
   static int num_aligned = 0, num_net1=0, net1_cnt=0;
   static int curr = 0;
   static unsigned long long int full, tlast;
   static int oldrecordlen = 0;

   int nonet[TOT_SEGS];	/* segments where charge collection not in time gate */
   int ti;

   int inlid;

   unsigned short int e0, e1;

  //printf("curr = %d\n", curr);
  if (sbuf == 0) {
    /* restart */
    printf("secondary init of preProcess()\n");
    num_xtal_evts = 0;
    num_seg_evts = 0;
    full = 0;
    tlast = 0;
    tr_len = evt_len - HDR_LEN;
    cur_tr = Calloc(37 * tr_len, sizeof(int));
    ebuf = Calloc(40 * evt_len, sizeof(short));
    ener_sp = Calloc(37 * 4096, sizeof(int));

    return 0;
  }

  /* block omitted which prev processed inlist */

  /* comments refer to the behavior tree diagram in preprocessing.isf */
  /* pp2 */
  if (0 == sbuf) {
    printf("preProcess calling phase error 2 \n");
    return -1;
  }
  //printf("inlentotal = %d, evt_len = %d\n", inlentotal, evt_len);

  if (new) {	/* new sbuf */
    curr = 0;
  }

  while (curr + evt_len <= inlentotal) {
    /* pp23 */
    if (sbuf[curr] != 0xaaaa || sbuf[curr + 1] != 0xaaaa) {
      curr++;
      printf("!");
      continue;	/* up to pp2 */
    }
    diagcnt->numMode3Evts++;
    /* pp10 */
    /*      locn = sbuf[0];*/
    locn = sbuf[3];
    vec = get.vec(sbuf + curr + 2);
    ch = get.ch(sbuf + curr + 2);
    //printf("---vec = %d, ch = %d\n", vec, ch);
    id = (vec << 4) + ch;
    e0 = *(sbuf + curr + 2);
    e1 = *(sbuf + curr + 3);
    //printf("e0 = 0x%x, e1 = 0x%x\n", e0, e1);
    //printf("vec = %d, ch = %d, id = %d, curr = %d\n", vec, ch, id, curr);
/*
      if (curr == 0) {
         char logbuf[80];
         t=get.time(sbuf + curr + 2);
         sprintf(logbuf, " %24lld first event in buf", t);
         sendLogMsg("Prep", logbuf);
      }
*/
    /* pp11 */
    if (id > 4095) {
      /* pp12 */
      curr += evt_len;
      printf("illegal id %d\n", id);
      continue; /* up to pp2 */
    }
    /* pp14 */
    slot = map[id].a;
    //printf("slot = %d\n", slot);
    if (slot == -1) {
      /* pp12 equivalent */
      curr += evt_len;
      printf("illegal slot with id %d\n", id);
      continue; /* up to pp2 */
    }
    num_seg_evts++;
    t=get.time(sbuf + curr + 2);
    tlast =  tlast?tlast:t;
    /* pp22 */
    endofbuffer = 0;
    if (abs(t - tlast) <= gap) {	/* no gap */
      /* pp21 */
				 //printf("there's no gap\n");
      memcpy(ebuf + slot * evt_len, sbuf+ curr, evt_len * 2);
      mask = 1;
      mask <<= slot;
      full |= mask;
      tlast = t;
      curr += evt_len;
      //printf("!");
      if (curr < inlentotal) {
        continue;	/* up to pp2 */
      } else {
        printf("*** endofbuffer ***\n");
        endofbuffer = 1;
      }
    }
    /* pp27 */
    //printf("--gap\n");

    if (full != FULL) {	/* there is a time gap but the event is not complete, chuck prev evt and start new one */
/*
         char logbuf[80];
         sprintf(logbuf, " %24lld full is 0x%llx", t, full);
         sendLogMsg("Prep", logbuf);
*/
      if (full != 0) {
        diagcnt->numPartial++;
        (diagcnt->full)[cntbits(full)]++;
        notFull++;
        /*printf("preProcess bufferend is %d full is 0x%llx\n",
                                  endofbuffer, full);*/
      }
      /* pp28 */
      full = 0;
      /* pp21 equivalent */
      mask = 1;
      mask <<= slot;
      full |= mask;
      if (!endofbuffer) {
        memcpy(ebuf + slot * evt_len, sbuf+ curr, evt_len * 2);
        tlast = t;
        curr += evt_len;
      }
      continue;	/* up to pp2 */
    }
    /* pp19  - there is a gap and the event is complete in ebuf */
    diagcnt->numFull++;
    (diagcnt->full)[cntbits(full)]++;
    if (num_xtal_evts == 0) {
      char logbuf[80];
      sprintf(logbuf, "first built time: %llu CN %d DP %d", tlast,
                         (locn & 0xc0) >> 6, (locn & 0x1f00) >> 8);
      crystalId = locn >> 6;
      sendLogMsg("Prep", logbuf);
      sendToLog(0);
    }
    num_xtal_evts++;
    full = 0;
    /* pp31 seems redundant so omitted */
    /* pp32 */
    ener_cc = segevt_ener(ebuf + evt_len * 36, 7); /* MC - randomize float conversion */
    //printf("ener_cc = %d\n", ener_cc);
    lh_incr(diagcnt->ener_cc, ener_cc);
    cal_ener_cc = a0[36] + a1[36] *
                    ( ((double) ener_cc) + drand48() - 0.5);
    // printf("ener_cc = %d, cal_ener_cc = %f\n", ener_cc, cal_ener_cc);
    /* pp34 ignored in code */
    /* pp36 */

    for (i = 0, num_net = 0; i < TOT_SEGS; i++) {
      /* pp43   */
      a = cur_tr + i *tr_len;
      get.tr(a, tr_len, 1, ebuf + evt_len *i +2);
      adjoff(a, tr_len);
      ti = time_sef(a, tr_len);
      /* pp45 */
      // printf("*************** i = %d ***************\n", i);
      if (!net(a, tr_len)  || (ti < ((int) (0.1 * tr_len))) || (ti > ((int) (0.9 * tr_len)))) {	/* MC - hack for high rate, ext rigger in BGS runs */
        nonet[i] = 1;	/* mark segment with no net so segment energy can be zeroed later  */
        continue;	/* up to pp36 */
      }
      nonet[i] = 0;
      /* pp46 */
      segs[num_net] = i;
      num_net++;
    }

    /* pp42 */
    a = cur_tr + i * tr_len;
    get.tr(a, tr_len, 0, ebuf + evt_len *i +2);

    #ifdef SAMPLE25
    cur_tr_0 = Calloc(37 * (tr_len / 2), sizeof(int));
    for (i = 0; i < 37; i++) {
      a = cur_tr + i * tr_len;
      b = cur_tr_0 + i * tr_len / 2;
      for (j = 0; j < tr_len / 2; j++) {
        b[j] = (a[2 * j] + a[2 * j + 1]) / 2;
      }
    }
    saveref = cur_tr;
    cur_tr = cur_tr_0;  // cur_tr now points to compressed traces
    tr_len /= 2;
    #endif

    prestep = 0;
    poststep = 0;
    if (tr_len >= 160) {
      for (ii = 40; ii < 60; ii++) {
        prestep += a[ii];
      }
      for (ii = 140; ii < 160; ii++) {
        poststep += a[ii];
      }
    }

    /* pp46 */
    /* inl calc from Heather's mode 3 sort */
    if(WITHINL) {
      // inlid = inl_id(sbuf + curr + 2);
      inlid = inl_id(ebuf + (evt_len * 36) + 2);
      for (ii = 0; ii < tr_len; ii++) {
	waveform[ii] = a[ii];
      }
      iterations = LSFitExponential(1, 60, waveform[0], tau, &NstartFit[0], &baselineFit[0], &NstartFit[1], &baselineFit[1]);

      if ( (fabs(baselineFit[0] - restingBase[inlid]) <
	    fabs(sigmaBase[inlid]*3)) ||
	   (baseSamples[inlid] == 0) || (sigmaBase[inlid] == 0) ) {

	baseSamples[inlid]++;
	Q[inlid] = (Q[inlid] +
		    (((float)(baseSamples[inlid]-1)/(float)baseSamples[inlid])*
		     (baselineFit[0] - restingBase[inlid])*
		     (baselineFit[0] - restingBase[inlid])));
	restingBase[inlid] = (restingBase[inlid] +
			      ((baselineFit[0] - restingBase[inlid])/baseSamples[inlid]));
	sigmaBase[inlid] = sqrt(Q[inlid]/baseSamples[inlid]);
      }

      /* If the baseline has changed, recalculate the ENL values. */

      if ((int) restingBase[inlid] != lastCalcBase[inlid]) {
	calculateENL(inlid / 40, inlid % 10,
		     restingBase[inlid], 20, 20,
		     prestep/20., poststep/20.);
	lastCalcBase[inlid] = (int) restingBase[inlid];
      }

      unshifted_ener_cc = segevt_ener(ebuf + evt_len * 36, 0);

      /* Remember, prestep is the base, poststep is the flat-top on the trace. */
      // printf("%d %f ", poststep, a0[36] + a1[36] * ( (((double)unshifted_ener_cc/128.)) + drand48() - 0.5));

      unshifted_ener_cc += (enl[inlid/40][0][(int)poststep/20+2000] -
			    enl[inlid/40][1][(int)prestep/20+2000]);
      ener_cc = unshifted_ener_cc / 128;

      cal_ener_cc = a0[36] + a1[36] *
	( ((double) ener_cc) + drand48() - 0.5);
    }

    adjoff(a, tr_len);
    /* pp47 */
    x = ener_cc;
    if (x >= fltr->xmin && x <= fltr->xmax && num_net == 2) {
         e_0 = segevt_ener(ebuf + evt_len * segs[0], 7);
         e_1 = segevt_ener(ebuf + evt_len * segs[1], 7);
         if (e_0 >= e_1) {
            net2[segs[0]][segs[1]]++;
         } else {
            net2[segs[1]][segs[0]]++;
         }
      }
      /* pp58 */
      if ( ((cal_ener_cc < fltr->emin ||
             cal_ener_cc > fltr->emax)) ||
	   /* 2012-06-24 CMC removed x filter, Raw energy */
	   /*           ((x < fltr->xmin || x > fltr->xmax)) || */
	   /*            (num_net != 1  && num_net != 2) || */
	   (!(num_net>0)) || /* TODO: should there be any such cut here? */
           segs[0] < fltr->segmin || segs[0] > fltr->segmax) {
         //printf("bailed at pp58\n");
         //printf("cal_ener_cc = %f, fltr->emin = %f, fltr->emax = %f\n", cal_ener_cc, fltr->emin, fltr->emax);
         diagcnt->pp58++;
         continue; /* up to pp2 */
      }
      /* pp59 */
      num_net1++;
      t0 = align_cfd_1(cur_tr, tr_len, delay1, delay0);
      v = t_cfd(cur_tr + 36 * tr_len, tr_len);  /* DCR: Why this? What is v used for? */
      /* pp61 */
      /*if (t0 < 0) { TODO: Do we need this?
         net1_cnt++;
         continue; up to pp2
      }*/
      num_aligned++;
      fact = ((float)x)/3.3639; /* TODO: this depends on integration time! */
      /* comment out for petcat
      event = calloc(1, sizeof(Event_Signal));
      if (!event) return -1;
      */
#ifdef TIMEOFFSET
      /* DCR: save offset coming from time alignment of signals */
      if (t0 > 0) {
	event->time_offset = 0.01f * (float) t0;
      } else {
	event->time_offset = 0.0f;
      }
#endif

      event->prestep = ((float) prestep) / 20.0;
      event->poststep = ((float) poststep) / 20.0;

      /* pp89 */
      for (i = 0;  i < TOT_SEGS; i++) {
         /* pp71 */
         /* cc_avg not used ..
         for (v=45,cc_avg=0; v < 50; v++) {
            cc_avg += (float) (*(cur_tr + 36 * tr_len + v));
         }
         */
         /* pp72 */
         segei = segevt_ener(ebuf + evt_len * i, 7);	/* MC randomize float conversion */
         event->seg_energy[i]  = a0[i] + a1[i] * ((double) segei);
	 /* segment time gate failed, zero segment energy to
	    prevent decomp fitting segment */
         if (nonet[i] == 1) { event->seg_energy[i] = 0.0; }
         /*
         event->seg_energy[i] = a0[i] + a1[i] *
              ((float) segevt_ener(ebuf + evt_len * i, 7));
         */
         ie  = (int) (event->seg_energy[i] - 0.5 + drand48());
         /* pp73 */
         if (ie >= 0 && ie < 4096) {
           ener_sp[i * 4096 + ie]++;
           (diagcnt->seg)[i][ie]++;

         }
         /* pp76 */
         for (j = 0; j < TIME_STEPS; j++) {
	   event->signal[i][j] = (1.0 /
				  (tr_gain[i] * fact)) * ((float) (*(cur_tr + i * tr_len + j)));
         }
         /* pp88 */
     }
     /* pp82 */
     /*event->total_energy = a0[36] + a1[36] * ((float) segevt_ener(ebuf + evt_len * 36, 7));*/
     event->total_energy = cal_ener_cc;
     ie  = (int) (event->total_energy - 0.5 + drand48());
     if (ie >= 0 && ie < 4096) {
       ener_sp[36 * 4096 + ie]++;
     }
     /* pp84 */
     for (k = 0; k < TIME_STEPS; k++) {
       event->signal[36][k] = - (1.0 / (tr_gain[36] * fact)) *
	 ((float) (*(cur_tr + 36 * tr_len + k)));
     }

     /* copy 4 core energies, detmap slot 36 first followed by remaining 3 */
     for (k = 0; k < 4; k++) {
	event->core_e[k] = segevt_ener(ebuf + evt_len * (k + 36), 7);
     }

     /* pp87 */
     event->time = tlast;
     //printf("--- returned 0\n");

     #ifdef SAMPLE25
     cur_tr = saveref; // restore cur_tr to full buffer
     free(cur_tr_0);
     #endif

     return 0;
   }
/*
   {
         char logbuf[80];
         sprintf(logbuf, " %24lld last event in buf", t);
         sendLogMsg("Prep", logbuf);
         sprintf(logbuf, "full 0x%llx of 0x%llx ", full, FULL);
         sendLogMsg("Prep", logbuf);
   }
*/

   printf("--- returned -1\n");
   return -1;
}

int preProcessMario(Mario *mario, Event_Signal *event, preprocCnt *diagcnt) {

   struct inbuf *inlist=0, *gone;
   long long t = 0;
   int  x, v, e_0, e_1, t0, i, j, k;  /* remne e0, e1 to e_0, e_1 to not conflict w inl globals */
   int m;
   int ener_cc, unshifted_ener_cc;
   int num_net, segs[36], net2[36][36];
   int locn, vec, ch, id, slot, *a;
   float cal_ener_cc, fact, cc_avg;
   int ie, endofbuffer;
   int recordlen;
   unsigned long long mask;
   int segei;
   int prestep, poststep, ii;

   static int *ener_sp;
   static int num_aligned = 0, num_net1=0, net1_cnt=0;
   static int curr = 0;
   static unsigned long long int full, tlast;
   static int oldrecordlen = 0;

   int nonet[TOT_SEGS];	/* segments where charge collection not in time gate */
   int ti;

   int inlid;
   FILE *ftmp;
   /* DO I NEED TO DO THIS? */
   /*   printf("secondary init of preProcess()\n");
   num_xtal_evts = 0;
   num_seg_evts = 0;
   full = 0;
   tlast = 0;
   tr_len = evt_len - HDR_LEN;

   ebuf = Calloc(40 * evt_len, sizeof(short));
   ener_sp = Calloc(37 * 4096, sizeof(int)); */
   tr_len = 286;
   #ifdef SAMPLE25
   for (i = 0; i < 37; i++) {  // compress trace to 1/2 length by summing
     for (j = 0; j < tr_len / 2; j++) {
       mario->wf[i][j] = (mario->wf[i][2 * j] + mario->wf[i][2 * j + 1]) / 2;
     }
     for (j = tr_len / 2; j < tr_len; j++) {
       mario->wf[i][j] = 0;
     }
   }
   tr_len /= 2;
   #endif

/*
   ftmp = fopen("tr15.csv", "w");
   for (i = 0; i < tr_len; i++) {
     fprintf(ftmp, "%d,%d\n", i, mario->wf[15][i]);
   }
   exit(1);
*/
   cur_tr = Calloc(37 * tr_len, sizeof(int));

   diagcnt->numMode3Evts++;

   diagcnt->numFull++;
   num_xtal_evts++;
   full = 0;

   ener_cc = mario->ccEnergy;
   cal_ener_cc = mario->ccEnergy;

   for (i = 0, num_net = 0; i < TOT_SEGS; i++) {

     a = cur_tr + i *tr_len;
     for (m=0; m<tr_len; m++) {
       a[m] = (int)mario->wf[i][m];
     }
     adjoff(a, tr_len);
     ti = time_sef(a, tr_len);
     /*
     if (i == 15) {
     printf("ti = %d\n", ti);
     exit(1);
     }
     */
     if (!net(a, tr_len)  || (ti < ((int) (0.1 * tr_len))) || (ti > ((int) (0.9 * tr_len)))) {	/* MC - hack for high rate,
												   ext rigger in BGS runs */
       nonet[i] = 1;	/* mark segment with no net so segment energy can be zeroed later  */
       continue;	/* up to pp36 */
     }
     nonet[i] = 0;
     segs[num_net] = i;
     num_net++;
   }

   a = cur_tr + i * tr_len;
   for (m=0; m<tr_len; m++) {
     a[m] = (int)mario->wf[i][m];
   }

   prestep = 0;
   poststep = 0;
   if (tr_len >= 160) {
     for (ii = 40; ii < 60; ii++) {
       prestep += a[ii];
     }
     for (ii = 140; ii < 160; ii++) {
       poststep += a[ii];
     }
   }
   adjoff(a, tr_len);

   x = (int)(((ener_cc*128) - a0[36])/(a1[36]));
   if ( ((cal_ener_cc < fltr->emin ||
	  cal_ener_cc > fltr->emax)) ||
	(!(num_net>0)) || /* TODO: should there be any such cut here? */
	segs[0] < fltr->segmin || segs[0] > fltr->segmax) {
     diagcnt->pp58++;
   }

   num_net1++;
   t0 = align_cfd_1(cur_tr, tr_len, delay1, delay0);
   //printf("t0=%d, num_net=%d %d %d\n", t0, num_net, segs[0], segs[1]);
   //exit(1);
   num_aligned++;
   fact = ((float)x)/3.3639; /* TODO: this depends on integration time! */

#ifdef TIMEOFFSET
   if (t0 > 0) {
     event->time_offset = 0.01f * (float) t0;
   } else {
     event->time_offset = 0.0f;
   }
#endif
   event->prestep = ((float) prestep) / 20.0;
   event->poststep = ((float) poststep) / 20.0;

   for (i = 0;  i < TOT_SEGS; i++) {
     /* cc_avg not used
     for (v=45,cc_avg=0; v < 50; v++) {
       cc_avg += (float) (*(cur_tr + 36 * tr_len + v));
     }
     */

     event->seg_energy[i]  = mario->segEnergy[i];
     if (nonet[i] == 1) { event->seg_energy[i] = 0.0; }
     for (j = 0; j < TIME_STEPS; j++) {
       event->signal[i][j] = (1.0 /
			      (tr_gain[i] * fact)) * ((float) (*(cur_tr + i * tr_len + j)));
     }
   }

   /* pp82 */
   event->total_energy = cal_ener_cc;
   for (k = 0; k < TIME_STEPS; k++) {
     event->signal[36][k] = - (1.0 / (tr_gain[36] * fact)) *
       ((float) (*(cur_tr + 36 * tr_len + k)));
   }

   /* copy 4 core energies, detmap slot 36 first followed by remaining 3 */
   for (k = 0; k < 4; k++) {
     event->core_e[k] = cal_ener_cc;
   }
/*
      ftmp = fopen("tr15post.csv", "w");
      for (i = 0; i < TIME_STEPS; i++) {
        fprintf(ftmp, "%d,%f\n", i, event->signal[15][i]);
      }
      exit(1);
*/
   /* pp87 */
   //event->time = tlast;
   return 1;
}

int startPreProcess(int evt_len, char *detMapFilename, char *filterFilename, char *trGainFilename, char*xTalkParsFilename) {

  int num;
  float avg = 0.0;
  FILE *f_tr_gain;
  int i, j, k;

  if (evt_len > 0) {
    ebuf = Calloc(40 * evt_len, sizeof(short));
  } else {
    evt_len = 100;
  }

  map = read_mapfile(detMapFilename, a0, a1);

  if (map == 0) {
    fprintf(stderr, "cannot read detector map file %s\n", detMapFilename);
    return -1;
  }

  fltr = read_filterfile(filterFilename);
  if (0 == fltr) {fprintf(stderr, "cannot open filter file\n"); return -1;}
  f_tr_gain = fopen(trGainFilename, "r");
  if (f_tr_gain == 0) {fprintf(stderr, "cannot open trgain file\n"); return -1;}
  for (i = 0; i < 37; i++) {
    num = fscanf(f_tr_gain, "%f", tr_gain + i);
    if (num != 1) {fprintf(stderr, "too few lines in trgain file\n"); return -1;}
  }
  num = read_param(xTalkParsFilename, "delay0", delay0, TOT_SEGS);
  if (num <  0) {
    fprintf(stderr, "error: cannot read delay0 from %s\n", xTalkParsFilename);
    return -1;
  }
  /* subtract mean value of delay0 from individual values */
  for (i = 0; i < 36; i++) {
    avg += delay0[i];
  }
  avg /= 36.0;
  for (i = 0; i < 36; i++) {
    delay0[i] -= avg;
  }

  num = read_param(xTalkParsFilename, "delay1", delay1, TOT_SEGS);
  if (num <  0) {
    fprintf(stderr, "error: cannot read delay1 from %s\n", xTalkParsFilename);
    return -1;
  }
  delay0[36] = delay1[36] = 0.0;

  /* init arrays for inl */
  /* if(WITHINL) { */
  /*   for (i = 0; i < 1120; i++){ */
  /*     restingBase[i] = 0.0; */
  /*     sigmaBase[i] = 0.0; */
  /*     baseSamples[i] = 0; */
  /*     Q[i] = 0.0; */
  /*     lastCalcBase[i] = -10000; */
  /*   } */
  /*   for (i = 0; i < 28; i++) { */
  /*     for (j = 0; j < 10; j++) { */
  /* 	for (k = 0; k < 8192; k++) { */
  /* 	  inl[i][j][k] = 0.0; */
  /* 	} */
  /*     } */
  /*   } */
  /*   for (i = 0; i < 28; i++) { */
  /*     for (j = 0; j < 3; j++) { */
  /* 	for (k = 0; k < 8192; k++) { */
  /* 	  enl[i][j][k] = 0.0; */
  /* 	} */
  /*     } */
  /*   } */
  /*   ReadRawINL(); */
  /* } 	 */

  preProcess(0,0,evt_len,0, 0,0);	/* reset internals */
  /* queueFlush(0); */

  notFull = 0;

  return 0;
}

void ReadRawINL() {

  FILE *digMap;
  FILE *inl_in;
  int dorder, dsn, i1, dhole;
  char *st, str[256];
  int crystalSN[28];
  char serialN[20];
  char inlName[256];
  int i, j, m;

  for (i = 0; i < 28; i++) {
    crystalSN[i] = 0;
  }

  if ((digMap = fopen("/home/users/crawford/PetCat/gretina-petcat/Digitizers.2012-07-11.MSUJuly9_Commissioning_I", "r")) == NULL) {
    printf("Failed to open digitizer SN map.\n");
  } else {
    printf("Reading digitizer map \"Digitizers.2012-07-11.MSUJuly9_Commissioning_I\".\n");
  }

  st = fgets(str, 200, digMap);
  while (st!= NULL) {
    if (str[0] == 35) {
			;
    } else if (str[0] == 59) {
			;
    } else if (str[0] == 10) {
			;
    } else {
      sscanf(str, "%i %i %X %X", &dorder, &i1, &dsn, &dhole);
      if (dorder%4 == 1) {
	if (i1 == 0) { crystalSN[dorder/4] = dsn; }
	if (i1 == 2) { crystalSN[dorder/4] = dsn | (2 << 10); }
	if (i1 == 1) { crystalSN[dorder/4] = dsn | (1 << 10); }
	if (i1 == 3) { crystalSN[dorder/4] = dsn | (3 << 10); }
      }
    }
    st = fgets(str, 200, digMap);
  }

  for (i=0; i<28; i++) {

    sprintf(serialN, "%X", crystalSN[i]);

    sprintf(inlName, "%s%X%s", "/home/users/crawford/INL/INL-0x", crystalSN[i], ".dat");


    if ((inl_in = fopen(inlName, "r")) == NULL) {
      printf("Failed to open INL file %s\n", inlName);
      printf("Skipping it and moving to the next one...\n");
    } else {
      printf("Reading INL corrections from %s\n", inlName);
      fread(inl[i][0], sizeof(inl), 1, inl_in);
      fclose(inl_in);
    }

    for (j=0; j<10; j++) {
      for (m=0; m<8192; m++) {
	inl[i][j][m] = -2.0f*(inl[i][j][m] + 0.0f);
      }
    }

  }

}

int inl_id(unsigned short int *x) {
  /* calc id as required for INL correction - taken from Mode3 sort */
  unsigned short int board_id, chan_id, module, id;

  board_id = (x[0] >> 11);
  chan_id = (x[1] & 0xf);
  module = (x[1] >> 4);

    if (module/16 == Q1) {
      id = (module - (Q1-Q1E)*4*4)*10 + chan_id;
    } else if (module/16 == Q2) {
      id = (module - (Q2-Q2E)*4*4)*10 + chan_id;
    } else if (module/16 == Q3) {
      id = (module - (Q3-Q3E)*4*4)*10 + chan_id;
    } else if (module/16 == Q4) {
      id = (module - (Q4-Q4E)*4*4)*10 + chan_id;
    } else if (module/16 == Q5) {
      id = (module - (Q5-Q5E)*4*4)*10 + chan_id;
    } else if (module/16 == Q6) {
      id = (module - (Q6-Q6E)*4*4)*10 + chan_id;
    } else if (module/16 == Q7) {
      id = (module - (Q7-Q7E)*4*4)*10 + chan_id;
    } else {
      id = (module)*10 + chan_id;
    }
    return (int) id;
}

double ComputeChiSquare(int startIndex, int nPoints, double tau, double Nstart, double baseline) {
  double chiSq = 0;
	int i;

  for (i=startIndex; i<(startIndex+nPoints); i++) {
    chiSq += ((waveform[i] - (Nstart*exp(-((i-startIndex)/tau))
				   + baseline))*
	      (waveform[i] - (Nstart*exp(-((i-startIndex)/tau))
				   + baseline)));
  }
  return chiSq;
}

int LSFitExponential(int startIndex, int nPoints, double Nstart, double tau, double *NstartFit, double *baselineFit, double *NstartFitError, double *baselineFitError) {

  /* Do a Levenberg-Marquardt minimization of a exponential decay + offset
     with a fixed tau, and get N0 and baseline(offset) value.

     Equation: N(t) = Nstart*exp(-(t-t0)/tau) + baseline

     We will define all matrices as:
        (  Nstart,Nstart     Nstart,baseline )
        ( baseline,Nstart   baseline,baseline)             */

  /* Set-up starting point for minimization...
     starting parameter, lambda, etc. */
  double baselineStart = Nstart - 150; /* Assume Nstart-10 is a decent guess... */
  double NstartStart = Nstart;
  double baselineNow, baselineLast;
  double NstartNow, NstartLast;

  double lambda = 0.001;
  double chiSqNow;
  double chiSqLast;

  double alpha[2][2] = {{0,0},{0,0}};
  double beta[2] = {0, 0};
  double error[2][2] = {{0,0},{0,0}};
  double stepping[2] = {0, 0};

  int notConverged = 1;
  int iterations = 0;

	int i;

  baselineNow = baselineStart;
  NstartNow = NstartStart;

  chiSqNow = ComputeChiSquare(startIndex, nPoints, tau, NstartNow, baselineNow);

  while (notConverged && iterations<10) {
    /* Clear the alpha, beta and error arrays... */
    alpha[0][0] = 0;    alpha[0][1] = 0;    alpha[1][0] = 0;    alpha[1][1] = 0;
    beta[0] = 0;    beta[1] = 0;
    error[0][0] = 0;    error[0][1] = 0;    error[1][0] = 0;    error[1][1] = 0;
    stepping[0] = 0;    stepping[1] = 0;

    /* Compute the curvature matrix and extremum vector first. */
    for (i=startIndex; i<(startIndex+nPoints); i++) {
      alpha[0][0] += exp(-2*(i-startIndex)/tau);
      alpha[0][1] += exp(-(i-startIndex)/tau);
      alpha[1][0] += exp(-(i-startIndex)/tau);
      alpha[1][1] += 1;

      beta[0] += ((waveform[i] -
		   (NstartNow*exp(-(i-startIndex)/tau)) -
		   baselineNow)*
		  (exp(-(i-startIndex)/tau)));
      beta[1] += ((waveform[i] -
		   (NstartNow*exp(-(i-startIndex)/tau)) -
		   baselineNow));
    }

    alpha[0][0] = alpha[0][0]*(1+lambda);
    alpha[1][1] = alpha[1][1]*(1+lambda);

    /* Now compute the error matrix and stepping vector. */
    error[0][0] = (1/((alpha[0][0]*alpha[1][1]) - (alpha[0][1]*alpha[0][1])))*alpha[1][1];
    error[0][1] = (1/((alpha[0][0]*alpha[1][1]) - (alpha[0][1]*alpha[0][1])))*(-alpha[0][1]);
    error[1][0] = (1/((alpha[0][0]*alpha[1][1]) - (alpha[0][1]*alpha[0][1])))*(-alpha[1][0]);
    error[1][1] = (1/((alpha[0][0]*alpha[1][1]) - (alpha[0][1]*alpha[0][1])))*alpha[0][0];

    stepping[0] = error[0][0]*beta[0] + error[0][1]*beta[1];
    stepping[1] = error[1][0]*beta[0] + error[1][1]*beta[1];

    /* Compute new values for variables. */
    NstartLast = NstartNow;
    NstartNow = NstartLast + stepping[0];
    baselineLast = baselineNow;
    baselineNow = baselineLast + stepping[1];

    chiSqLast = chiSqNow;
    chiSqNow = ComputeChiSquare(startIndex, nPoints, tau, NstartNow, baselineNow);

    /* Evaluate our progress... */
    if (chiSqNow < chiSqLast) { /* Success, decrease lambda */
      lambda /= 10;
    } else {
      lambda *= 10;
      NstartNow = NstartLast;
      baselineNow = baselineLast;
    }

    iterations++;

/*
    // cout << "Conditions: iter " << iterations << endl;
    // cout << "  Baseline now " << baselineNow << " Nstart now " << NstartNow << endl;
    // cout << "  Curvature:   " << alpha[0][0] << "     " << alpha[0][1] << endl;
    // cout << "               " << alpha[1][0] << "     " << alpha[1][1] << endl;
    // cout << "  Extremum:    " << beta[0] << "     " << beta[1] << endl;
    // cout << "  Error:       " << error[0][0] << "     " << error[0][1] << endl;
    // cout << "               " << error[1][0] << "     " << error[1][1] << endl;
    // cout << "  Stepping:    " << stepping[0] << "     " << stepping[1] << endl;
    // cout << "  Lambda       " << lambda      << " Chi2 " << chiSqNow << endl;
*/

    if ( (fabs(stepping[0]) < (0.00001*sqrt(error[0][0]))) &&
	 (fabs(stepping[1]) < (0.00001*sqrt(error[1][1]))) ) {
      notConverged = 0;
    }
  }

  /* We're done, get the final answer... */
  alpha[0][0] = 0;    alpha[0][1] = 0;    alpha[1][0] = 0;    alpha[1][1] = 0;
  beta[0] = 0;    beta[1] = 0;
  error[0][0] = 0;    error[0][1] = 0;    error[1][0] = 0;    error[1][1] = 0;
  stepping[0] = 0;    stepping[1] = 0;

  /* Compute the curvature matrix and extremum vector first. */
  for (i=startIndex; i<(startIndex+nPoints); i++) {
    alpha[0][0] += exp(-2*(i-startIndex)/tau);
    alpha[0][1] += exp(-(i-startIndex)/tau);
    alpha[1][0] += exp(-(i-startIndex)/tau);
    alpha[1][1] += 1;

    beta[0] += ((waveform[i] -
		 (NstartNow*exp(-(i-startIndex)/tau)) -
		 baselineNow)*
		(exp(-(i-startIndex)/tau)));
    beta[1] += ((waveform[i] -
		 (NstartNow*exp(-(i-startIndex)/tau)) -
		 baselineNow));
  }

  /* Now compute the error matrix and stepping vector. */
  error[0][0] = (1/((alpha[0][0]*alpha[1][1]) - (alpha[0][1]*alpha[0][1])))*alpha[1][1];
  error[0][1] = (1/((alpha[0][0]*alpha[1][1]) - (alpha[0][1]*alpha[0][1])))*(-alpha[0][1]);
  error[1][0] = (1/((alpha[0][0]*alpha[1][1]) - (alpha[0][1]*alpha[0][1])))*(-alpha[1][0]);
  error[1][1] = (1/((alpha[0][0]*alpha[1][1]) - (alpha[0][1]*alpha[0][1])))*alpha[0][0];

  *NstartFit = NstartNow;
  *NstartFitError = sqrt(error[0][0]);
  *baselineFit = baselineNow;
  *baselineFitError = sqrt(error[1][1]);

  return iterations;

}

void calculateENL(int xtal, int chn, float base, int s1_width, int s2_width, float s1, float s2) {

  float ss1=0., ss2=0., ss3=0., ss4=0.;
  float sum1=0., sum2=0., sum3=0., sum4=0.;
	int i, j;

  /* Clear the ENL correction from before. */
  for (i=0; i<2; i++) {
    for (j=0; j<8192; j++) {
      enl[xtal][i][j] = 0;
    }
  }

  for (i=0; i<8192; i++) {
    e1 = i*0.5 - 500.;
    if (chn != 9) { e1 = -e1; }
    sum1 = sum2 = sum3 = sum4 = 0;

    for (j=0; j<500; j++) {
      sum1 += e1 + base;
      sum2 += inl[xtal][9][(int)((8192.0 + e1 + base)/2.0)];

      if (j < s2_width) {
	sum3 += (e1 + base - inl[xtal][9][(int)((8192.0 + e1 + base)/2.0)]);
      }
      if (j >= (500-s1_width)) {
	sum4 += (e1 + base - inl[xtal][9][(int)((8192.0 + e1 + base)/2.0)]);
      }
      e1 *= (1.0 - (tau));
    }

    if (chn != 9) {
      sum3 /= -s2_width;
      sum4 /= -s1_width;
    } else {
      sum3 /= s2_width;
      sum4 /= s1_width;
    }

    sum1 += 2000.0 - sum2;
    sum3 += 2000.0;
    sum4 += 2000.0;

    if (i>0 && ss3>0 && ss3<sum3 && sum3<8190) {
      for (j = 1 + (int) ss3; (e1 = j) <= sum3; j++) {
	enl[xtal][0][j] = ss2 + (sum2 - ss2)*(e1 - ss3)/(sum3 - ss3);
      }
    }
    if (i>0 && ss4>0 && ss4<sum4 && sum4<8190) {
	for (j = 1 + (int) ss4; (e1 = j) <= sum4; j++) {
	  enl[xtal][1][j] = ss2 + (sum2 - ss2)*(e1 - ss4)/(sum4 - ss4);
	}
    }
    ss1 = sum1;
    ss2 = sum2;
    ss3 = sum3;
    ss4 = sum4;
  }

  //for (i=0; i<8192; i++) {
  // printf("%d, %f %f \n", i, enl[xtal][0][i], enl[xtal][1][i]);
  //}

}
