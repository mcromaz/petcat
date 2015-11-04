
/*** petcat.h: general definitions for offline preproc/decomp */

#include <stdio.h>
#include "diagcnt.h"

#define TOT_SEGS 36
#define MEAS_SEGS 37
#define MAX_TIME_STEPS 50
#define MAX_TR_LEN 1024
#define FULL 0xffffffffffll   /* 40 bits set - build all digitizer ch */
#define NSIGS 37
#define MAX_TRACE_LEN 1024
#define CFD_INT_LEN 4
#define CFD_DELAY 4
#define CFD_FRACTION 4
#define HDR_LEN 16

#define RBUF_SIZE (1024 * 1024)
/*#define OBUF_LEN (1024 * 1024)*/

#define MAX_3EVT_LEN 2048
#define MAX_EVT_LEN  16384	/* maximum payload length of an evt w a gebData hdr */

/* decomp definitions */
#define MAX_GRID_PTS   250000  /* max number of grid points in a basis */
#define MAX_GRID_SEGS  37 /* max of signals calculated for each basis point */
#define MEAS_SEGS  37      /* number of signals measured for each event */
#define MAX_TIME_STEPS 50      /* max number of time steps calculated/measured for each segment */
#define TOT_SEGS   36      /* total number of segments in crystal */
#define MAX_AGS    1000  /* CHANGEME max. no. of points in coarse grid for AGS */
#define MAX_SEGS   8              /* max. number of segments to take in events */
#define MAX_PARS   (8*MAX_SEGS+1) /* max. number of parameters to fit */
#define MAX_INTPTS (2*MAX_SEGS)   /* max. number of interaction points */
#define SSEG       TOT_SEGS       /* range of seg in basis grid */
#define MAX_SRAD       40         /* max range of r in basis grid */
#define MAX_SPHI       20        /* max range of phi in basis grid */
#define MAX_SZZZ       30        /* max range of z in basis grid */

typedef unsigned short int uint2;

typedef struct {
  int a;
  int b;
} pair;

struct filter {
  float emin, emax;
  int xmin, xmax;
  int segmin, segmax;
};

typedef struct {
  float total_energy;                    /* total energy in crystal */
  float seg_energy[TOT_SEGS];            /* net energy in each segment */
  float signal[MEAS_SEGS][MAX_TIME_STEPS];   /* actual measured signals */
#ifdef TIMESTAMP
  long long int time;                    /* timestamp for this crystal event */
#endif
#ifdef TIMEOFFSET
  float time_offset;  /* offset coming from time alignment of signals */
#endif
  int core_e[4];              			/* 4 core energy signals */
  float prestep;
  float poststep;
} Event_Signal;

extern struct fmt_defn {
  int sentinnel;
  int (*ener)(unsigned short int *);
  long long int (*time)(unsigned short int *);
  int (*vec)(unsigned short int *);
  int (*ch)(unsigned short int *);
  int (*tr)(int *, int, int, unsigned short int *);
  int (*tr_len)(unsigned short int *);
  int (*evt_len)(unsigned short int *);
} fmt_defn[];

extern struct fmt_defn get;     /* global struct of accessor functions */

typedef struct {
  int    seg;                 /* segment */
  int    pos;                 /* basis signal position, if applicable, else -1 */
  double r, p, z, e;          /* parameters */
  double dr, dp, dz, de;      /* uncertainties */
} Interaction;

typedef struct {
  float x, y, z, e;
  float diff;
} simpleInteraction;

typedef struct {
  char  iseg, ir, ip, iz;             /* integer cylindrical coordinates of grid point */
  float x, y, z;                      /* cartesian coordinates of grid point */
  float signal[MAX_GRID_SEGS][50];        /* actual basis signals */
  int   lo_time[MAX_GRID_SEGS], hi_time[MAX_GRID_SEGS];     /* limits for non-zero signal */
} Basis_Point;

typedef struct {
  int    npts;                /* number of AGS points for this segment */
  int    grid_pos[MAX_AGS];   /* pointer to basis-signal-ID for each AGS point */
  double *da;                 /* precalculated sums */
  float  *db;                 /* precalculated sums */
} Adaptive_Grid;

struct crys_intpts {
  int type;          /* defined as abcd5678 */
  int crystal_id;
  int num;           /* # of int pts from decomp, or # of nets on decomp error */
  float tot_e;       /* dnl corrected */
  int core_e[4];     /* 4 raw core energies from FPGA filter (no shift) */
  long long int timestamp;
  long long trig_time;    /* not yet impl */
  float t0;
  float cfd;
  float chisq;
  float norm_chisq;
  float baseline;
  float prestep;    /* avg trace value before step */
  float poststep;   /* avg trace value following step */
  int pad;          /* non-0 on decomp error, value gives error type */
  struct {
    float x, y, z, e;       /* here e refers to the fraction */
    int seg;                /* segment hit */
    float seg_ener;         /* energy of hit segment */
  } intpts[MAX_INTPTS];
};

struct decomp_state {
  Event_Signal *bsig;
  Interaction  *ints;
  int   cnt;
  float coal_dist;
  struct decomp_errcnt *err;
  /* next two are for diagnostic output */
  struct crys_intpts pos;
  short mat[4096];
};

struct decomp_errcnt {
  int nonet;
  int toomanynet;
  int sumener;
  int badchisq;
};

#define EBUF_SIZE 3
#define MAX_LINE_LEN 120

/* detector hole #'s - S800 */
#define Q1 (15)
#define Q2 (6)
#define Q3 (16)
#define Q4 (8)
#define Q5 (9)
#define Q6 (7)
#define Q7 (17)

#define Q1E (0)
#define Q2E (1)
#define Q3E (2)
#define Q4E (3)
#define Q5E (4)
#define Q6E (5)
#define Q7E (6)

extern float inl[28][10][8192];

#define COAL_DIST_DEFAULT 2.0

/* decomp globals */
int GRID_PTS;   /* actual number of grid points in the basis */
int GRID_SEGS;  /* actual number of signals calculated for each basis point */
int TIME_STEPS;  /*  actual time steps calculated/measured for each segment */
int SRAD;                 /* actual range of r in basis grid */
int SPHI;                 /* actual range of phi in basis grid */
int SZZZ;                 /* actual range of z in basis grid */
Basis_Point   *basis;                 /* basis-signal data */
int           grid_pos_lu[SSEG][MAX_SRAD][MAX_SPHI][MAX_SZZZ];    /* basis-grid position look-up table */
int           maxir[SSEG], maxip[SSEG], maxiz[SSEG];  /* max. values of ir, ip, iz for each segment */
Adaptive_Grid ags1[SSEG];             /* Adaptive Grid Search coarse grid, 1 segment */
int           quiet;                  /* set to 1 to prevent output of diagnostic info */
int           average_sigs;           /* set to 1 to output averaged observed and
						fitted signals for single-segment (net=1) events */
int           *bad_segs;              /* list of segment numbers that should be disabled */

/* prototypes */
int startPreProcess(int evt_len, char *detMapFileName, char *filterFilename, char *trGainFilename, char *xTalkParsFilename);
int preProcess(unsigned short int *sbuf, int inlentotal, int evt_len, Event_Signal *evt, int new, preprocCnt *diagcnt);

int read3nhdr(char *filename, int holenum, int xtal, unsigned short int *obuf);
int get_evt_len(FILE *fin);

pair *read_mapfile(char *filename, float *a0, float *a1);
struct filter *read_filterfile(char *filename);
int read_param(const char *filename, const char *label, float *x, int len);

void ReadRawINL();
int inl_id(unsigned short int *x);
double ComputeChiSquare(int startIndex, int nPoints, double tau, double Nstart, double baseline);
int LSFitExponential(int startIndex, int nPoints, double Nstart, double tau, double *NstartFit, double *baselineFit, double *NstartFitError, double *baselineFitError);
void calculateENL(int xtal, int chn, float base, int s1_width, int s2_width, float s1, float s2);

/* decomp routines */
struct decomp_state *dl_decomp_init(char *basis_file_name, int set_quiet);
struct crys_intpts *dl_decomp(struct decomp_state *di, Event_Signal *asig, postprocCnt *postCnt);
int read_basis(char *basisfile);
int grid_init(void);
char *dl_crys_intpts_2s(struct crys_intpts *x);
int cyl_to_cart(int seg, double *pars, double*x, double *y, double*z, double *e);
double fitter(const Event_Signal *asig, Event_Signal *bsig, Interaction *ints, double *t0, int nints, int final);
int eval(const Event_Signal *asig, Event_Signal *bsig, Interaction *ints, double t0, int nints, double *chisq_out, double beta[MAX_PARS], double alpha[MAX_PARS][MAX_PARS], int calc_deriv, int *ssel);
double eval_int_pos(const Event_Signal *asig, Event_Signal *bsig, Interaction *ints, double t0, int nints);
int matinv(double *array, int norder, int dim);
int decompose_1(const Event_Signal *asig, Event_Signal *bsig, int seg, Interaction *ints, double *t0, double *chisq_out, int grid2, int fit0, int fit1, int fit2, int fit3, int fit4, int final_fit, int coalesce, double min_e_fraction);
int decompose_n(const Event_Signal *asig, Event_Signal *bsig, int nseg, int *seg, int coalesce, Interaction *ints, double *t0, double *chisq_out);
double coarse_grid_1(const Event_Signal *asig, int seg, Interaction *ints, double *chisq0, double min_e_fraction);
double refine_grid_1(const Event_Signal *asig, double chisq, double chisq0, double min_e_fraction, Interaction *ints);
int interpolate(int seg, double rin, double pin, double zin, Basis_Point *signal, Basis_Point sigderiv[3], int calc_deriv, int *ssel);

/* fal layer */
int ener_nd(unsigned short int *e);
long long int time_nd(unsigned short int *e);
int vec_nd(unsigned short int *e);
int ch_nd(unsigned short int *e);
int tr_nd(int *tr, int tr_len, int, unsigned short int *e);
int tr_len_nd(unsigned short int *e);
int evt_len_nd(unsigned short int *e);

/* util routines */
int segevt_id(uint2 *e);
int segevt_ener(uint2 *e, int shift);
void adjoff(int *tr, int tr_len);
int net(int *tr, int tr_len);

/* time alignment routines */
int time_sef(int *tr, int len);
float t_cfd(int *buf, int buf_len);
int align_cfd_1(int *traces, int tr_len, float *delay1, float *delay0);

/* logging routines - not implemented  - from Carl's logSend */
int sendToLog(char *msg);
int sendLogMsg(char *from, char *what);

/* mode 2 logging routines */
void log2intpts(struct crys_intpts *a);
void logintpt(float x, float y, float z, float e, FILE *flog);
void logMsg(int verboseFlag, const char *fmt, ...);
void errMsg(const char *fmt, ...);

/* routines to support one interaction search */
void numGridPts(int stride);
Basis_Point *oneIntSearch0(const Event_Signal *asig, Event_Signal *bsig, int seg, double *chisq);
int oneIntSearch1(const Event_Signal *asig, Event_Signal *bsig, int numNet, int *netSeg, simpleInteraction *ints, double *chisq);
int oneIntSearch2(const Event_Signal *asig, Event_Signal *bsig, int numNet, int *netSeg, simpleInteraction *ints, double *chisq);
int oneIntSearchIter(const Event_Signal *asig, Event_Signal *bsig, int numNet, int *netSeg, simpleInteraction *ints, double *chisq);
void scaleToOne(const Event_Signal *asig, Event_Signal *bsig, int netSeg);
void subOverlapImage(const Event_Signal *asig, Event_Signal *xsig, int overlap, int net, int othernet, Basis_Point *otherpoint);
void asigDiag(const Event_Signal *asig, int seg);
int normDiag(const Event_Signal *asig, int nseg, int *seg);
int badTrPred(const Event_Signal *asig, int nseg, int *seg);
