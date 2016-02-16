/*
Feb 5 2016
Code to make signal from given bases.
Optional function
* select interested segmant
* shift signal with time
* add random noise

Contact: Ryo Taniuchi taniuchi@nucl.phys.s.u-tokyo.ac.jp / rtaniuchi@lbl.gov
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <time.h>
#include "petcat.h"
#include "vegcat.h"
#include "neigh.h"

int i, j, k, ii, dummy;
int nseg = 15, tOffset=140, tShiftFlag = 0, NoiseFlag = 0, numevt =1;
char *outfilename = "basis/seg15.dat";
double Shift = 0., Noise = 0., Rand;
Mario *evts, *oneeve;
FILE *fgeo, *fdat;
struct gebData defaulthdr = {100, sizeof(Mario), 0ll};

double genrand(unsigned int argint){
  srand((unsigned)time(NULL) + argint );
  //printf("%ld",(unsigned)time(NULL) * argint);
  //return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
  return ((double)(rand()%10000)+1.0)/(10000.+2.0);
}
	     

int main(){
  evts = malloc (numevt * sizeof(Mario));
  oneeve = malloc ( sizeof(Mario));

  fgeo = fopen("inputs/geo_basis.txt", "w");
  fprintf(fgeo, "runn, x, y, z \n");
  fdat = fopen(outfilename,"w");
  
  struct decomp_state *a;
  a = dl_decomp_init(cfg.basisName, 1); // 1 to suppress diag info
  if (a == 0) { fprintf(stderr, "decomp init failed!\n"); exit(1); }
  
  for(i=0; i<GRID_PTS; i++){
    //printf("%i %i %f \n",i, basis[i].iseg,  basis[i].signal[0][0]);
    if(basis[i].iseg == nseg){
      //Generate signal
      int tshift = tOffset + (int)(Shift * (genrand(i)-0.5) );
      Rand = genrand( j * i );
      for(j=0; j<numevt; j++){
	for(k=0; k<tshift; k++){
	  for(ii=0; ii<37; ii++){
	    oneeve->wf[ii][k] = 0;
	  }
	}
	for(k=tshift; k<tshift+50; k++){
	  for(ii=0; ii<37; ii++){
	    oneeve->wf[ii][k] = (short)(Rand * 10000.*basis[i].signal[ii][k - tshift] + 0.5); //Float (0, 1) to short(Assuming 14bit digitizer)
	    }	  
	  }
	}
	for(k=tshift+50; k<300; k++){
	  for(ii=0; ii<37; ii++){
	    oneeve->wf[ii][k] =(short)(Rand *10000. * (basis[i].signal[ii][49] * exp(-((double)k - tshift- 49.)/800.)
						      //  +(basis[i].signal[ii][49] - basis[i].signal[ii][48]) * (k - 49 - tshift) 
						      ) +0.5);
	  }
	}

	for (k = 0; k < 36; k++) {
	  oneeve->segEnergy[k] = 0.;
	}
	oneeve->segEnergy[nseg] = 10000.;
	oneeve->ccEnergy = 10000.;
	
	
	//Write output information
	fprintf(fgeo, "%i, %f, %f, %f \n", i, basis[i].x, basis[i].y, basis[i].z);
	//output dat
	assert( 1 == fwrite(&defaulthdr, sizeof(struct gebData), 1, fdat));
	assert( 1 == fwrite(oneeve, sizeof(Mario), 1, fdat));
	//set energy
	//break;
    }
  }
  fclose(fgeo);
  fclose(fdat);



  /*
    struct gebData ghdr, defaulthdr = {100, sizeof(Mario), 0ll};
    assert( 1 == fwrite(&defaulthdr, sizeof(struct gebData), 1, fou));
    assert( 1 == fwrite(pList[i].rawEvts, sizeof(Mario), 1, fou));
  */


  return 0;

}
