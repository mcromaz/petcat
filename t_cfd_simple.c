#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "petcat.h"
//#include "vegcat.h"

//float cft(Mario *evt, int seg, double cf, int cfdelay, int threshold){ // cf: constant fraction, i:channel(normally 36) cfdelay: delay channnels, threshold: should be negative

typedef struct {
  float ccEnergy;
  float segEnergy[36];
  float pad;
  short wf[37][300];
} Mario;


float t_cfd_simple(Mario *evt, int seg, double cf, int cfdelay, int threshold){
  double sumwf[300], zerocross;
  float s0, s1;
  int j, k, thflg = 0;

  for (j = 0; j<cfdelay; j++) {
    s1 = cf * (double)evt->wf[seg][j];
    sumwf[j] = s1;
  }
  for (j = cfdelay; j < 300; j++) {
    s0 = (double)evt->wf[seg][j - cfdelay];
    s1 = cf * (double)evt->wf[seg][j];
    sumwf[j] = s0 + s1;
  }

  for(j= 10; j < 300 - 10; j++) {
    if(thflg == 1) {
      s0 = 0.;
      s1 = 0.;
      for(k = 0; k < 5; k++){
	s0 += sumwf[j - k]/5.;
	s1 += sumwf[j + 1 + k]/5.;
      }
      if(s0 < 0. && s1 > 0. && s1 - s0 > sumwf[j+1] - sumwf[j])  break;
    }else if(thflg == 0 && sumwf[j] < threshold){
      thflg = 1;
    }
  }
  //zerocross =  (double) j + sumwf[j] / (sumwf[j]-sumwf[j+1]); 
  zerocross =  (double) j -2. + 5. * s0 / (s0 - s1); 

  if(thflg == 1){
    return zerocross;
  } else {
    return -1000.;
    }
}


float cft(Mario *evt) {
  return t_cfd_simple(evt, 36, -0.5, 10, -50);
}
