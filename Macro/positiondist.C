#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;

const int NUMRUNS = 17;

ifstream aveposlist("./vegcatout/runlist.txt_avepos.csv");
ifstream poslist("./inputs/geometry.txt");
char dummy, outname[] = "./vegcatout/positiondist";
char *index[4] = {"x", "y", "z", "pos"};
string dummystring;
int runn[NUMRUNS], runn_poslist[NUMRUNS], segn[NUMRUNS],  dumint;
double par[2][NUMRUNS], parerr[2][NUMRUNS], x[NUMRUNS], y[NUMRUNS], z[NUMRUNS], r[NUMRUNS], x_list[NUMRUNS], y_list[NUMRUNS], z_list[NUMRUNS], x_meas[NUMRUNS], y_meas[NUMRUNS], z_meas[NUMRUNS], deltapos[NUMRUNS], dumdouble;

int i, j, k, numpos, draw = 1;
TCanvas *c;//, *c_gr;
TGraph *g[4]; //Delta x, y, z, pos
void positiondist(void) {

  if(!poslist||!aveposlist){
    cerr<<"file is not found"<<endl;
    break;
  }
  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(111);
  c = new TCanvas("c", "Deviation of measured positions", 900, 600);

  ////////////////////// Get poslist
  i = 0;
  getline(poslist, dummystring); //runn, x, y, z
  while(true){
    poslist >> runn_poslist[i] >> dummy >> x_list[i] >> dummy >> y_list[i] >> dummy >> z_list[i];
    if(i == 16 || !poslist) break; 
    i++;
  }
  numpos = i;
  //////////////////////

  getline(aveposlist, dummystring); //runnum, segnum, evtnum, avex, avey, avez, aveE
  //cout<<dummystring<<endl;
  for(k = 0; k < NUMRUNS;  k++){
    // for(k = 0; k < 1;  k++){
    aveposlist >> runn[k] >> dummy >> segn[k] >> dummy >>dumint >> dummy >> x_meas[k] >> dummy >> y_meas[k] >> dummy >> z_meas[k] >> dummy >> dumdouble;
    //cout << runn[k] << dummy << segn[k] << dummy <<dumint << dummy << x_meas[k] << dummy << y_meas[k] << dummy << z_meas[k] << dummy << dumdouble <<endl;
    if (!aveposlist) break;

    /// serch matched geometrical information from poslist ///
    for (i = 0; i < numpos + 1; i++){
      if(i == numpos){
	cout<<"no position information for run"<<runn[k]<<" is not found"<<endl;
	x[k] = 0;
	y[k] = 0;
	z[k] = 0;
	r[k] = 0;
	break;
      }else if(runn[k] == runn_poslist[i]){
	x[k] = x_list[i];
	y[k] = y_list[i];
	z[k] = z_list[i];
	r[k] = sqrt(x[k]*x[k]+y[k]*y[k]);
	x_meas[k] -= x[k]; //Delta x
	y_meas[k] -= y[k];
	z_meas[k] -= z[k];
	deltapos[k] = sqrt(pow(x_meas[k], 2) +  pow(y_meas[k], 2) +  pow(z_meas[k], 2));
	break;
      }
    }
  }

  //  for (i=0; i<4; i++){
  g[0]= new TGraph(NUMRUNS, r, x_meas);
  g[1]= new TGraph(NUMRUNS, r, y_meas);
  g[2]= new TGraph(NUMRUNS, r, z_meas);
  g[3]= new TGraph(NUMRUNS, r, deltapos);

  if(draw) c->Print(Form("%s/posdist.pdf[",outname));
  for (i=0; i<4; i++){
    //h[i] = new TH2F(Form("h%i", i), Form("#Delta%s vs Radius",index[i]), 100, 0, 30, );
    g[i] -> SetTitle(Form("#Delta%s vs Radius",index[i]));
    g[i] -> GetXaxis() -> SetTitle("Radius / mm");
    g[i] -> GetYaxis() -> SetTitle(Form("#Delta%s / mm",index[i]));
    g[i] -> GetXaxis() -> SetLimits(0,30.);
    //g[i] -> SetMarkerStyle(21);
    g[i] -> Draw("AP*");
    if(draw) c->Print(Form("%s/posdist.pdf",outname));
  }
  if(draw) c->Print(Form("%s/posdist.pdf]",outname));

}


