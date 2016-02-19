#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;


const int NUMRUNS=6000;
char name[] ="basis";
ifstream poslist("./inputs/geo_basis.txt");
//ifstream aveposlist(Form("./vegcatout/%s.txt_avepos.csv",name));
ifstream aveposlist(Form("./vegcatout/%s.txt_out.csv",name));

char dummy, outname[] = Form("./vegcatout/positiondist/%s",name);
char *index[4] = {"x", "y", "z", "pos"};
string dummystring;
int runn[NUMRUNS], runn_poslist[NUMRUNS], segn[NUMRUNS],  dumint;
double par[2][NUMRUNS], parerr[2][NUMRUNS], x[NUMRUNS], y[NUMRUNS], z[NUMRUNS], r[NUMRUNS], x_list[NUMRUNS], y_list[NUMRUNS], z_list[NUMRUNS], x_meas[NUMRUNS], y_meas[NUMRUNS], z_meas[NUMRUNS], deltapos[NUMRUNS], dumdouble;

TH1F *hx[5], *hy[5], *hz[5], *hpos[4][6];

int i, j, k, numpos, draw = 1;
TCanvas *c;//, *c_gr;
TGraph *g[4]; //Delta x, y, z, pos
void positiondist_basis(void) {

  if(!poslist||!aveposlist){
    cerr<<"file is not found"<<endl;
    break;
  }
  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(111);
  c = new TCanvas("c", "Deviation of measured positions", 900, 600);

  for(i=0; i<6; i++){
    hpos[0][i] = new TH1F(Form("hx%i",i),Form("hx_%ito%i",5*i+5,5*i+5+5),100,-10,10);
    hpos[1][i] = new TH1F(Form("hy%i",i),Form("hy_%ito%i",5*i+5,5*i+5+5),100,-10,10);
    hpos[2][i] = new TH1F(Form("hz%i",i),Form("hz_%ito%i",5*i+5,5*i+5+5),100,-10,10);
    hpos[3][i] = new TH1F(Form("hpos%i",i),Form("hpos_%ito%i",5*i+5,5*i+5+5),100,0,20);
  }
  

  ////////////////////// Get poslist
  i = 0;
  getline(poslist, dummystring); //runn, x, y, z
  while(true){
    poslist >> runn_poslist[i] >> dummy >> x_list[i] >> dummy >> y_list[i] >> dummy >> z_list[i];
    if(i == NUMRUNS-1 || !poslist) break; 
    i++;
  }
  numpos = i;
  //////////////////////

  getline(aveposlist, dummystring); //runnum, segnum, evtnum, avex, avey, avez, aveE
  //cout<<dummystring<<endl;
  for(k = 0; k < NUMRUNS;  k++){
    // for(k = 0; k < 1;  k++){
    //aveposlist >> runn[k] >> dummy >> segn[k] >> dummy >>dumint >> dummy >> x_meas[k] >> dummy >> y_meas[k] >> dummy >> z_meas[k] >> dummy >> dumdouble;
    aveposlist >> runn[k] >> dummy >> dumint >> dummy >> dumint >>dummy >> segn[k] >> dummy >> x_meas[k] >> dummy >> y_meas[k] >> dummy >> z_meas[k] >> dummy >> dumdouble >> dummy >> dumdouble;
    //cout << runn[k] << dummy << segn[k] << dummy <<dumint << dummy << x_meas[k] << dummy << y_meas[k] << dummy << z_meas[k] << dummy << dumdouble <<endl;
    if(dumint!=1) aveposlist >> runn[k] >> dummy >> dumint >> dummy >> dumint >>dummy >> segn[k] >> dummy >> x_meas[k] >> dummy >> y_meas[k] >> dummy >> z_meas[k] >> dummy >> dumdouble >> dummy >> dumdouble;
    if (!aveposlist) break;
    //cout << runn[k] << dummy << dumint << dummy << dumint <<dummy << segn[k] << dummy << x_meas[k] << dummy << y_meas[k] << dummy << z_meas[k] << dummy << dumdouble<<endl;

    /// serch matched geometrical information from poslist ///
    /*for (i = 0; i < numpos + 1; i++){
      if(i == numpos){
	cout<<"no position information for run"<<runn[k]<<" is not found"<<endl;
	x[k] = 0;
	y[k] = 0;
	z[k] = 0;
	r[k] = 0;
	break;
	}else if(runn[k] == runn_poslist[i]){*/
    i=k;
	x[k] = x_list[i];
	y[k] = y_list[i];
	z[k] = z_list[i];
	r[k] = sqrt(x[k]*x[k]+y[k]*y[k]);
	x_meas[k] -= x[k]; //Delta x
	y_meas[k] -= y[k];
	z_meas[k] -= z[k];
	deltapos[k] = sqrt(pow(x_meas[k], 2) +  pow(y_meas[k], 2) +  pow(z_meas[k], 2));
	cout<<k<<" "<<r[k]<<" "<<deltapos[k]<<endl;
	
	for(i=0;i<6;i++){
	  if(r[k]>5.*((double)i+1.)&&r[k]<=5.*((double)i+2.)){
	    hpos[0][i]->Fill(x_meas[k]);
	    hpos[1][i]->Fill(y_meas[k]);
	    hpos[2][i]->Fill(z_meas[k]);
	    hpos[3][i]->Fill(deltapos[k]);
	    break;
	  }
	  else if(i==5){
	    cout<<"Error: this point is out of r=35"<<endl;
	  }
	}

	/*break;
      }
      }*/
  }

  cout<<endl<<" Numpos = "<<numpos<<" "<<k<<endl;

  //  for (i=0; i<4; i++){
  g[0]= new TGraph(numpos, r, x_meas);
  g[1]= new TGraph(numpos, r, y_meas);
  g[2]= new TGraph(numpos, r, z_meas);
  g[3]= new TGraph(numpos, r, deltapos);
  /*g[0]= new TGraph(NUMRUNS, r, x_meas);
  g[1]= new TGraph(NUMRUNS, r, y_meas);
  g[2]= new TGraph(NUMRUNS, r, z_meas);
  g[3]= new TGraph(NUMRUNS, r, deltapos);
*/
  if(draw) c->Print(Form("%s.pdf[",outname));
  for (i=0; i<4; i++){
    //h[i] = new TH2F(Form("h%i", i), Form("#Delta%s vs Radius",index[i]), 100, 0, 30, );
    g[i] -> SetTitle(Form("#Delta%s vs Radius",index[i]));
    g[i] -> GetXaxis() -> SetTitle("Radius / mm");
    g[i] -> GetYaxis() -> SetTitle(Form("#Delta%s / mm",index[i]));
    g[i] -> GetXaxis() -> SetLimits(0,35.);
    //g[i] -> SetMarkerStyle(21);
    g[i] -> Draw("AP*");
    if(draw) c->Print(Form("%s.pdf",outname));
  }
  //c->Clear();
  delete c;
  c = new TCanvas("c", "Deviation of measured positions hist", 900, 600);
  c->Divide(3,2);
  gStyle -> SetOptFit(111);
   
  for (i=0; i<4; i++){
    for (j=0; j<6; j++){
      c->cd(j+1);
      hpos[i][j]->Draw("");
      if(i!=3) hpos[i][j]->Fit("gaus","Q");
    }
    if(draw) c->Print(Form("%s.pdf",outname));
  }
  
  if(draw) c->Print(Form("%s.pdf]",outname));

}


