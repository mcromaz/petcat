#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TStyle.h"
using namespace std;

const int NUMRUNS=20000;
char *filedir="Mar14";
char name[] ="basis_noise";
char *index_xyz[4] = {"x", "y", "z", "pos"};
string dummystring;
int runn[NUMRUNS], runn_poslist[NUMRUNS], segn[NUMRUNS],  dumint;
double par[2][NUMRUNS], parerr[2][NUMRUNS], x[NUMRUNS], y[NUMRUNS], z[NUMRUNS], r[NUMRUNS], x_list[NUMRUNS], y_list[NUMRUNS], z_list[NUMRUNS], x_meas[NUMRUNS], y_meas[NUMRUNS], z_meas[NUMRUNS], deltapos[NUMRUNS], efrac[NUMRUNS], dumdouble;

TH1F *hx[5], *hy[5], *hz[5], *hpos[36][4][6], *hpos2hit[3], *hratio2hit;
TH2F *hxy;

int i, j, k, numpos, draw = 0;
TCanvas *c, *cxyz, *cdist;//, *c_gr;
TGraph *g[4]; //Delta x, y, z, pos
TLegend *leg, *leg2hit;

ofstream table;
int num1mm, num2mm, num5mm, num2hits;

int positiondist_basis_id(int segid);


void positiondist_basis(void) {
  cxyz = new TCanvas("cxyz","cxyz",1050,800);
  cxyz->Print(Form("./%s/positiondist/%s_xyzdist.pdf[",filedir, name));
  cxyz -> Divide(2,2);
  cdist = new TCanvas("cdist","cdist",1050,800);
  cdist->Print(Form("./%s/positiondist/%s_dist.pdf[",filedir, name));
  cdist ->Divide(3,2);
  leg = new TLegend(0.82,0.7,0.92,0.9);

  for(i=0; i<3; i++){
    hpos2hit[i] = new TH1F(Form("hpos2hit%i",i), "position distribution of weighed average position for 2-hits-decomposition",100,0,100);
    hpos2hit[i]->SetLineColor(i+1);
    hpos2hit[i]->GetXaxis()->SetTitle("deviation / mm");
    hpos2hit[i]->GetYaxis()->SetTitle("counts / mm");
  }
  hratio2hit = new TH1F("hratio2hit", "ratio of the fraction energy for 2-hits-decomposition",100,0,10);
  hratio2hit->GetXaxis()->SetTitle("deviation / mm");
  hratio2hit->GetYaxis()->SetTitle("counts / mm");

  hxy = new TH2F("xydist", "xy distribution for 2 hits-decomposition", 80, -40, 40, 80, -40 ,40);
  hxy->GetXaxis()->SetTitle("x / mm");
  hxy->GetYaxis()->SetTitle("y / mm");


  table.open(Form("./%s/table/%s.csv",filedir, name),ios::out);
  table<<"segid, numpos, num1mm, num2mm, num5mm, num2hits, ";
  table<<"ratio numpos, ratio num1mm, ratio num2mm, ratio num5mm, ratio num2hits"<<endl;

  for(int segid=0; segid<36; segid++) {
    //for(int segid=0; segid<1; segid++) {
    num1mm=0;
    num2mm=0;
    num5mm=0;
    num2hits=0;
    positiondist_basis_id(segid);
    table<<segid<<", "<<numpos<<", "<<num1mm<<", "<<num2mm<<", "<<num5mm<<", "<<num2hits;
    table<<", "<<(double)numpos/(double)numpos *100.<<", "<<(double)num1mm/(double)numpos *100.<<", "<<(double)num2mm/(double)numpos *100.<<", "<<(double)num5mm/(double)numpos *100.<<", "<<(double)num2hits/(double)numpos *100.<<endl;

    delete cxyz;
    cxyz = new TCanvas("cxyz","cxyz",1050,800);
    cxyz -> Divide(2,2);
    if(segid%6 == 5){ 
      cdist->cd(0);
      leg->SetFillColor(0);
      leg->Draw("");
      cdist->Print(Form("./%s/positiondist/%s_dist.pdf",filedir, name));
      delete cdist;
      cdist = new TCanvas("cdist","cdist",1050,800);
      cdist ->Divide(3,2);
      delete leg;
      leg = new TLegend(0.82,0.7,0.92,0.9);
    }
  }
  table.close();
  cdist->cd(0);
  gStyle -> SetOptStat(110110);
  hxy->Draw("colz");
  cdist->Print(Form("./%s/positiondist/%s_dist.pdf",filedir, name));
  hratio2hit->Draw("");
  cdist->Print(Form("./%s/positiondist/%s_dist.pdf",filedir, name));

  gStyle -> SetOptStat(0);
  for(i=1; i<3; i++){
    if(hpos2hit[0]->GetMaximum()<hpos2hit[i]->GetMaximum()*1.1)
      hpos2hit[0]->SetMaximum(hpos2hit[i]->GetMaximum()*1.1);
  }
  hpos2hit[0]->Draw("");
  hpos2hit[1]->Draw("same");
  hpos2hit[2]->Draw("same");
  leg2hit = new TLegend(0.7,0.7,0.92,0.9);
  leg2hit->SetFillColor(0);
  leg2hit->AddEntry(hpos2hit[0],"Weighed","l");
  leg2hit->AddEntry(hpos2hit[1],"Higher energy","l");
  leg2hit->AddEntry(hpos2hit[2],"Lower energy","l");
  leg2hit->Draw("same");
  cdist->Print(Form("./%s/positiondist/%s_dist.pdf",filedir, name));
 
  cxyz->Print(Form("./%s/positiondist/%s_xyzdist.pdf]",filedir, name));
  cdist->Print(Form("./%s/positiondist/%s_dist.pdf]",filedir, name));
}


int positiondist_basis_id(int segid) {

  ifstream poslist(Form("./inputs/geo_basis%i.txt",segid)); //in principle same
  ifstream decoposlist(Form("./%s/%s%i_out.csv",filedir, name,segid));
  
  char dummy, outname[50];
  sprintf(outname, Form("./%s/positiondist/%s%i",filedir, name, segid));


  if(!poslist||!decoposlist){
    cerr<<"file is not found"<<endl;
    return 1;
  }
  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(111);
  c = new TCanvas("c", "Deviation of measured positions", 900, 600);

  for(i=0; i<6; i++){
    hpos[segid][0][i] = new TH1F(Form("hx%i_seg%i",i,segid),Form("hx_%i to %i",5*i+5,5*i+5+5),100,-10,10);
    hpos[segid][1][i] = new TH1F(Form("hy%i_seg%i",i,segid),Form("hy_%i to %i",5*i+5,5*i+5+5),100,-10,10);
    hpos[segid][2][i] = new TH1F(Form("hz%i_seg%i",i,segid),Form("hz_%i to %i",5*i+5,5*i+5+5),100,-10,10);
    hpos[segid][3][i] = new TH1F(Form("hpos%i_seg%i",i,segid),Form("hpos_%i to %i",5*i+5,5*i+5+5),100,0,20);
    for(j=0; j<4; j++){
      hpos[segid][j][i]->GetXaxis()->SetTitle("deviation / mm");
      hpos[segid][j][i]->GetYaxis()->SetTitle("counts / mm");
    }
  }
  

  ////////Get Pos list
  i = 0;
  getline(poslist, dummystring); //runn, x, y, z
  while(true){
    poslist >> runn_poslist[i] >> dummy >> x_list[i] >> dummy >> y_list[i] >> dummy >> z_list[i];
    if(i == NUMRUNS-1 || !poslist) {
      break; 
    }
    i++;
  }
  numpos = i;
  //////////////////////

  getline(decoposlist, dummystring); //runnum, segnum, evtnum, avex, avey, avez, aveE
  i=0;
  k=0;
  while(true){
    decoposlist >> runn[k] >> dummy >> dumint >> dummy >> dumint >>dummy >> segn[k] >> dummy >> x_meas[k] >> dummy >> y_meas[k] >> dummy >> z_meas[k] >> dummy >> efrac[k] >> dummy >> dumdouble >> dummy >> dumdouble >> dummy >> dumdouble >> dummy >> dumdouble >> dummy >> dumdouble;
    if(dumint!=1) {
      if(dumint==2) {
	hxy->Fill(x_meas[k],y_meas[k]);
	hxy->Fill(x_meas[k-1],y_meas[k-1]);
	hratio2hit->Fill(efrac[k]>efrac[k-1]?efrac[k]/efrac[k-1]:efrac[k-1]/efrac[k]);
	//
	x[k] = x_list[i-1];
	y[k] = y_list[i-1];
	z[k] = z_list[i-1];
	deltapos[k-1] = sqrt(pow(x_meas[k-1]-x[k], 2) +  pow(y_meas[k-1]-y[k], 2) +  pow(z_meas[k-1]-z[k], 2));
	deltapos[k] = sqrt(pow(x_meas[k]-x[k], 2) +  pow(y_meas[k]-y[k], 2) +  pow(z_meas[k]-z[k], 2));
	hpos2hit[1] ->Fill(efrac[k]>efrac[k-1]?deltapos[k]:deltapos[k-1]);
	hpos2hit[2] ->Fill(efrac[k]>efrac[k-1]?deltapos[k-1]:deltapos[k]);

	x_meas[k] = (efrac[k]*x_meas[k] + efrac[k-1]*x_meas[k-1])/(efrac[k-1]+efrac[k]);
	y_meas[k] = (efrac[k]*y_meas[k] + efrac[k-1]*y_meas[k-1])/(efrac[k-1]+efrac[k]);
	z_meas[k] = (efrac[k]*z_meas[k] + efrac[k-1]*z_meas[k-1])/(efrac[k-1]+efrac[k]);
	x_meas[k] -= x[k]; //Delta x
	y_meas[k] -= y[k];
	z_meas[k] -= z[k];
	deltapos[k] = sqrt(pow(x_meas[k], 2) +  pow(y_meas[k], 2) +  pow(z_meas[k], 2));
	hpos2hit[0]->Fill(deltapos[k]);
	//
	k--; // Overwrite the previous measured values with several intaractions
	num2hits++;
      }else{
	cerr<<"more than 2 events error"<<endl;
	return 1;
      }
      continue;
    }
    if (!decoposlist) {
      cout<<"End of decomp list of seg"<<segid<<": numpos ="<<numpos<<", numdecomp ="<<k<<endl;
      break;
    }
    x[k] = x_list[i];
    y[k] = y_list[i];
    z[k] = z_list[i];
    r[k] = sqrt(x[k]*x[k]+y[k]*y[k]);
    x_meas[k] -= x[k]; //Delta x
    y_meas[k] -= y[k];
    z_meas[k] -= z[k];
    deltapos[k] = sqrt(pow(x_meas[k], 2) +  pow(y_meas[k], 2) +  pow(z_meas[k], 2));
    if(deltapos[k]<1){
      num1mm++;
      num2mm++;
      num5mm++;
    }else if(deltapos[k]<2.){
      num2mm++;
      num5mm++;
    }else if(deltapos[k]<5.){
      num5mm++;
    }
    for(j=0;j<6;j++){
      if(j==5 || (r[k]>5.*((double)j+1.)&&r[k]<=5.*((double)j+2.))){
	hpos[segid][0][j]->Fill(x_meas[k]);
	hpos[segid][1][j]->Fill(y_meas[k]);
	hpos[segid][2][j]->Fill(z_meas[k]);
	hpos[segid][3][j]->Fill(deltapos[k]);
	break;
      }
    }
    i++;
    k++;
  }
  cout<<endl<<" Numpos = "<<numpos<<" num 1hit="<<k<<endl<<endl;

  ////////////////////////  
  g[0]= new TGraph(k, r, x_meas);
  g[1]= new TGraph(k, r, y_meas);
  g[2]= new TGraph(k, r, z_meas);
  g[3]= new TGraph(k, r, deltapos);

  if(draw) c->Print(Form("%s.pdf[",outname));
  for (i=0; i<4; i++){
    g[i] -> SetTitle(Form("#Delta%s vs Radius",index_xyz[i]));
    g[i] -> GetXaxis() -> SetTitle("Radius / mm");
    g[i] -> GetYaxis() -> SetTitle(Form("#Delta%s / mm",index_xyz[i]));
    g[i] -> GetXaxis() -> SetLimits(0,40.);
    g[i] -> SetMinimum(-30.);
    g[i] -> SetMaximum(30.);
    g[i] -> Draw("AP*");
    if(draw) c->Print(Form("%s.pdf",outname));
  }
  delete c;
  c = new TCanvas("c", "Deviation of measured positions hist", 900, 600);
  c->Divide(3,2);
  gStyle -> SetOptFit(111);
   
  for (i=0; i<4; i++){
    for (j=0; j<6; j++){
      c->cd(j+1);
      hpos[segid][i][j]->Draw("");
      if(i!=3) hpos[segid][i][j]->Fit("gaus","Q");
    }
    if(draw) c->Print(Form("%s.pdf",outname));
  }
  
  if(draw) c->Print(Form("%s.pdf]",outname));

  delete c;

  for(i=0; i<4; i++){
    cxyz->cd(i+1);
    g[i]->Draw("AP*");
  }
  i=3;
  for (j=0; j<6; j++){
    cdist->cd(j+1);
    gPad->SetLogy(1);
    hpos[segid][i][j]->SetLineColor(segid%6+1);
    if (segid%6==0){
      hpos[segid][i][j]->Draw("");
    }else hpos[segid][i][j]->Draw("same");
  }
  leg->AddEntry(hpos[segid][i][0],Form("seg%i",segid),"l");
    
  //if(draw) 
  cxyz->Print(Form("./%s/positiondist/%s_xyzdist.pdf",filedir, name));
   

  for (i=0; i<4; i++){
    delete g[i];
  }
  return 0;
}


