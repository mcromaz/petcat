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
const int seg =15;
char *filedir="Mar14";
char name[] ="basis_2sig15";
char *index_xyz[4] = {"x", "y", "z", "pos"};
string dummystring;
int runn[2][NUMRUNS], runn_poslist[2][NUMRUNS], segn[2][NUMRUNS],  dumint, ratio[6]={1,2,4,6,8,9};
double par[2][NUMRUNS], parerr[2][NUMRUNS], x[2][NUMRUNS], y[2][NUMRUNS], z[2][NUMRUNS], r[2][NUMRUNS], x_list[2][NUMRUNS], y_list[2][NUMRUNS], z_list[2][NUMRUNS], x_meas[2][NUMRUNS], y_meas[2][NUMRUNS], z_meas[2][NUMRUNS], deltapos[2][NUMRUNS], efrac[2][NUMRUNS], posdist[NUMRUNS], dumdouble;

TH1F *hx[5], *hy[5], *hz[5], *hpos[36][4][6], *hpos2hit[3], *hratio2hit[6], *h_effdist[6], *h_totevts[6];
TH2F *hxy[6];

int i, j, k, numpos, draw = 0, invertflag=0;
TCanvas *c, *cxyz, *cdist;//, *c_gr;
TGraph *g[4]; //Delta x, y, z, pos
TLegend *leg, *leg2hit;

ofstream table;
int num1mm, num2mm, num5mm, num2hits;

int positiondist_basis_id(int ratioindex);


void positiondist_2sig(void) {
  cxyz = new TCanvas("cxyz","cxyz",1050,800);
  cxyz->Print(Form("./%s/positiondist/%s_xyzdist.pdf[",filedir, name));
  cxyz -> Divide(2,2);
  cdist = new TCanvas("cdist","cdist",1050,800);
  cdist->Print(Form("./%s/positiondist/%s_dist.pdf[",filedir, name));
  cdist ->Divide(3,2);
  leg = new TLegend(0.82,0.7,0.92,0.9);
  /*
  for(i=0; i<3; i++){
    hpos2hit[i] = new TH1F(Form("hpos2hit%i",i), "position distribution of weighed average position ",100,0,100);
    hpos2hit[i]->SetLineColor(i+1);
    hpos2hit[i]->GetXaxis()->SetTitle("deviation / mm");
    hpos2hit[i]->GetYaxis()->SetTitle("counts / mm");
  }*/
  for(i=0; i<6; i++){
    hratio2hit[i] = new TH1F(Form("hratio2hit%i",ratio[i]), Form("ratio of the fraction energy %i:%i (%2.1f)",ratio[i],10-ratio[i],ratio[i]<5?(float)ratio[i]/10.:1.-(float)ratio[i]/10.),200,0,.5);
    hratio2hit[i]->GetXaxis()->SetTitle("ratio");
    hratio2hit[i]->GetYaxis()->SetTitle("counts / a.u.");
    hxy[i] = new TH2F(Form("xydist%i",ratio[i]), Form("xy distribution %i:%i", ratio[i],10-ratio[i]), 400, -40, 40, 400, -40 ,40);
    hxy[i]->GetXaxis()->SetTitle("x / mm");
    hxy[i]->GetYaxis()->SetTitle("y / mm");
    h_totevts[i] = new TH1F(Form("totevts%i",ratio[i]), Form("totevts%i",ratio[i]), 35, 0, 35);
    h_totevts[i]->GetXaxis()->SetTitle("distance of two points / mm");
    h_totevts[i]->GetYaxis()->SetTitle("total counts / mm");
    h_effdist[i] = new TH1F(Form("2-hits decomposition efficiencies for %i:%i", ratio[i],10-ratio[i]), Form("2 hit decomposition efficiencies for %i:%i", ratio[i],10-ratio[i]), 35, 0, 35);
    h_effdist[i]->GetXaxis()->SetTitle("distance of two points / mm");
    h_effdist[i]->GetYaxis()->SetTitle("total counts / mm");
  }

  table.open(Form("./%s/table/%s.csv",filedir, name),ios::out);
  table<<"ratioindex, numpos, num1mm, num2mm, num5mm, num2hits, ";
  table<<"ratio numpos, ratio num1mm, ratio num2mm, ratio num5mm, ratio num2hits"<<endl;

  //for(int ratioindex=0; ratioindex<36; ratioindex++) {
  for(int ratioindex=0; ratioindex<6; ratioindex++) { //ratio[9]
    num1mm=0;
    num2mm=0;
    num5mm=0;
    num2hits=0;
    positiondist_basis_id(ratioindex);
    table<<ratioindex<<", "<<numpos<<", "<<num1mm<<", "<<num2mm<<", "<<num5mm<<", "<<num2hits;
    table<<", "<<(double)numpos/(double)numpos *100.<<", "<<(double)num1mm/(double)numpos *100.<<", "<<(double)num2mm/(double)numpos *100.<<", "<<(double)num5mm/(double)numpos *100.<<", "<<(double)num2hits/(double)numpos *100.<<endl;

    delete cxyz;
    cxyz = new TCanvas("cxyz","cxyz",1050,800);
    cxyz -> Divide(2,2);
    if(ratioindex%6 == 5){ 
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
  
  for(i=0; i<6; i++){
    cdist->cd(i+1);
    gStyle -> SetOptStat(110110);
    hxy[i]->Draw("colz");
  }
  cdist->Print(Form("./%s/positiondist/%s_dist.pdf",filedir, name));

  for(i=0; i<6; i++){
    cdist->cd(i+1);
    hratio2hit[i]->Draw("");
  }
  cdist->Print(Form("./%s/positiondist/%s_dist.pdf",filedir, name));
  
  for(i=0; i<6; i++){
    cdist->cd(i+1);
    h_effdist[i]->SetMaximum(0.55);
    h_effdist[i]->SetMinimum(0);
    h_effdist[i]->Draw("");
  }
  cdist->Print(Form("./%s/positiondist/%s_dist.pdf",filedir, name));
  
  /*gStyle -> SetOptStat(0);
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
  */
  cxyz->Print(Form("./%s/positiondist/%s_xyzdist.pdf]",filedir, name));
  cdist->Print(Form("./%s/positiondist/%s_dist.pdf]",filedir, name));
}


int positiondist_basis_id(int ratioindex) {

  //ifstream poslist(Form("./inputs/geo_basis%i.txt",ratioindex)); //in principle same
  ifstream poslist(Form("./basis_wevt/seg%i_geo.txt",seg));
  ifstream decoposlist(Form("./%s/%s%i_out.csv",filedir, name,ratio[ratioindex]));
  char dummy, outname[50];

  sprintf(outname, Form("./%s/positiondist/%s_%i",filedir, name, ratio[ratioindex]));


  if(!poslist||!decoposlist){
    cerr<<"file is not found"<<endl;
    return 1;
  }
  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(111);
  //c = new TCanvas("c", "Deviation of measured positions", 900, 600);

  for(i=0; i<6; i++){
    hpos[ratioindex][0][i] = new TH1F(Form("hx%i_seg%i",i,ratioindex),Form("hx_%i to %i",5*i+5,5*i+5+5),100,-10,10);
    hpos[ratioindex][1][i] = new TH1F(Form("hy%i_seg%i",i,ratioindex),Form("hy_%i to %i",5*i+5,5*i+5+5),100,-10,10);
    hpos[ratioindex][2][i] = new TH1F(Form("hz%i_seg%i",i,ratioindex),Form("hz_%i to %i",5*i+5,5*i+5+5),100,-10,10);
    hpos[ratioindex][3][i] = new TH1F(Form("hpos%i_seg%i",i,ratioindex),Form("position dispersion: %i mm to %i mm",5*i+5,5*i+5+5),100,0,20);
    
    for(j=0; j<4; j++){
      hpos[ratioindex][j][i]->GetXaxis()->SetTitle("deviation / mm");
      hpos[ratioindex][j][i]->GetYaxis()->SetTitle("counts / 0.2 mm");
    }
  }
  

  ////////Get Pos list
  i = 0;
  getline(poslist, dummystring); //runn, x, y, z
  while(true){
    poslist >> runn_poslist[0][i] >> dummy >> x_list[0][i] >> dummy >> y_list[0][i] >> dummy >> z_list[0][i];
    poslist >> runn_poslist[1][i] >> dummy >> x_list[1][i] >> dummy >> y_list[1][i] >> dummy >> z_list[1][i];

    posdist[i] = sqrt(pow(x_list[0][i]-x_list[1][i], 2) +  pow(y_list[0][i]-y_list[1][i], 2) +  pow(z_list[0][i]-z_list[1][i], 2));

    if(i == NUMRUNS-1 || !poslist) {
      break; 
    }
    if(runn_poslist[0][i]!=runn_poslist[1][i]) cerr<<"strange geometry information"<<endl;
    i++;
  }
  numpos = i;
  //////////////////////

  getline(decoposlist, dummystring); //runnum, segnum, evtnum, avex, avey, avez, aveE
  i=0;//index for the geometry list
  k=0;//index for the decomposed parameters
  while(decoposlist){
    decoposlist >> runn[0][k] >> dummy >> dumint >> dummy >> dumint >>dummy >> segn[0][k] >> dummy >> x_meas[0][k] >> dummy >> y_meas[0][k] >> dummy >> z_meas[0][k] >> dummy >> efrac[0][k] >> dummy >> dumdouble >> dummy >> dumdouble >> dummy >> dumdouble >> dummy >> dumdouble >> dummy >> dumdouble;
    if(dumint!=1) {
      cerr<<"more than 2 events or some other error"<<endl;
      return 1;
    }
    //for second event
    while(decoposlist){
      decoposlist >> runn[1][k] >> dummy >> dumint >> dummy >> dumint >>dummy >> segn[1][k] >> dummy >> x_meas[1][k] >> dummy >> y_meas[1][k] >> dummy >> z_meas[1][k] >> dummy >> efrac[1][k] >> dummy >> dumdouble >> dummy >> dumdouble >> dummy >> dumdouble >> dummy >> dumdouble >> dummy >> dumdouble;
      if(dumint==2){
	h_effdist[ratioindex]->Fill(posdist[i]);
	h_totevts[ratioindex]->Fill(posdist[i]);
	break;
      }else if(dumint==1){// 0 was only 1 event decomposition. Discard.
	runn[0][k] = runn[1][k];
	segn[0][k] = segn[1][k];
	x_meas[0][k] = x_meas[1][k];
	y_meas[0][k] = y_meas[1][k];
	z_meas[0][k] = z_meas[1][k];
	efrac[0][k] = efrac[1][k];
	h_totevts[ratioindex]->Fill(posdist[i]);
	i++;
	//cout<<"one eve"<<endl;
      }else{
	cerr<<"err."<<endl;
	return 1;
      }
    }
    if (!decoposlist) {
      cout<<"End of decomp list of seg"<<ratioindex<<": numpos ="<<numpos<<", numdecomp ="<<k<<endl;
      break;
    }


    /*
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
    */

    //2D plot
    hxy[ratioindex]->Fill(x_meas[0][k],y_meas[0][k]);
    hxy[ratioindex]->Fill(x_meas[1][k],y_meas[1][k]);
    hratio2hit[ratioindex]->Fill(efrac[0][k]<efrac[1][k]?efrac[0][k]/(efrac[0][k]+efrac[1][k]):efrac[1][k]/(efrac[0][k]+efrac[1][k]));
    
    for(j=0; j<2; j++){
      x[j][k] = x_list[j][i];
      y[j][k] = y_list[j][i];
      z[j][k] = z_list[j][i];
      r[j][k] = sqrt(x[j][k]*x[j][k]+y[j][k]*y[j][k]);
      //cout<<"inv:"<<invertflag<<" r:"<<x[j][k]<<" r1:"<<y[j][k]<<endl;
      deltapos[j][k] = sqrt(pow(x_meas[0][k]-x[j][k], 2) +  pow(y_meas[0][k]-y[j][k], 2) +  pow(z_meas[0][k]-z[j][k], 2));
    }
    if(deltapos[0][k]<deltapos[1][k]){
      invertflag=0;
      deltapos[1][k] = sqrt(pow(x_meas[1][k]-x[1][k], 2) +  pow(y_meas[1][k]-y[1][k], 2) +  pow(z_meas[1][k]-z[1][k], 2));      
    }else{
      invertflag=1;
      //deltapos[0][k] = deltapos[1][k];
      deltapos[0][k] = sqrt(pow(x_meas[1][k]-x[0][k], 2) +  pow(y_meas[1][k]-y[0][k], 2) +  pow(z_meas[1][k]-z[0][k], 2));      
    }
    
    //cout<<"inv:"<<invertflag<<" r:"<<r[0][k]<<endl;

    for(int ii=0; ii<2; ii++){
      if(deltapos[ii][k]<1){
	num1mm++;
	num2mm++;
	num5mm++;
      }else if(deltapos[ii][k]<2.){
	num2mm++;
	num5mm++;
      }else if(deltapos[ii][k]<5.){
	num5mm++;
      }
      for(j=0;j<6;j++){
	if(j==5 || (r[ii][k]>5.*((double)j+1.)&&r[ii][k]<=5.*((double)j+2.))){
	  hpos[ratioindex][0][j]->Fill(x_meas[ii][k]);
	  hpos[ratioindex][1][j]->Fill(y_meas[ii][k]);
	  hpos[ratioindex][2][j]->Fill(z_meas[ii][k]);
	  hpos[ratioindex][3][j]->Fill(deltapos[ii][k]);
	  break;
	}
      }
    }
    i++;
    k++;
    num2hits++;
  }
  cout<<endl<<" Numpos = "<<numpos<<" num 2-hits-evt="<<k<<endl<<endl;

  i=3;
  for (j=0; j<6; j++){
    cdist->cd(j+1);
    gPad->SetLogy(0);
    hpos[ratioindex][i][j]->SetLineColor(ratioindex%6+1);
    if (ratioindex%6==0){
      hpos[ratioindex][i][j]->Draw("");
    }else{
      hpos[ratioindex][i][j]->Draw("same");
      if(hpos[0][i][j]->GetMaximum()<1.1 * hpos[ratioindex][i][j]->GetMaximum()){
	hpos[0][i][j]->SetMaximum(1.1 * hpos[ratioindex][i][j]->GetMaximum());
      }
    }
  }
  leg->AddEntry(hpos[ratioindex][i][0],Form("%i:%i",ratio[ratioindex],10-ratio[ratioindex]),"l");

  //if(draw) 
  cxyz->Print(Form("./%s/positiondist/%s_xyzdist.pdf",filedir, name));
   

  /*for (i=0; i<4; i++){
    delete g[i];
    }*/

  h_effdist[ratioindex]->Divide(h_effdist[ratioindex],h_totevts[ratioindex]);

  return 0;
}
