#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;


const int NUMRUNS=20000;
char name[] ="basis_PF0";
char *index[4] = {"x", "y", "z", "pos"};
string dummystring;
int runn[NUMRUNS], runn_poslist[NUMRUNS], segn[NUMRUNS],  dumint;
double par[2][NUMRUNS], parerr[2][NUMRUNS], x[NUMRUNS], y[NUMRUNS], z[NUMRUNS], r[NUMRUNS], x_list[NUMRUNS], y_list[NUMRUNS], z_list[NUMRUNS], x_meas[NUMRUNS], y_meas[NUMRUNS], z_meas[NUMRUNS], deltapos[NUMRUNS], dumdouble;

TH1F *hx[5], *hy[5], *hz[5], *hpos[36][4][6];

int i, j, k, numpos, draw = 0;
TCanvas *c, *cxyz, *cdist;//, *c_gr;
TGraph *g[4]; //Delta x, y, z, pos
TLegend *leg;
void positiondist_basis(void) {
  cxyz = new TCanvas("cxyz","cxyz",1050,800);
  cxyz->Print(Form("./vegcatout/positiondist/%s_xyzdist.pdf[",name));
  cxyz -> Divide(2,2);
  cdist = new TCanvas("cdist","cdist",1050,800);
  cdist->Print(Form("./vegcatout/positiondist/%s_dist.pdf[",name));
  cdist ->Divide(3,2);
  leg = new TLegend(0.82,0.7,0.92,0.9);
  
  for(int segid=0; segid<36; segid++) {
    //for(int segid=0; segid<6; segid++) {
    positiondist_basis_id(segid);


    delete cxyz;
    cxyz = new TCanvas("cxyz","cxyz",1050,800);
    cxyz -> Divide(2,2);
    if(segid%6 == 5){ 
      cdist->cd(0);
      leg->SetFillColor(0);
      leg->Draw("");
      cdist->Print(Form("./vegcatout/positiondist/%s_dist.pdf",name));
      delete cdist;
      cdist = new TCanvas("cdist","cdist",1050,800);
      cdist ->Divide(3,2);
      delete leg;
      leg = new TLegend(0.82,0.7,0.92,0.9);
    }
  }
  cxyz->Print(Form("./vegcatout/positiondist/%s_xyzdist.pdf]",name));
  cdist->Print(Form("./vegcatout/positiondist/%s_dist.pdf]",name));
}

void positiondist_basis_id(int segid) {

  ifstream poslist(Form("./inputs/geo_basis%i.txt",segid));
  ifstream decoposlist(Form("./vegcatout/%s%i_out.csv",name,segid));
  
  char dummy, outname[] = Form("./vegcatout/positiondist/%s%i",name,segid);


  if(!poslist||!decoposlist){
    cerr<<"file is not found"<<endl;
    break;
  }
  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(111);
  c = new TCanvas("c", "Deviation of measured positions", 900, 600);

  for(i=0; i<6; i++){
    hpos[segid][0][i] = new TH1F(Form("hx%i_seg%i",i,segid),Form("hx_%i to %i",5*i+5,5*i+5+5),100,-10,10);
    hpos[segid][1][i] = new TH1F(Form("hy%i_seg%i",i,segid),Form("hy_%i to %i",5*i+5,5*i+5+5),100,-10,10);
    hpos[segid][2][i] = new TH1F(Form("hz%i_seg%i",i,segid),Form("hz_%i to %i",5*i+5,5*i+5+5),100,-10,10);
    hpos[segid][3][i] = new TH1F(Form("hpos%i_seg%i",i,segid),Form("hpos_%i to %i",5*i+5,5*i+5+5),100,0,20);
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

  getline(decoposlist, dummystring); //runnum, segnum, evtnum, avex, avey, avez, aveE
  for(k = 0; k < NUMRUNS;  k++){

    decoposlist >> runn[k] >> dummy >> dumint >> dummy >> dumint >>dummy >> segn[k] >> dummy >> x_meas[k] >> dummy >> y_meas[k] >> dummy >> z_meas[k] >> dummy >> dumdouble >> dummy >> dumdouble >> dummy >> dumdouble >> dummy >> dumdouble >> dummy >> dumdouble >> dummy >> dumdouble;
    if(dumint!=1) {
      k--;
      continue;
    }
    if (!decoposlist) break;
    i=k;
    x[k] = x_list[i];
    y[k] = y_list[i];
    z[k] = z_list[i];
    r[k] = sqrt(x[k]*x[k]+y[k]*y[k]);
    x_meas[k] -= x[k]; //Delta x
    y_meas[k] -= y[k];
    z_meas[k] -= z[k];
    deltapos[k] = sqrt(pow(x_meas[k], 2) +  pow(y_meas[k], 2) +  pow(z_meas[k], 2));
    //cout<<k<<" "<<r[k]<<" "<<deltapos[k]<<endl;
	
    for(i=0;i<6;i++){
      if(i==5 || r[k]>5.*((double)i+1.)&&r[k]<=5.*((double)i+2.)){
	hpos[segid][0][i]->Fill(x_meas[k]);
	hpos[segid][1][i]->Fill(y_meas[k]);
	hpos[segid][2][i]->Fill(z_meas[k]);
	hpos[segid][3][i]->Fill(deltapos[k]);
	break;
      }
      /*else if(i==5){
	cout<<"Error: this point is out of r=35"<<endl;
	}*/
    }
  }

  cout<<endl<<" Numpos = "<<numpos<<" "<<k<<endl;

  g[0]= new TGraph(numpos, r, x_meas);
  g[1]= new TGraph(numpos, r, y_meas);
  g[2]= new TGraph(numpos, r, z_meas);
  g[3]= new TGraph(numpos, r, deltapos);

  if(draw) c->Print(Form("%s.pdf[",outname));
  for (i=0; i<4; i++){
    g[i] -> SetTitle(Form("#Delta%s vs Radius",index[i]));
    g[i] -> GetXaxis() -> SetTitle("Radius / mm");
    g[i] -> GetYaxis() -> SetTitle(Form("#Delta%s / mm",index[i]));
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
    hpos[segid][i][j]->SetLineColor(segid%6+1);
    if (segid%6==0){
      hpos[segid][i][j]->Draw("");
    }else hpos[segid][i][j]->Draw("same");
  }
  leg->AddEntry(hpos[segid][i][0],Form("seg%i",segid),"l");
    
  //if(draw) 
  cxyz->Print(Form("./vegcatout/positiondist/%s_xyzdist.pdf",name));
   

  for (i=0; i<4; i++){
    delete g[i];
  }
}


