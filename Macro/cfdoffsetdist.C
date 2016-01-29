#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;

const int NUMRUNS = 17;

ifstream runlist("./runlists/runlist.txt");
ifstream poslist("./inputs/geometry.txt");
char dummy, outname[] = "./cevtout/offsetdist";
string dummystring;
int runn[NUMRUNS], runn_poslist[NUMRUNS], segn[NUMRUNS], evtn, deltat, dumint;
double timing, par[2][NUMRUNS], parerr[2][NUMRUNS], x[NUMRUNS], y[NUMRUNS], z[NUMRUNS], r[NUMRUNS], x_list[NUMRUNS], y_list[NUMRUNS], z_list[NUMRUNS];

int i, j, k, numpos, draw = 1;
TCanvas *c, *c_gr;
TH1F *h[NUMRUNS]; // num of Runs -> NUMRUNS
TGraphErrors *gr;
TF1 *f;
TF1 *fpar = new TF1("fpar","[0]+[1]*(x-[2])*(x-[2])");
TF1 *fpar2 = new TF1("fpar2","pol2", 0, 30);
TF1 *fexp = new TF1("fexp","pol0(0)+expo(1)", 0, 30);

ofstream fitresult (Form("%s/fitresult.csv",outname));
ofstream TvsRresult (Form("%s/TvsRfitresult.csv",outname));

void cfdoffsetdist(void) {

  if(!runlist){
    cerr<<"runlist is not found"<<endl;
    break;
  }
  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(111);
  c = new TCanvas("c", "gaus fit", 900, 600);
    

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

  fitresult << "Run number, Segmant number, average timing(ns), timing distribution (ns) FWHM" <<endl;

  if(draw) c->Print(Form("%s/fitresult.pdf[",outname));

  for(k = 0; k < NUMRUNS;  k++){
    // for(k = 0; k < 1;  k++){
    runlist >> dummystring >> runn[k] >> segn[k];
    cout <<endl<< dummystring << " run:"<<runn[k] <<" seg:"<<segn[k]<<endl;
    h[k] = new TH1F(Form("h%i", k), Form("run: %i, seg: %i", runn[k], segn[k]), 100, 130, 180);
    h[k] -> GetXaxis() -> SetTitle("Trigger timing / 0.5 ns");
    h[k] -> GetYaxis() -> SetTitle("Count / ns");
    /*    h[k] -> SetMaximum(10); 
	  h[k] -> SetMinimum(-10);*/
    ifstream offsetdata(Form("./cevtout/offset/WFOUT-Run0%iSegment%i.dat.csv", runn[k], segn[k]));
    if(!offsetdata) {
      cerr << "File was not found." <<endl;
      break;
    }
    getline(offsetdata, dummystring);
    cout << dummystring <<endl;
    getline(offsetdata, dummystring);
    cout << dummystring <<endl ;
    getline(offsetdata, dummystring);
    cout << dummystring <<endl;
    //while(offsetdata)
    while(true){
      offsetdata >> evtn >> dummy >> timing >> dummy >> deltat;
      if(!offsetdata) break;
      //cout << evtn << dummy << timing << dummy << deltat <<endl;
      h[k] -> Fill (timing);
    }
    h[k]->Draw("");
    h[k]->Fit("gaus", "L");
    f = func = h[k]->GetFunction("gaus");
    par[0][k] = f->GetParameter(1);
    par[1][k] = f->GetParameter(2) * 2.355;
    parerr[0][k] = f->GetParError(1);
    parerr[1][k] = f->GetParError(2) * 2.355;
    fitresult << runn[k] <<", " <<segn[k] <<", "<<f->GetParameter(1)<<", "<< f->GetParameter(2) * 2.355 <<endl;
    if(draw) c->Print(Form("%s/fitresult.pdf",outname));
 

    /// 2nd step ///
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
	//cout<<r[k]<<endl;
	break;
      }
    }
    /// ///
  }
  if(draw) c->Print(Form("%s/fitresult.pdf]",outname));

  ///////
  c_gr = new TCanvas("c_gr", "Timing Distributions", 900, 600);
  c_gr->cd();
  gStyle -> SetOptFit(0);
  gr= new TGraphErrors(NUMRUNS, r, &par[0][0], 0 , &parerr[0][0]);
  //gr= new TGraphErrors(NUMRUNS, r, &par[0][0], 0 , &par[1][0] );  
  gr -> SetTitle("Timing Distributions");
  gr -> GetXaxis() -> SetTitle("Radius / mm");
  gr -> GetYaxis() -> SetTitle("Trigger timing / ns");
  gr -> GetXaxis() -> SetLimits(0,30.);
  //gr -> SetMarkerStyle(21);
  gr -> Draw("AP*");

  fpar -> SetParameter(2, 5.);
  fpar -> FixParameter(2, 5.);
  gr -> Fit("fpar", "");
  gr -> Fit("fpar", "M");
  gr->SetMinimum((int)(fpar->GetParameter(0)));

  TvsRresult << "name, function, par0, par1, par2, parerr1, 2, 3"<<endl;
  TvsRresult << "fpar, [0]+[1]*(x-[2])*(x-[2]), "<< fpar->GetParameter(0) <<", "<< fpar->GetParameter(1) <<", "<< fpar->GetParameter(2)
	     <<", "<< fpar->GetParError(0)  <<", "<< fpar->GetParError(1)  <<", "<< fpar->GetParError(2) <<endl;
  if(draw) c_gr->Print(Form("%s/TvsR_parabolic_5mmCC.pdf", outname));
  
  gr -> Fit("fpar2", "");
  gr -> Fit("fpar2", "M");
  TvsRresult << "fpar2, pol2, "<< fpar2->GetParameter(0) <<", "<< fpar2->GetParameter(1) <<", "<< fpar2->GetParameter(2)
	     <<", "<< fpar2->GetParError(0)  <<", "<< fpar2->GetParError(1)  <<", "<< fpar2->GetParError(2) <<endl;
  if(draw) c_gr->Print(Form("%s/TvsR_para2.pdf", outname));
   
  gr -> Fit("fexp", "");
  gr -> Fit("fexp", "M");
  TvsRresult << "fexp, pol0(0)+expo(1), "<< fexp->GetParameter(0) <<", "<< fexp->GetParameter(1) <<", "<< fexp->GetParameter(2)
	     <<", "<< fexp->GetParError(0)  <<", "<< fexp->GetParError(1)  <<", "<< fexp->GetParError(2) <<endl;
  if(draw) c_gr->Print(Form("%s/TvsR_exp.pdf", outname));
  
  //Create several fittings compere ?

  //gr -> Fit("pol0+gaus[1]");
  
}
