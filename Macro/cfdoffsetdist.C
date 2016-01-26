#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;

ifstream runlist("./runlist.txt");
char dummy, outname[] = "./cevtout/offsetdist";
string dummystring;
int runn, seg, evtn, deltat;
double timing;

int i, j, k;
TCanvas *c;
TH1F *h[17]; // num of Runs -> 17

ofstream fitresult (Form("%s/fitresult.csv",outname));

void cfdoffsetdist(void) {
  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(111);
  c = new TCanvas("c", "gaus fit", 900, 600);

  fitresult << "Run number, Segmant number, average timing(ns), timing distribution (ns) FWHM" <<endl;

  c->Print(Form("%s/fitresult.pdf[",outname));

  for(k = 0; k < 17;  k++){
    // for(k = 0; k < 1;  k++){
    runlist >> dummystring >> runn >> seg;
    if(!runlist) break;
    cout << dummystring << " run:"<<runn <<" seg:"<<seg<<endl;
    h[k] = new TH1F(Form("h%i", k), Form("run: %i, seg: %i", runn, seg), 50, 130, 180);
    h[k] -> GetXaxis() -> SetTitle("time / 1 ns");
    h[k] -> GetYaxis() -> SetTitle("Count / 1 ns");
    /*    h[k] -> SetMaximum(10); 
	  h[k] -> SetMinimum(-10);*/
    ifstream offsetdata(Form("./cevtout/offset/WFOUT-Run0%iSegment%i.dat.csv", runn, seg));
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
      cout << evtn << dummy << timing << dummy << deltat <<endl;
      h[k] -> Fill (timing);
    }
    h[k]->Draw("");
    h[k]->Fit("gaus");
    TF1 *f = func = h[k]->GetFunction("gaus");
    fitresult << runn <<", " <<seg <<", "<< f->GetParameter(1)<<", "<< f->GetParameter(2) * 2.355 <<endl;
    c->Print(Form("%s/fitresult.pdf",outname));
  }
  c->Print(Form("%s/fitresult.pdf]",outname));
}
