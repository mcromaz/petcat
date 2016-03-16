#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;

//ifstream trcsv("./tr_cfd.csv");
//char outname[] = "./output/withCFD";

ifstream trcsv("./tr.csv");
//char outname[] = "./output/withoutCFD";
char outname[] = "./basis_noise/tr15";

char dummy, evt, seg;
int ch, val, dumint = 0;
int i, j, k;
char evtlabel[] = {'a', 'b', 's'}, seglabel[] = {'n', 'l', 'r', 'u', 'd', 'c'};
string dummystring;
TCanvas *c[3];
//TGraph *g[16];
TH1I *h[18];

void fillhists(int i, int dumint, TH1I *h_temp[]){
  if(dumint < 300){
    j = 0;
    k = 6 * i + j;
    h_temp[k]->SetBinContent(dumint, val);
  } else if(dumint < 600) {
    j = 1;
    k = 6 * i + j;
    h_temp[k]->SetBinContent(dumint - 300 * j, val);
  } else if(dumint < 900) {
    j = 2;
    k = 6 * i + j;
    h_temp[k]->SetBinContent(dumint - 300 * j, val);
  } else if(dumint < 1200) {
    j = 3;
    k = 6 * i + j;
    h_temp[k]->SetBinContent(dumint - 300 * j, val);
  } else if(dumint < 1500) {
    j = 4;
    k = 6 * i + j;
    h_temp[k]->SetBinContent(dumint - 300 * j, val);
  } else if(dumint < 1800) {
    j = 5;
    k = 6 * i + j;
    h_temp[k]->SetBinContent(dumint - 300 * j, val);
  }
}

void tr2hist(void) {
  gStyle -> SetOptStat(0);
  for (i = 0; i<3; i++){
    c[i] = new TCanvas(Form("c%i",i), Form("Canvas of evt %c", evtlabel[i]), 900, 600);
    c[i] -> Divide(3,2);
  }
 
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 6; j++) {
      k = 6 * i + j;
      h[k] = new TH1I(Form("h%i", k), Form("evt: %c, seg: %c", evtlabel[i], seglabel[j]), 300, 0, 300);
      h[k] -> GetXaxis() -> SetTitle("time / 10 ns");
      h[k] -> GetYaxis() -> SetTitle("ADC val / a.u.");
      if (j == 0 || j == 5) {
	h[k] -> SetMaximum(5500);
	h[k] -> SetMinimum(-100);
      } else {
	h[k] -> SetMaximum(1000);
	h[k] -> SetMinimum(-200);
      }
    }
  }

  getline(trcsv, dummystring);
  //cout<<dummystring<<endl;
  while (trcsv >> evt >> dummy >> seg >> dummy >> ch >> dummy >> val){
    dumint = dumint % 1800;
    //cout << evt << dummy << seg << dummy << ch << dummy << val << endl;
    switch (evt) {
    case 'a':
      fillhists(0, dumint, h);
      break;
    case 'b':
      fillhists(1, dumint, h);
      break;
    case 's':
      fillhists(2, dumint, h);
     break;
    }
    dumint++;
  }

  for (i = 2; i > -1 ; i--) {
    for (j = 0; j < 6; j++) {
      k = 6 * i + j;
      c[i] -> cd(j + 1);
      h[k]->Draw();
      if (i == 0) {
	h[k]->SetLineColor(2);
	c[2] -> cd(j + 1);
	h[k]->Draw("same");
      }
      if (i == 1) {
	h[k]->SetLineColor(6);
	c[2] -> cd(j + 1);
	h[k]->Draw("same");
      }
    }
  }

  c[0] -> Print(Form("%s.pdf[", outname));
  for (i = 0; i < 3; i++){
    c[i] -> Print(Form("%s.pdf", outname));
  }
  c[2] -> Print(Form("%s.pdf]", outname));
  
}
