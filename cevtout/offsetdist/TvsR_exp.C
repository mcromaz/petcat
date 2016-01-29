{
//=========Macro generated from canvas: c_gr/Timing Distributions
//=========  (Fri Jan 29 10:22:46 2016) by ROOT version5.32/00
   TCanvas *c_gr = new TCanvas("c_gr", "Timing Distributions",1602,54,900,600);
   gStyle->SetOptStat(0);
   c_gr->Range(-3.75,152.9586,33.75,163.3724);
   c_gr->SetFillColor(0);
   c_gr->SetBorderMode(0);
   c_gr->SetBorderSize(2);
   c_gr->SetFrameBorderMode(0);
   c_gr->SetFrameBorderMode(0);
   
   TGraphErrors *gre = new TGraphErrors(17);
   gre->SetName("Graph");
   gre->SetTitle("Timing Distributions");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(3);
   gre->SetPoint(0,27.59801,161.0221);
   gre->SetPointError(0,0,0.1381229);
   gre->SetPoint(1,9.60833,155.0833);
   gre->SetPointError(1,0,0.2066489);
   gre->SetPoint(2,26.73275,161.4875);
   gre->SetPointError(2,0,0.1658159);
   gre->SetPoint(3,26.73275,160.0948);
   gre->SetPointError(3,0,0.1358656);
   gre->SetPoint(4,17.95355,156.3827);
   gre->SetPointError(4,0,0.1134219);
   gre->SetPoint(5,17.42699,156.5);
   gre->SetPointError(5,0,0.1477271);
   gre->SetPoint(6,17.42699,156.1919);
   gre->SetPointError(6,0,0.1454249);
   gre->SetPoint(7,17.95689,156.4423);
   gre->SetPointError(7,0,0.1097124);
   gre->SetPoint(8,26.16448,160.6974);
   gre->SetPointError(8,0,0.1312992);
   gre->SetPoint(9,26.16448,159.3125);
   gre->SetPointError(9,0,0.2081702);
   gre->SetPoint(10,17.39799,156.3017);
   gre->SetPointError(10,0,0.1623658);
   gre->SetPoint(11,17.39799,155.7115);
   gre->SetPointError(11,0,0.1696566);
   gre->SetPoint(12,17.92261,156.2955);
   gre->SetPointError(12,0,0.1458999);
   gre->SetPoint(13,18.22224,156.5705);
   gre->SetPointError(13,0,0.1489429);
   gre->SetPoint(14,8.709191,155.05);
   gre->SetPointError(14,0,0.1619362);
   gre->SetPoint(15,8.709191,155.4583);
   gre->SetPointError(15,0,0.1499902);
   gre->SetPoint(16,9.113726,155.2908);
   gre->SetPointError(16,0,0.112556);
   
   TH1F *Graph_Graph_Graph_Graph123 = new TH1F("Graph_Graph_Graph_Graph123","Timing Distributions",100,0,30);
   Graph_Graph_Graph_Graph123->SetMinimum(154);
   Graph_Graph_Graph_Graph123->SetMaximum(162.331);
   Graph_Graph_Graph_Graph123->SetDirectory(0);
   Graph_Graph_Graph_Graph123->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph_Graph123->SetLineColor(ci);
   Graph_Graph_Graph_Graph123->GetXaxis()->SetTitle("Radius / mm");
   Graph_Graph_Graph_Graph123->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph_Graph123->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph_Graph123->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph_Graph123->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph_Graph123->GetYaxis()->SetTitle("Trigger timing / ns");
   Graph_Graph_Graph_Graph123->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph_Graph123->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph_Graph123->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph_Graph123->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph_Graph123->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph_Graph123->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph_Graph123->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph_Graph123->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph_Graph_Graph123);
   
   
   TF1 *fexp = new TF1("fexp","pol0(0)+expo(1)",0,30);
   fexp->SetFillColor(19);
   fexp->SetFillStyle(0);
   fexp->SetLineColor(2);
   fexp->SetLineWidth(2);
   fexp->SetChisquare(102.4473);
   fexp->SetNDF(14);
   fexp->GetXaxis()->SetLabelFont(42);
   fexp->GetXaxis()->SetLabelSize(0.035);
   fexp->GetXaxis()->SetTitleSize(0.035);
   fexp->GetXaxis()->SetTitleFont(42);
   fexp->GetYaxis()->SetLabelFont(42);
   fexp->GetYaxis()->SetLabelSize(0.035);
   fexp->GetYaxis()->SetTitleSize(0.035);
   fexp->GetYaxis()->SetTitleFont(42);
   fexp->SetParameter(0,154.7549);
   fexp->SetParError(0,0.1553176);
   fexp->SetParLimits(0,0,0);
   fexp->SetParameter(1,-2.067463);
   fexp->SetParError(1,0.2800499);
   fexp->SetParLimits(1,0,0);
   fexp->SetParameter(2,0.1431893);
   fexp->SetParError(2,0.009575046);
   fexp->SetParLimits(2,0,0);
   gre->GetListOfFunctions()->Add(fexp);
   gre->Draw("app");
   
   TPaveText *pt = new TPaveText(0.333192,0.9364634,0.666808,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *text = pt->AddText("Timing Distributions");
   pt->Draw();
   c_gr->Modified();
   c_gr->cd();
   c_gr->SetSelected(c_gr);
}
