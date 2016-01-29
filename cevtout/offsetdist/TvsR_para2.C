{
//=========Macro generated from canvas: c_gr/Timing Distributions
//=========  (Fri Jan 29 10:22:45 2016) by ROOT version5.32/00
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
   
   TH1F *Graph_Graph_Graph12 = new TH1F("Graph_Graph_Graph12","Timing Distributions",100,0,30);
   Graph_Graph_Graph12->SetMinimum(154);
   Graph_Graph_Graph12->SetMaximum(162.331);
   Graph_Graph_Graph12->SetDirectory(0);
   Graph_Graph_Graph12->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph12->SetLineColor(ci);
   Graph_Graph_Graph12->GetXaxis()->SetTitle("Radius / mm");
   Graph_Graph_Graph12->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph12->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph12->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph12->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph12->GetYaxis()->SetTitle("Trigger timing / ns");
   Graph_Graph_Graph12->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph12->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph12->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph12->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph12->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph12->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph12->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph12->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph_Graph12);
   
   
   TF1 *fpar2 = new TF1("fpar2","pol2",0,30);
   fpar2->SetFillColor(19);
   fpar2->SetFillStyle(0);
   fpar2->SetLineColor(2);
   fpar2->SetLineWidth(2);
   fpar2->SetChisquare(93.62784);
   fpar2->SetNDF(14);
   fpar2->GetXaxis()->SetLabelFont(42);
   fpar2->GetXaxis()->SetLabelSize(0.035);
   fpar2->GetXaxis()->SetTitleSize(0.035);
   fpar2->GetXaxis()->SetTitleFont(42);
   fpar2->GetYaxis()->SetLabelFont(42);
   fpar2->GetYaxis()->SetLabelSize(0.035);
   fpar2->GetYaxis()->SetTitleSize(0.035);
   fpar2->GetYaxis()->SetTitleFont(42);
   fpar2->SetParameter(0,157.3081);
   fpar2->SetParError(0,0.2757017);
   fpar2->SetParLimits(0,0,0);
   fpar2->SetParameter(1,-0.4072944);
   fpar2->SetParError(1,0.03227173);
   fpar2->SetParLimits(1,0,0);
   fpar2->SetParameter(2,0.01984656);
   fpar2->SetParError(2,0.0008770588);
   fpar2->SetParLimits(2,0,0);
   gre->GetListOfFunctions()->Add(fpar2);
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
