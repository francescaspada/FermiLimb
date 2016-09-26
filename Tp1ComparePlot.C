void Tp1ComparePlot() {

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","Gammas from models",200,10,700,800);
  c1->SetFillColor(0);
  c1->cd();
  pad1 = new TPad("pad1","Spectrum from Kamae",0.05,0.50,0.95,0.95,21);
  pad2 = new TPad("pad2","Spectrum from Ostap",0.05,0.05,0.95,0.45,21);
  pad1->Draw();
  pad2->Draw();

  /* Get spectrum from files */
  pad1->cd();
  pad1->SetFillColor(0);
  TGraph *RawGraphK = new TGraph("kamae/gammaspectrum-Tp1TeV.csv");
  RawGraphK->SetTitle("[Kamae] gamma-rays from 1 TeV pp");
  RawGraphK->SetMarkerStyle(21);
  RawGraphK->SetMarkerSize(.6);
  RawGraphK->SetMarkerColor(kAzure-4);
  RawGraphK->GetXaxis()->SetTitle("log_{10}(R [GV])");
  RawGraphK->GetYaxis()->SetTitle("log(E #times d#sigma/dlogE) [E in GV]");
  RawGraphK->Draw("AP");

  pad2->cd();
  pad2->SetFillColor(0);
  TGraph *RawGraphO = new TGraph("ostap/p+p-gamma-Tp1TeV.dat");
  RawGraphO->SetTitle("[Ostap] gamma-rays from 1 TeV pp");
  RawGraphO->SetMarkerStyle(20);
  RawGraphO->SetMarkerSize(.6);
  RawGraphO->SetMarkerColor(kMagenta-4);
  RawGraphO->GetXaxis()->SetTitle("log(R [GV])");
  RawGraphO->GetYaxis()->SetTitle("log(E #times d#sigma/dlogE) [E in GV]");
  RawGraphO->Draw("AP");


}
