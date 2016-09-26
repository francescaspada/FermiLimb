void OKComparePlot() {

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
  TGraph *RawGraphK = new TGraph("kamae/gammaspectrum.csv");
  RawGraphK->SetTitle("Gamma-rays from BPL protons");
  RawGraphK->SetFillColor(0);
  RawGraphK->SetLineColor(0);
  RawGraphK->SetMarkerStyle(25);
  RawGraphK->SetMarkerSize(.5);
  RawGraphK->SetMarkerColor(kCyan+2);
  RawGraphK->GetXaxis()->SetTitle("log_{10}(R [GV])");
  RawGraphK->GetYaxis()->SetTitle("log(E #times d#sigma/dlogE) [E in GV]");
  RawGraphK->Draw("AP");

  pad2->cd();
  pad2->SetFillColor(0);
  TGraph *RawGraphO = new TGraph("ostap/p+p-gamma.dat");
  RawGraphO->SetTitle("Gamma-rays from BPL protons");
  RawGraphO->SetFillColor(0);
  RawGraphO->SetLineColor(0);
  RawGraphO->SetMarkerStyle(24);
  RawGraphO->SetMarkerSize(.5);
  RawGraphO->SetMarkerColor(kAzure-5);
  RawGraphO->GetXaxis()->SetTitle("log(R [GV])");
  RawGraphO->GetYaxis()->SetTitle("log(E #times d#sigma/dlogE) [E in GV]");
  RawGraphO->Draw("AP");

  TCanvas *c2 = new TCanvas("c2","Gammas from K and O models",300,10,700,400);
  RawGraphK->Draw("AP");
  RawGraphO->Draw("P");
  TLegend *legR = new TLegend(0.2,0.31,0.5,0.49);
  legR->SetTextSize(0.035);
  legR->SetLineColor(0);
  legR->SetFillColor(0);
  legR->AddEntry(RawGraphO,"AMS-02 BLP p + Kachelriess/Ostapchenko (2011)");
  legR->AddEntry(RawGraphK,"AMS-02 BLP p + Kamae (2006)");
  legR->Draw();


  TCanvas *c3 = new TCanvas("c3","Gammas from K and O models [10GeV - 10TeV]",400,10,700,400);
  c3->SetFillColor(0);
  c3->SetLogx();
  Int_t nK =  RawGraphK->GetN();
  Int_t nO =  RawGraphO->GetN();
  TVectorD logEK(nK), EdSigmadlogEK(nK), logEO(nO), EdSigmadlogEO(nO);
  TVectorD EneK(nK), FluxK(nK), EneO(nO), FluxO(nO);
  Double_t ind = 2.75;

  Int_t j=0;
  for(Int_t i=0; i<nK; i++) {
    RawGraphK->GetPoint(i,logEK[i],EdSigmadlogEK[i]);
    //printf("%3dth line: logE:%6.2f EdSigmadlogE:%7.3f ", i+1, logE[i], EdSigmadlogE[i]);
    if (10.**(logEK[i]) > 1. && 10.**(logEK[i]) < 10000.) {
      EneK[j] = 10.**(logEK[i]);
      FluxK[j] = 10.**(EdSigmadlogEK[i]) / EneK[j]**2. * EneK[j]**ind;
      //printf(" [Ene:%11.4f Flux:%11.4f]", Ene[j], Flux[j]);
      j++;
    }
    printf("\n");
  }

  j=0;
  for(Int_t i=0; i<nO; i++) {
    RawGraphO->GetPoint(i,logEO[i],EdSigmadlogEO[i]);
    //printf("%3dth line: logE:%6.2f EdSigmadlogE:%7.3f ", i+1, logE[i], EdSigmadlogE[i]);
    if (10.**(logEO[i]) > 1. && 10.**(logEO[i]) < 10000.) {
      EneO[j] = 10.**(logEO[i]);
      FluxO[j] = 10.**(EdSigmadlogEO[i]) / EneO[j]**2. * EneO[j]**ind;
      //printf(" [Ene:%11.4f Flux:%11.4f]", Ene[j], Flux[j]);
      j++;
    }
    printf("\n");
  }

  TGraph *SpecGraphK = new TGraph(EneK, FluxK);
  TGraph *SpecGraphO = new TGraph(EneO, FluxO);
  SpecGraphO->SetTitle("Gamma-rays from BPL protons - 10GV->10TV");
  SpecGraphO->SetFillColor(0);
  SpecGraphO->SetLineColor(0);
  SpecGraphO->SetMarkerStyle(24);
  SpecGraphO->SetMarkerSize(.5);
  SpecGraphO->SetMarkerColor(kAzure-5);
  SpecGraphO->GetXaxis()->SetTitle("R [GV]");
  SpecGraphO->GetYaxis()->SetTitle("Flux #times E^{2.75}");
  SpecGraphO->Draw("AP");

  SpecGraphK->SetTitle("Gamma-rays from BPL protons - 10GV->10TV");
  SpecGraphK->SetFillColor(0);
  SpecGraphK->SetLineColor(0);
  SpecGraphK->SetMarkerStyle(25);
  SpecGraphK->SetMarkerSize(.5);
  SpecGraphK->SetMarkerColor(kCyan+2);
  SpecGraphK->GetXaxis()->SetTitle("R [GV]");
  SpecGraphK->GetYaxis()->SetTitle("Flux #times E^{2.75}");
  SpecGraphK->Draw("P");

  TLegend *leg = new TLegend(0.2,0.71,0.5,0.89);
  leg->SetTextSize(0.035);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->AddEntry(SpecGraphO,"AMS-02 BLP p + Kachelriess/Ostapchenko (2011)");
  leg->AddEntry(SpecGraphK,"AMS-02 BLP p + Kamae (2006)");
  leg->Draw();


}
