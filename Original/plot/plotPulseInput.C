//---- plot output of multifit

void plotPulseInput (std::string nameInputFile = "output.root", int nEvent = 1){
  
  TFile *file = new TFile(nameInputFile.c_str());
  
  TTree* tree = (TTree*) file->Get("Samples");
  
  int    nWF;
  std::vector<double>* pulse_signal = new std::vector<double>;
  std::vector<double>* samples = new std::vector<double>;
  
  tree->SetBranchAddress("nWF",      &nWF);
  tree->SetBranchAddress("pulse_signal", &pulse_signal);
  tree->SetBranchAddress("samples", &samples);
  
  tree->GetEntry(nEvent);
  
//   std::cout << " nWF = " << nWF << std::endl;
  
  TCanvas* ccpulse_signal = new TCanvas ("ccpulse_signal","",800,600);
  TGraph *gr = new TGraph();
  for(int i=0; i<nWF; i++){
    gr->SetPoint(i, i, pulse_signal->at(i));
//     std::cout << " i = " << pulse_signal->at(i) << std::endl;
  }
  gr->Draw("AL");
  gr->SetLineColor(kMagenta);
  gr->SetLineWidth(2);
  gr->GetXaxis()->SetTitle("time [ns]");
  
  
  TCanvas* ccPulse = new TCanvas ("ccPulse","",800,600);
  TGraph *grPulse = new TGraph();
  for(int i=0; i<samples->size(); i++){
    grPulse->SetPoint(i, i, samples->at(i));
  }
  grPulse->SetMarkerSize(2);
  grPulse->SetMarkerStyle(21);
  grPulse->SetMarkerColor(kRed);
  grPulse->Draw("ALP");
  grPulse->GetXaxis()->SetTitle("BX");
  
  
  
  
}

