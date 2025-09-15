{
  
  TChain* tree = new TChain("RecoAndSim");
  //   tree->Add("output/advanced.multifit.mysample_100_-13.000_0.000_10_25.00_10.00_80.00_0.000_*.00_CRRC43.root");    //---- no noise
  //   tree->Add("output/advanced.multifit.mysample_100_-13.000_0.000_10_25.00_10.00_80.00_1.000_*.00_CRRC43.root");    //---- noise x1
  //     tree->Add("output/advanced.multifit.mysample_100_-13.000_0.000_10_25.00_10.00_80.00_10.000_*.00_CRRC43.root");    //---- noise x10
  //   tree->Add("output/advanced.multifit.mysample_100_-13.000_0.000_10_25.00_10.00_80.00_0.000_*.00_CRRC43_-1.00.root");   //---- pedestal -1
//     tree->Add("output/advanced.multifit.mysample_100_-13.000_0.000_10_25.00_10.00_80.00_0.000_*.00_CRRC43_1.00.root");    //---- pedestal +1


//     tree->Add("outputExternal/advanced.multifit.mysample_200_-13.000_0.000_10_25.00_10.00_80.00_0.000_*.00_CRRC43.root");    //---- no noise
//     tree->Add("outputExternal/advanced.multifit.mysample_200_-13.000_0.000_10_25.00_10.00_80.00_1.000_*.00_CRRC43.root");    //---- noise x1
//       tree->Add("outputExternal/advanced.multifit.mysample_200_-13.000_0.000_10_25.00_10.00_80.00_10.000_*.00_CRRC43.root");    //---- noise x10
//     tree->Add("outputExternal/advanced.multifit.mysample_200_-13.000_0.000_10_25.00_10.00_80.00_0.000_*.00_CRRC43_-1.00.root");   //---- pedestal -1
      tree->Add("outputExternal/advanced.multifit.mysample_200_-13.000_0.000_10_25.00_5.00_80.00_0.000_*.00_CRRC43_-1.00.root");   //---- pedestal -1   pulse 5 GeV
  //     tree->Add("outputExternal/advanced.multifit.mysample_200_-13.000_0.000_10_25.00_10.00_80.00_0.000_*.00_CRRC43_1.00.root");    //---- pedestal +1
  
  TH2F* histo = new TH2F ("histo", "", 39, 1, 40, 1600, -8, 8);
  
  tree->Draw("best_pedestal:puFactor >> histo","","goff");
  histo->GetXaxis()->SetTitle("pileup");
  histo->GetYaxis()->SetTitle("Pedestal [GeV]");
  
  histo->Draw("colz");
  
  profileX = histo->ProfileX();
  
  gPad->SetGrid();
  
  TCanvas* cc = new TCanvas("cc", "", 800, 600);
  
  histo->Draw("AXIS");
  
  profileX->SetMarkerSize(1);
  profileX->SetMarkerStyle(20);
  profileX->SetMarkerColor(kCyan);
  profileX->SetLineColor(kCyan);
  profileX->SetLineWidth(2);
  
  profileX->Draw("PL same");
  profileX->GetXaxis()->SetTitle("pileup");
  profileX->GetYaxis()->SetTitle("Pedestal [GeV]");
  
  gPad->SetGrid();
  
  
  
  
  
  
  
  TCanvas* ccEnergy = new TCanvas("ccEnergy", "Advanced", 800, 600);
  
  
  TH2F* histoEnergy = new TH2F ("histoEnergy", "", 39, 1, 40, 1600, -8, 8);
  
  tree->Draw("(samplesReco[5]-amplitudeTruth):puFactor >> histoEnergy","","goff");
  histoEnergy->GetXaxis()->SetTitle("pileup");
  histoEnergy->GetYaxis()->SetTitle("Reco Energy - True energy [GeV]");
  
  histoEnergy->Draw("colz");
  
  profileX_energy = histoEnergy->ProfileX();
  
  gPad->SetGrid();
  
  TCanvas* ccEnergyProfile = new TCanvas("ccEnergyProfile", "Advanced", 800, 600);
  
  histoEnergy->Draw("colz");
  
  profileX_energy->SetMarkerSize(1);
  profileX_energy->SetMarkerStyle(23);
  profileX_energy->SetMarkerColor(kRed);
  profileX_energy->SetLineColor(kRed);
  profileX_energy->SetLineWidth(2);
  
  profileX_energy->Draw("PL same");
  profileX_energy->GetXaxis()->SetTitle("pileup");
  profileX_energy->GetYaxis()->SetTitle("Reco Energy - True energy [GeV]");
  
  
  // Fit slices projected along Y fron bins in X [0,40] with more than -1 bins in Y filled
  histoEnergy->FitSlicesY(0,0,40,0, "QN");
  
  histoEnergy_1->SetMarkerSize(1);
  histoEnergy_1->SetMarkerStyle(32);
  histoEnergy_1->SetMarkerColor(97);
  histoEnergy_1->SetLineColor(97);
  histoEnergy_1->SetLineWidth(2);
  
  histoEnergy_2->SetMarkerSize(1);
  histoEnergy_2->SetMarkerStyle(32);
  histoEnergy_2->SetMarkerColor(97);
  histoEnergy_2->SetLineColor(97);
  histoEnergy_2->SetLineWidth(2);
  
  histoEnergy_1->Draw("PL same");
  
  
  
  gPad->SetGrid();
  
  
  



  TCanvas* ccEnergySimple = new TCanvas("ccEnergySimple", "Simple", 800, 600);
  
  
  TH2F* histoEnergySimple = new TH2F ("histoEnergySimple", "", 39, 1, 40, 1600, -8, 8);
  
  //---- 200 = iMAX_pedestals/2  --> 0 shift
  
  tree->Draw("(complete_samplesReco[200][5]-amplitudeTruth):puFactor >> histoEnergySimple","","goff");
  histoEnergySimple->GetXaxis()->SetTitle("pileup");
  histoEnergySimple->GetYaxis()->SetTitle("Reco Energy - True energy [GeV]");
  
  histoEnergySimple->Draw("colz");
  
  profileX_energySimple = histoEnergySimple->ProfileX();
  
  profileX_energySimple->SetMarkerSize(1);
  profileX_energySimple->SetMarkerStyle(21);
  profileX_energySimple->SetMarkerColor(kCyan+1);
  profileX_energySimple->SetLineColor(kCyan+1);
  profileX_energySimple->SetLineWidth(2);
  
  profileX_energySimple->Draw("APL same");
  
  
  
  // Fit slices projected along Y fron bins in X [0,40] with more than -1 bins in Y filled
  histoEnergySimple->FitSlicesY(0,0,40,0, "QN");
  
  histoEnergySimple_1->SetMarkerSize(1);
  histoEnergySimple_1->SetMarkerStyle(25);
  histoEnergySimple_1->SetMarkerColor(8);
  histoEnergySimple_1->SetLineColor(8);
  histoEnergySimple_1->SetLineWidth(2);

  histoEnergySimple_2->SetMarkerSize(1);
  histoEnergySimple_2->SetMarkerStyle(25);
  histoEnergySimple_2->SetMarkerColor(8);
  histoEnergySimple_2->SetLineColor(8);
  histoEnergySimple_2->SetLineWidth(2);
  
  histoEnergySimple_1->Draw("PL same");
  
  gPad->SetGrid();


  
  TCanvas* ccEnergyProfileSimple = new TCanvas("ccEnergyProfileSimple", "Simple vs Advanced", 1200, 1200);
  ccEnergyProfileSimple->Divide(2,2);
  
  ccEnergyProfileSimple->cd(1);
  histoEnergySimple->Draw("AXIS");
  //   histoEnergySimple->GetYaxis()->SetRangeUser(9.8, 10.4);
  
  profileX_energySimple->Draw("APL same");
  
  profileX_energy->Draw("PL same");
  
  profileX_energySimple->GetXaxis()->SetTitle("pileup");
  profileX_energySimple->GetYaxis()->SetTitle("Reco Energy - True energy [GeV]");
   
  gPad->SetGrid();
  
  
  ccEnergyProfileSimple->cd(2);
  histoEnergySimple->Draw("AXIS");
  //   histoEnergySimple->GetYaxis()->SetRangeUser(9.8, 10.4);
  
  histoEnergySimple_1->Draw("APL same");
  
  histoEnergy_1->Draw("PL same");
  
  histoEnergySimple_1->GetXaxis()->SetTitle("pileup");
  histoEnergySimple_1->GetYaxis()->SetTitle("Reco Energy - True energy [GeV]");
  
  
  gPad->SetGrid();
  



  ccEnergyProfileSimple->cd(4);

  histoEnergySimple_2->Draw("PL");
  
  histoEnergy_2->Draw("PL same");
  
  histoEnergySimple_2->GetXaxis()->SetTitle("pileup");
  histoEnergySimple_2->GetYaxis()->SetTitle("#sigma_{Energy} [GeV]");
  
  
  gPad->SetGrid();
  
  
  
  
  
  
  ccEnergyProfileSimple->cd(3);
  
  //---- sigma_energy/energy
  
  profileX_rms_energy = histoEnergy->ProfileX("_rms",0, -1, "s"); 
  TGraphErrors* gr_profileX_rms_energy = new TGraphErrors();
  
  for (int i=0; i<profileX_rms_energy->GetNbinsX(); i++) {
    float rms = profileX_rms_energy->GetBinError(i+1);
    float x = profileX_rms_energy->GetXaxis()->GetBinCenter(i+1);
    
    gr_profileX_rms_energy->SetPoint(i, x, rms / profileX_energy->GetBinContent(i+1));
    gr_profileX_rms_energy->SetPointError(i, profileX_energy->GetBinError(i+1)/rms*profileX_energy->GetBinError(i+1) / profileX_energy->GetBinContent(i+1));
  }
  
  gr_profileX_rms_energy->SetMarkerSize(1);
  gr_profileX_rms_energy->SetMarkerStyle(23);
  gr_profileX_rms_energy->SetMarkerColor(kRed);
  gr_profileX_rms_energy->SetLineColor(kRed);
  gr_profileX_rms_energy->SetLineWidth(2);
  
  
  profileX_rms_energySimple = histoEnergySimple->ProfileX("_rms2",0, -1, "s");  
  TGraphErrors* gr_profileX_rms_energySimple = new TGraphErrors();
  
  for (int i=0; i<profileX_rms_energySimple->GetNbinsX(); i++) {
    float rms = profileX_rms_energySimple->GetBinError(i+1);
    float x = profileX_rms_energySimple->GetXaxis()->GetBinCenter(i+1);
    
    gr_profileX_rms_energySimple->SetPoint(i, x, rms/profileX_energySimple->GetBinContent(i+1));
    gr_profileX_rms_energySimple->SetPointError(i, profileX_energySimple->GetBinError(i+1)/rms*profileX_energySimple->GetBinError(i+1)/profileX_energySimple->GetBinContent(i+1));

  }
  
  gr_profileX_rms_energySimple->SetMarkerSize(1);
  gr_profileX_rms_energySimple->SetMarkerStyle(21);
  gr_profileX_rms_energySimple->SetMarkerColor(kCyan+1);
  gr_profileX_rms_energySimple->SetLineColor(kCyan+1);
  gr_profileX_rms_energySimple->SetLineWidth(2);
  
  
  
  gr_profileX_rms_energy->Draw("APL");
  gr_profileX_rms_energy->GetXaxis()->SetTitle("pileup");
  gr_profileX_rms_energy->GetYaxis()->SetTitle("#sigma_{Energy}/Energy");
  gr_profileX_rms_energy->GetXaxis()->SetRangeUser(0,40);
  
  gr_profileX_rms_energySimple->Draw("PL");
  
  gPad->SetGrid();
  
  
  
  
  
//   histo->Draw("AXIS");
//   
//   profileX->SetMarkerSize(1);
//   profileX->SetMarkerStyle(20);
//   profileX->SetMarkerColor(kCyan);
//   profileX->SetLineColor(kCyan);
//   profileX->SetLineWidth(2);
//   
//   profileX->Draw("PL same");
//   profileX->GetXaxis()->SetTitle("pileup");
//   profileX->GetYaxis()->SetTitle("Pedestal [GeV]");
//   
//   gPad->SetGrid();
  
  
  
  
  
}




