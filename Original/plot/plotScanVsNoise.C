{
  
  TChain* tree = new TChain("RecoAndSim");
//   tree->Add("output/advanced.multifit.mysample_11_-13.000_0.000_10_25.00_10.00_0.00_*.000_0.00_CRRC43.root");
  tree->Add("output/advanced.multifit.mysample_100_-13.000_0.000_10_25.00_10.00_0.00_*.000_0.00_CRRC43.root");
  
  TH2F* histo = new TH2F ("histo", "", 20, 0, 20*0.044, 200, -2, 2);
  
  tree->Draw("best_pedestal:sigmaNoise >> histo","","goff");
  histo->GetXaxis()->SetTitle("noise rms [GeV]");
  histo->GetYaxis()->SetTitle("Pedestal [GeV]");
  
  histo->Draw("colz");
  
  profileX = histo->ProfileX();
  
  gPad->SetGrid();
  
  TCanvas* cc = new TCanvas("cc", "", 800, 600);
  
  histo->Draw("AXIS");
  
  profileX->SetMarkerSize(1);
  profileX->SetMarkerStyle(20);
  profileX->SetMarkerColor(kBlue);
  profileX->SetLineColor(kBlue);
  profileX->SetLineWidth(2);
  
  profileX->Draw("PL same");
  profileX->GetXaxis()->SetTitle("noise rms [GeV]");
  profileX->GetYaxis()->SetTitle("Pedestal [GeV]");
  
  gPad->SetGrid();
  
  
  
  
  
  
  
  TCanvas* ccEnergy = new TCanvas("ccEnergy", "Advanced", 800, 600);
  
  
  TH2F* histoEnergy = new TH2F ("histoEnergy", "", 20, 0, 20*0.044, 100, 8, 12);
  
  tree->Draw("samplesReco[4]:sigmaNoise >> histoEnergy","","goff");
  histoEnergy->GetXaxis()->SetTitle("noise rms [GeV]");
  histoEnergy->GetYaxis()->SetTitle("Energy [GeV]");
  
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
  profileX_energy->GetXaxis()->SetTitle("noise rms [GeV]");
  profileX_energy->GetYaxis()->SetTitle("Energy [GeV]");
  
  gPad->SetGrid();
  
  
  



  TCanvas* ccEnergySimple = new TCanvas("ccEnergySimple", "Simple", 800, 600);
  
  
  TH2F* histoEnergySimple = new TH2F ("histoEnergySimple", "", 20, 0, 20*0.044, 100, 8, 12);
  
  //---- 200 = iMAX_pedestals/2  --> 0 shift
  
  tree->Draw("complete_samplesReco[200][4]:sigmaNoise >> histoEnergySimple","","goff");
  histoEnergySimple->GetXaxis()->SetTitle("noise rms [GeV]");
  histoEnergySimple->GetYaxis()->SetTitle("Energy [GeV]");
  
  histoEnergySimple->Draw("colz");
  
  profileX_energySimple = histoEnergySimple->ProfileX();
  
  profileX_energySimple->SetMarkerSize(1);
  profileX_energySimple->SetMarkerStyle(25);
  profileX_energySimple->SetMarkerColor(kBlue);
  profileX_energySimple->SetLineColor(kBlue);
  profileX_energySimple->SetLineWidth(2);
  
  profileX_energySimple->Draw("APL same");
  
  gPad->SetGrid();


  
  TCanvas* ccEnergyProfileSimple = new TCanvas("ccEnergyProfileSimple", "Simple vs Advanced", 800, 600);
  
  histoEnergySimple->Draw("AXIS");
//   histoEnergySimple->GetYaxis()->SetRangeUser(9.8, 10.4);
  
  profileX_energySimple->Draw("APL same");
  
  profileX_energy->Draw("PL same");
  
  profileX_energySimple->GetXaxis()->SetTitle("noise rms [GeV]");
  profileX_energySimple->GetYaxis()->SetTitle("Energy [GeV]");
  
  gPad->SetGrid();
  
  
  
}




