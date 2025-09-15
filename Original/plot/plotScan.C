//---- plot output of multifit

void plotScan (std::string nameInputFile = "output.root", int nEvent = 10){
  
  
  float pedestal = 0;
  
  Color_t* color = new Color_t [200];
  color[0] = kAzure; //kRed ;
  color[1] = kAzure + 10 ;
  color[2] = kYellow + 2 ;
  color[3] = kGreen ;
  color[4] = kGreen + 4 ;
  color[5] = kBlue ;
  color[6] = kCyan ;
  color[7] = kPink + 1 ;
  color[8] = kBlack ;
  color[9] = kYellow + 4 ;
  for (int i=0; i<30; i++) {
    color[i+10] = kBlue + i;
  }
  
  
  TFile *file = new TFile(nameInputFile.c_str());
  
  TTree* tree = (TTree*) file->Get("RecoAndSim");
  
  int    nWF;
  std::vector<double>* pulse_signal    = new std::vector<double>;
  std::vector<double>* pileup_signal   = new std::vector<double>;
  std::vector<double>* samplesReco = new std::vector<double>;
  std::vector<double>* samples     = new std::vector<double>;
  std::vector<double>* samples_noise = new std::vector<double>;
  std::vector<int>*    activeBXs   = new std::vector<int>;
  std::vector<double>* pulseShapeTemplate     = new std::vector<double>;
  
  float NFREQ;
  
  float best_pedestal = 0;
  float return_chi2;
  std::vector <double>* complete_chi2   = new std::vector<double>;
  std::vector <double>* complete_pedestal   = new std::vector<double>;
  std::vector < std::vector<double> >* complete_samplesReco   = new std::vector< std::vector <double> >;
  
  tree->SetBranchAddress("complete_chi2",          &complete_chi2);
  tree->SetBranchAddress("complete_pedestal",      &complete_pedestal);
  tree->SetBranchAddress("complete_samplesReco",   &complete_samplesReco);
//   tree->SetBranchAddress("best_pedestal",   &best_pedestal);
  tree->SetBranchAddress("chi2",   &return_chi2);
  
  tree->SetBranchAddress("nWF",      &nWF);
  tree->SetBranchAddress("pulse_signal", &pulse_signal);
  tree->SetBranchAddress("pileup_signal", &pileup_signal);
  tree->SetBranchAddress("samplesReco", &samplesReco);
  tree->SetBranchAddress("samples",   &samples);
  tree->SetBranchAddress("samples_noise",   &samples_noise);
  tree->SetBranchAddress("activeBXs", &activeBXs);
  tree->SetBranchAddress("nFreq",   &NFREQ);
  tree->SetBranchAddress("pulseShapeTemplate",   &pulseShapeTemplate);
  
  
  tree->GetEntry(nEvent);
  std::cout << " NFREQ = " << NFREQ << std::endl;
  
  TCanvas* ccpulse_signal = new TCanvas ("ccpulse_signal","",800,600);
  TGraph *gr = new TGraph();
  for(int i=0; i<nWF; i++){
    gr->SetPoint(i, i, pulse_signal->at(i));
  }
  gr->Draw("AL");
  gr->SetLineColor(kMagenta);
  gr->SetLineWidth(2);
  gr->GetXaxis()->SetTitle("time [ns]");
  
  TGraph *grPUall = new TGraph();
  for(int i=0; i<nWF; i++){
    grPUall->SetPoint(i, i, pileup_signal->at(i));
  }
  grPUall->Draw("L");
  grPUall->SetLineColor(kMagenta+3);
  grPUall->SetLineWidth(2);
  grPUall->SetLineStyle(3);
  grPUall->GetXaxis()->SetTitle("time [ns]");
  
  
  int IDSTART = 6*25;
  int WFLENGTH = 500*4; // step 1/4 ns in waveform
  int NSAMPLES = samplesReco->size();
  float shift_for_pu_plot = IDSTART + NSAMPLES * NFREQ;
  shift_for_pu_plot =  shift_for_pu_plot + 4*25*2 ;
  
  std::cout << " shift_for_pu_plot = " << shift_for_pu_plot << std::endl;
  std::cout << " NSAMPLES = " << NSAMPLES << std::endl;
  
  
  TGraph *grPU = new TGraph();
  for(int i=0; i<nWF; i++){
    grPU->SetPoint(i, (i-shift_for_pu_plot)/4, pileup_signal->at(i));
  }
  grPU->SetLineColor(kMagenta+3);
  grPU->SetLineWidth(2);
  grPU->SetLineStyle(3);
  grPU->GetXaxis()->SetTitle("time [ns]");
  
  ccpulse_signal->SetGrid();
  
  
  
  
  
  
  
  
  
  
  
  TCanvas* ccReco = new TCanvas ("ccReco","Reco Advanced",800,600);
//   ccReco->SetRightMargin(0.3);
  TGraph *grReco = new TGraph();
  for(int i=0; i<samplesReco->size(); i++){
    std::cout << " i, activeBXs->at(i), samplesReco->at(i) = " << i << "::" << samplesReco->size() << " -> " << activeBXs->at(i) << " , " << samplesReco->at(i) << std::endl;
    grReco->SetPoint(i, activeBXs->at(i), samplesReco->at(i));
  }
  grReco->SetMarkerSize(1);
  grReco->SetMarkerStyle(22);
  grReco->SetLineColor(kBlue);
  grReco->SetMarkerColor(kBlue);
  grReco->Draw("ALP");
  grReco->GetXaxis()->SetTitle("BX");
  
  
  TCanvas* ccPulse = new TCanvas ("ccPulse","",800,600);
  
  TGraph *grPulse_noise = new TGraph();
  for(int i=0; i<samples->size(); i++){
    grPulse_noise->SetPoint(i, i * NFREQ , samples_noise->at(i));
  }
  grPulse_noise->SetMarkerSize(1.5);
  grPulse_noise->SetMarkerStyle(21);
  grPulse_noise->SetMarkerColor(kGray);
  grPulse_noise->SetLineColor(kGray);
  grPulse_noise->SetLineStyle(3);
  
  
  
  TGraph *grPulse = new TGraph();
  for(int i=0; i<samples->size(); i++){
    grPulse->SetPoint(i, i * NFREQ , samples->at(i));
  }
  grPulse->SetMarkerSize(2);
  grPulse->SetMarkerStyle(21);
  grPulse->SetMarkerColor(kRed);
  grPulse->Draw("ALP");
  grPulse->GetYaxis()->SetRangeUser(0 > (best_pedestal-1) ? (best_pedestal-1)  : 0 , samplesReco->at(5)+10+best_pedestal);
  grPulse->GetXaxis()->SetTitle("time [ns]");
  
  grPU->Draw("L");
  
  ccPulse->SetGrid();
  
  std::cout << " end " << std::endl;
  
  
  TCanvas* ccPulseAndReco = new TCanvas ("ccPulseAndReco","Advanced",800,1200);
  ccPulseAndReco->Divide(1,2);
  ccPulseAndReco->cd(1);
  TGraph *grPulseRecoAll = new TGraph();
  TGraph *grPulseReco_OOT_All = new TGraph();
  TGraph *grPulseReco[samplesReco->size()];
  std::cout << " samplesReco->size() = " << samplesReco->size() << std::endl;
  std::cout << " activeBXs->size() = " << activeBXs->size() << std::endl;
  std::cout << " pulseShapeTemplate->size() = " << pulseShapeTemplate->size() << std::endl;
  std::cout << " samples->size() = " << samples->size() << std::endl;
  
  TLegend* leg = new TLegend(0.91,0.10,0.99,0.90);
  
  std::vector<float> totalRecoSpectrum;
  for(int i=0; i<samples->size(); i++){
    totalRecoSpectrum.push_back(0);
  }
  
  std::vector<float> totalRecoSpectrum_OOT;
  for(int i=0; i<samples->size(); i++){
    totalRecoSpectrum_OOT.push_back(0);
  }
  
  
  //  for(int iBx=0; iBx<3; iBx++){
  for(int iBx=0; iBx<samplesReco->size(); iBx++){
    std::cout << " iBx = " << iBx << " :: " << samplesReco->size() << std::endl;
    std::cout << " Energy = " << samplesReco->at(iBx) << std::endl;
    grPulseReco[iBx] = new TGraph();
    for(int i=0; i<samples->size(); i++){
      //    std::cout << "  >> i = " << i << std::endl;
      grPulseReco[iBx]->SetPoint(i, i * NFREQ + activeBXs->at(iBx)*NFREQ + 2 * 25, pulseShapeTemplate->at(i) * samplesReco->at(iBx));
      
      int iReco = (i * NFREQ + activeBXs->at(iBx)*NFREQ + 2 * 25) / NFREQ;
      if ( iReco >= 0 && iReco <samples->size() ) {
        totalRecoSpectrum.at(iReco) += pulseShapeTemplate->at(i) * samplesReco->at(iBx);
        if (iBx != 5) totalRecoSpectrum_OOT.at(iReco) += pulseShapeTemplate->at(i) * samplesReco->at(iBx);
      } 
      
    }
    grPulseReco[iBx]->SetMarkerColor(color[iBx]);
    grPulseReco[iBx]->SetLineColor(color[iBx]);
    grPulseReco[iBx]->SetMarkerSize(1);
    grPulseReco[iBx]->SetMarkerStyle(21+iBx);
    TString nameHistoTitle = Form ("BX %d", activeBXs->at(iBx));
    leg->AddEntry(grPulseReco[iBx],nameHistoTitle.Data(),"p");
  }
  
  for(int i=0; i<samples->size(); i++){
    totalRecoSpectrum.at(i) += best_pedestal;
    totalRecoSpectrum_OOT.at(i) += 0;
  }
    
  grPulse->Draw("ALP");
  grPulse->GetYaxis()->SetRangeUser(0 > (best_pedestal-1) ? (best_pedestal-1)  : 0 , samplesReco->at(5)+10+best_pedestal);
  
  //  for(int iBx=0; iBx<3; iBx++){
  for(int iBx=0; iBx<samplesReco->size(); iBx++){
    grPulseReco[iBx]->Draw("PL");
  }
  
  for(int i=0; i<samples->size(); i++){
    grPulseRecoAll->SetPoint(i, i * NFREQ, totalRecoSpectrum.at(i));
    grPulseReco_OOT_All->SetPoint(i, i * NFREQ, totalRecoSpectrum_OOT.at(i));
  }
  
  grPulseRecoAll->SetMarkerColor(kMagenta);
  grPulseRecoAll->SetLineColor(kMagenta);
  grPulseRecoAll->SetLineStyle(1);
  grPulseRecoAll->SetMarkerSize(2);
  grPulseRecoAll->SetMarkerStyle(24);
  grPulseRecoAll->Draw("PL");

  grPulseReco_OOT_All->SetMarkerColor(kMagenta);
  grPulseReco_OOT_All->SetLineColor(kMagenta);
  grPulseReco_OOT_All->SetLineStyle(2);
  grPulseReco_OOT_All->SetMarkerSize(2);
  grPulseReco_OOT_All->SetMarkerStyle(27);
  grPulseReco_OOT_All->Draw("PL");
  
  
  grPulseReco_OOT_All->GetXaxis()->SetTitle("time [ns]");
  grPulseReco_OOT_All->GetYaxis()->SetTitle("Out of time [GeV]");
    
  
  grPulse->GetXaxis()->SetTitle("time [ns]");
 
  
  TGraph *grPedestal = new TGraph();
  for(int i=0; i<samples->size(); i++){
    grPedestal->SetPoint(i, i * NFREQ , best_pedestal);
  }
  grPedestal->SetMarkerColor(kBlue);
  grPedestal->SetLineColor(kBlue);
  grPedestal->SetLineWidth(4);
  grPedestal->SetLineStyle(3);
  grPedestal->SetMarkerSize(1);
  grPedestal->Draw("PL");
  
  
  grPulse_noise->Draw("PL");
  
  grPU->Draw("L");
  
  leg->Draw();
  
  
   
  
  
  ccPulseAndReco->cd(2);
  grPulseReco_OOT_All->Draw("APL");
  grPU->Draw("L");
  
  gPad->SetGrid();
  
  ccPulseAndReco->cd(1);
  
  
  std::cout << " done " << std::endl;
  
  std::cout << " complete_pedestal = " << complete_pedestal->size() << std::endl;
  std::cout << " complete_chi2 = " << complete_chi2->size() << std::endl;
  
  TCanvas* ccScan = new TCanvas ("ccScan", "", 800, 600);
  
  TGraph* chi2_scan_graph =  new TGraph();
  for (int i=0; i<complete_pedestal->size(); i++) {
    chi2_scan_graph->SetPoint(i, complete_pedestal->at(i), complete_chi2->at(i));
  }
  
//   TF1* fitParabola = new TF1("fitParabola", "[1] * (x-[0])*(x-[0]) + [2]", -0.5, 0.5);
//   
//   chi2_scan_graph->Fit("fitParabola", "R");
//   
//   std::cout << " best_pedestal = " << best_pedestal << std::endl;
//   std::cout << " best_pedestal[fit] = " << fitParabola->GetParameter(0) << std::endl;
//   
//   std::cout << " return_chi2 = " << return_chi2 << std::endl;
//   std::cout << " return_chi2[fit] = " << fitParabola->GetParameter(2) << std::endl;
  
  chi2_scan_graph->SetMarkerColor(kMagenta);
  chi2_scan_graph->SetLineColor(kMagenta);
  chi2_scan_graph->SetLineStyle(2);
  chi2_scan_graph->SetMarkerSize(2);
  chi2_scan_graph->SetMarkerStyle(25);
  chi2_scan_graph->Draw("APL");
  chi2_scan_graph->GetXaxis()->SetTitle("pedestal [GeV]");
  chi2_scan_graph->GetYaxis()->SetTitle("#chi^{2}");
  
  ccScan->SetGrid();
  
  
  
  
  
  TCanvas* ccPedestals = new TCanvas ("ccPedestals", "", 800, 600);
  
  TH1F* histo_pedestal = new TH1F("histo_pedestal", "fitted pedestal", 100, -2, 2);
  tree->Draw("best_pedestal >> histo_pedestal", "", "goff");
  
  histo_pedestal->SetLineColor(kBlue);
  histo_pedestal->SetFillColor(kBlue);
  histo_pedestal->SetFillStyle(3001);
  histo_pedestal->Draw();
  histo_pedestal->GetXaxis()->SetTitle("pedestal [GeV]");
  
  
  ccPedestals->SetGrid();
  
  
  
  
  
  
  
  
  TCanvas* ccSimple = new TCanvas ("ccSimple", "Simple", 800, 1200);
  
  ccSimple->Divide(1,2);
  ccSimple->cd(1);
  
  TGraph *simple_grPulseRecoAll = new TGraph();
  TGraph *simple_grPulseReco_OOT_All = new TGraph();
  TGraph *simple_grPulseReco[samplesReco->size()];
 
  std::vector<float> simple_totalRecoSpectrum;
  for(int i=0; i<samples->size(); i++){
    simple_totalRecoSpectrum.push_back(0);
  }
  
  std::vector<float> simple_totalRecoSpectrum_OOT;
  for(int i=0; i<samples->size(); i++){
    simple_totalRecoSpectrum_OOT.push_back(0);
  }
  
  
  //  for(int iBx=0; iBx<3; iBx++){
  for(int iBx=0; iBx<samplesReco->size(); iBx++){
    std::cout << " iBx = " << iBx << " :: " << samplesReco->size() << std::endl;
    std::cout << " Energy = " << (complete_samplesReco->at( complete_pedestal->size()/2 )).at(iBx) << std::endl;
    simple_grPulseReco[iBx] = new TGraph();
    for(int i=0; i<samples->size(); i++){
      //    std::cout << "  >> i = " << i << std::endl;
      simple_grPulseReco[iBx]->SetPoint(i, i * NFREQ + activeBXs->at(iBx)*NFREQ + 2 * 25, pulseShapeTemplate->at(i) * ( (complete_samplesReco->at( complete_pedestal->size()/2 )).at(iBx) ));
      
      int iReco = (i * NFREQ + activeBXs->at(iBx)*NFREQ + 2 * 25) / NFREQ;
      if ( iReco >= 0 && iReco <samples->size() ) {
        simple_totalRecoSpectrum.at(iReco) += pulseShapeTemplate->at(i) * (complete_samplesReco->at( complete_pedestal->size()/2 )).at(iBx);
        if (iBx != 5) simple_totalRecoSpectrum_OOT.at(iReco) += pulseShapeTemplate->at(i) * (complete_samplesReco->at( complete_pedestal->size()/2 )).at(iBx);
      } 
      
    }
    simple_grPulseReco[iBx]->SetMarkerColor(color[iBx]);
    simple_grPulseReco[iBx]->SetLineColor(color[iBx]);
    simple_grPulseReco[iBx]->SetMarkerSize(1);
    simple_grPulseReco[iBx]->SetMarkerStyle(21+iBx);
    TString nameHistoTitle = Form ("BX %d", activeBXs->at(iBx));
//     leg->AddEntry(simple_grPulseReco[iBx],nameHistoTitle.Data(),"p");
  }
  
  
  //---- no additional pedestals here
  for(int i=0; i<samples->size(); i++){
    simple_totalRecoSpectrum.at(i) += 0;
    simple_totalRecoSpectrum_OOT.at(i) += 0;
  }
  
  grPulse->Draw("ALP");
  grPulse->GetYaxis()->SetRangeUser(0 > (best_pedestal-1) ? (best_pedestal-1)  : 0 , samplesReco->at(5)+10+best_pedestal);
  
  //  for(int iBx=0; iBx<3; iBx++){
  for(int iBx=0; iBx<samplesReco->size(); iBx++){
    simple_grPulseReco[iBx]->Draw("PL");
  }
  
  for(int i=0; i<samples->size(); i++){
    simple_grPulseRecoAll->SetPoint(i, i * NFREQ, simple_totalRecoSpectrum.at(i));
    simple_grPulseReco_OOT_All->SetPoint(i, i * NFREQ, simple_totalRecoSpectrum_OOT.at(i));
  }
  
  simple_grPulseRecoAll->SetMarkerColor(kMagenta);
  simple_grPulseRecoAll->SetLineColor(kMagenta);
  simple_grPulseRecoAll->SetLineStyle(1);
  simple_grPulseRecoAll->SetMarkerSize(2);
  simple_grPulseRecoAll->SetMarkerStyle(24);
  simple_grPulseRecoAll->Draw("PL");
  
  
  simple_grPulseReco_OOT_All->SetMarkerColor(kMagenta);
  simple_grPulseReco_OOT_All->SetLineColor(kMagenta);
  simple_grPulseReco_OOT_All->SetLineStyle(1);
  simple_grPulseReco_OOT_All->SetMarkerSize(2);
  simple_grPulseReco_OOT_All->SetMarkerStyle(27);
  simple_grPulseReco_OOT_All->Draw("PL");
  
  simple_grPulseReco_OOT_All->GetXaxis()->SetTitle("time [ns]");
  simple_grPulseReco_OOT_All->GetYaxis()->SetTitle("Out of time [GeV]");
  
  
  
  
  grPulse->GetXaxis()->SetTitle("time [ns]");
  grPulse_noise->Draw("PL");
  
  grPU->Draw("L");
  
  leg->Draw();
    
  
  ccSimple->cd(2);
  
  simple_grPulseReco_OOT_All->Draw("APL");
  grPU->Draw("L");
  
  gPad->SetGrid();
  
  ccSimple->cd(1);
  
  
  
  
  
  TCanvas* ccRecoSimple = new TCanvas ("ccRecoSimple","Reco Simple",800,600);
  ccRecoSimple->SetRightMargin(0.3);
  TGraph *grRecoSimple = new TGraph();
  for(int i=0; i<samplesReco->size(); i++){
    grRecoSimple->SetPoint(i, activeBXs->at(i), (complete_samplesReco->at( complete_pedestal->size()/2 )).at(i));
  }
  grRecoSimple->SetMarkerSize(1);
  grRecoSimple->SetMarkerStyle(24);
  grRecoSimple->SetMarkerColor(kRed);
  grRecoSimple->Draw("ALP");
  grRecoSimple->GetXaxis()->SetTitle("BX");
  grReco->Draw("LP");
  
  TLegend* leg2 = new TLegend(0.71,0.10,0.99,0.90);
  leg2->AddEntry(grReco,      "Advanced Multifit","p");
  leg2->AddEntry(grRecoSimple,"Multifit","p");
  leg2->Draw();
  
  ccRecoSimple->SetGrid();
  
  
  
}






