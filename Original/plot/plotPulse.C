//---- plot output of multifit

void plotPulse (std::string nameInputFile = "output.root", int nEvent = 10, float pedestal = 0.){
 
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
 std::vector<double>* samplesReco = new std::vector<double>;
 std::vector<double>* samples     = new std::vector<double>;
 std::vector<double>* samples_noise = new std::vector<double>;
 std::vector<int>*    activeBXs   = new std::vector<int>;
 std::vector<double>* pulseShapeTemplate     = new std::vector<double>;
 float best_pedestal = 0;
 
 float NFREQ;
 
 tree->SetBranchAddress("nWF",      &nWF);
 tree->SetBranchAddress("pulse_signal", &pulse_signal);
 tree->SetBranchAddress("samplesReco", &samplesReco);
 tree->SetBranchAddress("samples",   &samples);
 tree->SetBranchAddress("samples_noise",   &samples_noise);
 tree->SetBranchAddress("activeBXs", &activeBXs);
 tree->SetBranchAddress("nFreq",   &NFREQ);
 tree->SetBranchAddress("pulseShapeTemplate",   &pulseShapeTemplate);
 
 if (tree->GetListOfBranches()->FindObject("best_pedestal")) {
   tree->SetBranchAddress("best_pedestal", &best_pedestal);
 }
 
 float input_pedestal = 0;
 if (tree->GetListOfBranches()->FindObject("input_pedestal")) {
   tree->SetBranchAddress("input_pedestal", &input_pedestal);
 }
 
 
 
 
 tree->GetEntry(nEvent);
 std::cout << " NFREQ = " << NFREQ << std::endl;
 
 int shift_my_test = 13*4;
 
 TCanvas* ccpulse_signal = new TCanvas ("ccpulse_signal","",800,600);
 TGraph *gr = new TGraph();
 TGraph *gr_sampling = new TGraph();
 for(int i=0; i<nWF; i++){
  gr->SetPoint(i, i/4., pulse_signal->at(i));
  
  if (!((i+shift_my_test)%(4*25))) {
    gr_sampling->SetPoint((i+shift_my_test)/(4*25), i/4., pulse_signal->at(i));
  }
  
 }
 gr->Draw("AL"); 
 gr->SetLineColor(kMagenta);
 gr->SetLineWidth(2);
 gr->GetXaxis()->SetTitle("time [ns]");
 
 gr_sampling->Draw("P");
 gr_sampling->SetMarkerStyle(4);
 gr_sampling->SetMarkerColor(kBlue);
 
 ccpulse_signal->SetGrid();
 
 
 
 TCanvas* ccReco = new TCanvas ("ccReco","",800,600);
 ccReco->SetRightMargin(0.3);
 TGraph *grReco = new TGraph();
 for(int i=0; i<samplesReco->size(); i++){
  std::cout << " i, activeBXs->at(i), samplesReco->at(i) = " << i << "::" << samplesReco->size() << " -> " << activeBXs->at(i) << " , " << samplesReco->at(i) << std::endl;
  grReco->SetPoint(i, activeBXs->at(i), samplesReco->at(i));
 }
 grReco->SetMarkerSize(2);
 grReco->SetMarkerStyle(22);
 grReco->SetMarkerColor(kBlue);
 grReco->Draw("ALP");
 grReco->GetXaxis()->SetTitle("BX");
 
 
 TCanvas* ccPulse = new TCanvas ("ccPulse","",800,600);
 
 TGraph *grPulse_noise = new TGraph();
 for(int i=0; i<samples->size(); i++){
   grPulse_noise->SetPoint(i, i * NFREQ , samples_noise->at(i));
 }
 grPulse_noise->SetMarkerSize(2);
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
 grPulse->SetLineStyle(3);
 grPulse->SetLineColor(kRed);
 grPulse->SetLineWidth(2);
 grPulse->Draw("ALP");
 grPulse->GetYaxis()->SetRangeUser(0 > (input_pedestal-1) ? (input_pedestal-1)  : 0 , TMath::MaxElement(grPulse->GetN(),grPulse->GetY()) +2);
 grPulse->GetXaxis()->SetTitle("time [ns]");
 
 std::cout << " grPulse->GetMaximum() = " << grPulse->GetMaximum() << std::endl;
 std::cout << "  TMath::MaxElement(grPulse->GetN(),grPulse->GetY())  = " <<  TMath::MaxElement(grPulse->GetN(),grPulse->GetY())  << std::endl;
 
 
 
 TGraph *grPulse_pedestal = new TGraph();
 for(int i=0; i<samples->size(); i++){
   grPulse_pedestal->SetPoint(i, i * NFREQ , input_pedestal);
 }
 grPulse_pedestal->SetMarkerSize(1.5);
 grPulse_pedestal->SetMarkerStyle(21);
 grPulse_pedestal->SetMarkerColor(kRed);
 grPulse_pedestal->SetLineStyle(3);
 grPulse_pedestal->SetLineColor(kRed);
 grPulse_pedestal->SetLineWidth(1);
 grPulse_pedestal->Draw("LP");
 grPulse_pedestal->GetXaxis()->SetTitle("time [ns]");
 
 
 
 TGraph *grPulseReco_pedestal = new TGraph();
 for(int i=0; i<samples->size(); i++){
   grPulseReco_pedestal->SetPoint(i, i * NFREQ , best_pedestal);
 }
 grPulseReco_pedestal->SetMarkerSize(1.5);
 grPulseReco_pedestal->SetMarkerStyle(21);
 grPulseReco_pedestal->SetMarkerColor(kRed);
 grPulseReco_pedestal->SetLineStyle(3);
 grPulseReco_pedestal->SetLineColor(kRed);
 grPulseReco_pedestal->SetLineWidth(1);
 grPulseReco_pedestal->Draw("LP");
 grPulseReco_pedestal->GetXaxis()->SetTitle("time [ns]");
 
 
 grPulseReco_pedestal->SetMarkerColor(kMagenta);
 grPulseReco_pedestal->SetLineColor(kMagenta);
 grPulseReco_pedestal->SetLineStyle(1);
 grPulseReco_pedestal->SetMarkerSize(1);
 grPulseReco_pedestal->SetMarkerStyle(24);
 
 
 std::cout << " best_pedestal = " << best_pedestal << std::endl;
 
 
 std::cout << " end " << std::endl;
 
 
 TCanvas* ccPulseAndReco = new TCanvas ("ccPulseAndReco","",800,600);
 TGraph *grPulseRecoAll = new TGraph();
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
 
 
//  for(int iBx=0; iBx<3; iBx++){
 for(int iBx=0; iBx<samplesReco->size(); iBx++){
  std::cout << " iBx = " << iBx << std::endl;
  std::cout << " Energy = " << samplesReco->at(iBx) << std::endl;
  grPulseReco[iBx] = new TGraph();
  for(int i=0; i<samples->size(); i++){
//    std::cout << "  >> i = " << i << std::endl;
   grPulseReco[iBx]->SetPoint(i, i * NFREQ + activeBXs->at(iBx)*NFREQ + 2 * 25, pulseShapeTemplate->at(i) * samplesReco->at(iBx));
  
   int iReco = (i * NFREQ + activeBXs->at(iBx)*NFREQ + 2 * 25) / NFREQ;
   if ( iReco >= 0 && iReco <samples->size() ) {
    totalRecoSpectrum.at(iReco) += pulseShapeTemplate->at(i) * samplesReco->at(iBx);
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
   totalRecoSpectrum.at(i) += pedestal;
 }
 
 
 for(int i=0; i<samples->size(); i++){
   totalRecoSpectrum.at(i) += best_pedestal;
 }
 
 
 
 grPulse->Draw("ALP");
 grPulse_pedestal->Draw("LP");
 
//  for(int iBx=0; iBx<3; iBx++){
 for(int iBx=0; iBx<samplesReco->size(); iBx++){
  grPulseReco[iBx]->Draw("PL");
 }
 
 for(int i=0; i<samples->size(); i++){
  grPulseRecoAll->SetPoint(i, i * NFREQ, totalRecoSpectrum.at(i));
 }
 
 grPulseRecoAll->SetMarkerColor(kMagenta);
 grPulseRecoAll->SetLineColor(kMagenta);
 grPulseRecoAll->SetLineStyle(1);
 grPulseRecoAll->SetMarkerSize(2);
 grPulseRecoAll->SetMarkerStyle(24);
 grPulseRecoAll->Draw("PL");
 grPulseReco_pedestal->Draw("PL");
 grPulse->GetXaxis()->SetTitle("time [ns]");
 
 grPulse_noise->Draw("PL");
 
 
 leg->AddEntry(grPulse_noise,"noise input","p");
 
 
 
 
 
 leg->Draw();
 
 
 
 
}


