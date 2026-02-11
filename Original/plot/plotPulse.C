// on the tree output of multifit.cc, run plotPulse("test.root","RecoAndSim",3)

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TCanvas.h>

#include <iostream>
#include <iomanip>

//---- plot output of multifit

void plotPulse (std::string nameInputFile = "output.root", std::string treeName="Samples", int nEvent = 10){

 Color_t* color = new Color_t [200];
 color[0] = kAzure; //kRed ;
 color[1] = kAzure + 10 ;
 color[2] = kYellow + 2 ;
 color[3] = kGreen ;
 color[4] = kRed ; // kGreen + 4 ;
 color[5] = kBlue ;
 color[6] = kCyan ;
 color[7] = kPink + 1 ;
 color[8] = kBlack ;
 color[9] = kYellow + 4 ;
 for (int i=0; i<30; i++) {
  color[i+10] = kBlue + i;
 }
 
 
 TFile *file = new TFile(nameInputFile.c_str());
 
 TTree* tree = (TTree*) file->Get(treeName.c_str()); // for the file after multifit it should be "RecoAndSim"
  
 int    nWF;
 int    nTemplateBins;
 std::vector<double>* pulse_signal     = new std::vector<double>;
 std::vector<double>* pileup_signal    = new std::vector<double>;
 std::vector<double>* samples     = new std::vector<double>;
 std::vector<double>* samples_noise = new std::vector<double>;
 std::vector<double>* samplesReco = new std::vector<double>;
 std::vector<double>* pedestalsReco = new std::vector<double>;
 std::vector<int>*    activeBXs   = new std::vector<int>;
 float pulseShapeTemplate[9];
 
 float NFREQ;
 
 tree->SetBranchAddress("nWF",      &nWF);
 tree->SetBranchAddress("pulse_signal", &pulse_signal);
 tree->SetBranchAddress("pileup_signal", &pileup_signal);
 tree->SetBranchAddress("samples",   &samples);
 tree->SetBranchAddress("samples_noise",   &samples_noise);
 tree->SetBranchAddress("nFreq",   &NFREQ);
 tree->SetBranchAddress("samplesReco", &samplesReco);
 tree->SetBranchAddress("activeBXs", &activeBXs);
 tree->SetBranchAddress("nTemplateBins", &nTemplateBins);
 tree->SetBranchAddress("pulseShapeTemplate",   pulseShapeTemplate);
 tree->SetBranchAddress("pedestalsReco", &pedestalsReco);

 tree->GetEntry(nEvent);

 std::cout << "Got event " << nEvent << std::endl;

 TCanvas* ccPulse = new TCanvas ("ccPulse","",800,600);
 ccPulse->SetGrid();

 TGraph *grPulse_signal = new TGraph();
 for(int i=0; i<nWF/2; i++){
   grPulse_signal->SetPoint(i, i/4., pulse_signal->at(nWF/2+i));
 }

 grPulse_signal->SetMarkerSize(0.4);
 grPulse_signal->SetMarkerStyle(kFullCircle);
 grPulse_signal->SetMarkerColor(kRed);
 grPulse_signal->SetLineColor(kRed);
 grPulse_signal->SetLineWidth(1);

 TGraph *grPulse_pileup = new TGraph();
 for(int i=0; i<nWF/2; i++){
   grPulse_pileup->SetPoint(i, i/4., pileup_signal->at(nWF/2+i));
 }

 grPulse_pileup->SetMarkerSize(0.4);
 grPulse_pileup->SetMarkerStyle(kFullDiamond);
 grPulse_pileup->SetMarkerColor(kOrange+4);
 grPulse_pileup->SetLineColor(kOrange+4);
 grPulse_pileup->SetLineWidth(1);

 TGraph *grPulse_noise = new TGraph();
 for(int i=0; i<(int)samples->size(); i++){
   grPulse_noise->SetPoint(i, i * NFREQ, samples_noise->at(i) + pedestalsReco[0].at(0));
 }

 grPulse_noise->SetMarkerSize(1);
 grPulse_noise->SetMarkerStyle(kFullCircle);
 grPulse_noise->SetMarkerColor(kGray);
 grPulse_noise->SetLineColor(kGray);
 grPulse_noise->SetLineStyle(3);
  
 
 TGraph *grPulse = new TGraph();
 for(int i=0; i<(int)samples->size(); i++){
   grPulse->SetPoint(i, i * NFREQ, samples->at(i));
 }

 grPulse->SetMarkerSize(1);
 grPulse->SetMarkerStyle(kFullCircle);
 grPulse->SetMarkerColor(kBlack);
 grPulse->SetLineStyle(3);
 grPulse->SetLineColor(kBlack);
 grPulse->SetLineWidth(2);

 TGraph* pedestalGraph = new TGraph();
 TGraph *grPulse_sum = new TGraph();
 for(int i=0; i<(int)samples->size(); i++){
   grPulse_sum->SetPoint(i, i * NFREQ, grPulse_pileup->Eval(grPulse_noise->GetX()[i])+grPulse_signal->Eval(grPulse_noise->GetX()[i]) + pedestalsReco[0].at(0) );
   pedestalGraph->SetPoint(i, i*NFREQ, pedestalsReco[0].at(0));
 }

 grPulse_sum->SetMarkerSize(1);
 grPulse_sum->SetMarkerStyle(kFullSquare);
 grPulse_sum->SetMarkerColor(kViolet);
 grPulse_sum->SetLineStyle(3);
 grPulse_sum->SetLineColor(kViolet);
 grPulse_sum->SetLineWidth(2);


 grPulse_signal->GetXaxis()->SetTitle("time [ns]");
 grPulse_signal->GetYaxis()->SetTitle("amplitude [GeV]");
 grPulse_signal->GetXaxis()->SetRangeUser(0, (samples->size()-1) * NFREQ);
 grPulse_signal->GetYaxis()->SetRangeUser(-0.5, TMath::MaxElement(grPulse->GetN(),grPulse->GetY())*1.2);
 grPulse_signal->Draw("APL");
 grPulse_pileup->Draw("PL");
 grPulse_noise->Draw("PL");
 grPulse->Draw("LP");
 grPulse_sum->Draw("LP");
 pedestalGraph->Draw("PL");
 pedestalGraph->SetLineColor(kBlue);
 pedestalGraph->SetMarkerColor(kBlue);

 TLegend* leg = new TLegend(0.91,0.10,0.99,0.90);

 leg->AddEntry(grPulse_signal,"signal","p");
 leg->AddEntry(grPulse_pileup,"PU sum","pl");
 leg->AddEntry(grPulse_noise,"noise","p");
 leg->AddEntry(grPulse,"DIGIs","p");
 leg->Draw();
 
 ccPulse->SaveAs("ccpulse.pdf");



 // === Now plot reconstruction output ===
 
 TCanvas* ccPulseAndReco = new TCanvas ("ccPulseAndReco","",800,600);
 TGraph *grPulseRecoAll = new TGraph();
 const int nReco = (int)samplesReco->size();
 std::vector<TGraph*> grPulseReco;
 /* std::cout << " samplesReco->size() = " << samplesReco->size() << std::endl; */
 /* std::cout << " activeBXs->size() = " << activeBXs->size() << std::endl; */
 /* std::cout << " nTemplateBins = " << nTemplateBins << std::endl; */
 /* std::cout << " samples->size() = " << samples->size() << std::endl; */
 /* std::cout << " NFREQ = " << NFREQ << std::endl; */
 
 TLegend* leg2 = new TLegend(0.91,0.10,0.99,0.90);
 
 std::vector<float> totalRecoSpectrum;
 for(int i=0; i<(int)samples->size(); i++){
   totalRecoSpectrum.push_back(0);
 }
 
 for(int iBx=0; iBx<(int)samplesReco->size(); iBx++){
  std::cout << " iBx = " << iBx << std::endl;
  std::cout << " Energy = " << samplesReco->at(iBx) << std::endl;
  grPulseReco.push_back(new TGraph());
  for(int i=0; i<(int)samples->size(); i++){
    float templateVal = i < 9 ? pulseShapeTemplate[i] : 0;
    grPulseReco[iBx]->SetPoint(i, i*NFREQ + activeBXs->at(iBx)*25 + 6*NFREQ, templateVal * samplesReco->at(iBx));
    int iReco = i + activeBXs->at(iBx) * int(25./NFREQ) + 6;
    if(iReco >= 0 && iReco < (int)samples->size()) {
      totalRecoSpectrum.at(iReco) += templateVal * samplesReco->at(iBx);
    }
  }
  grPulseReco[iBx]->SetMarkerColor(color[iBx]);
  grPulseReco[iBx]->SetLineColor(color[iBx]);
  grPulseReco[iBx]->SetMarkerSize(1);
  grPulseReco[iBx]->SetMarkerStyle(21+iBx);
  TString nameHistoTitle = Form ("BX %d", activeBXs->at(iBx));
  leg2->AddEntry(grPulseReco[iBx],nameHistoTitle.Data(),"p");
 }

  for(int i=0; i<(int)samples->size(); i++){
   grPulseRecoAll->SetPoint(i, i*NFREQ, totalRecoSpectrum.at(i) + pedestalsReco[0].at(0));
 }

 grPulseRecoAll->SetMarkerColor(kBlack);
 grPulseRecoAll->SetLineColor(kBlack);
 grPulseRecoAll->SetLineStyle(1);
 grPulseRecoAll->SetMarkerSize(1.5);
 grPulseRecoAll->SetMarkerStyle(kFullCircle);
 leg2->AddEntry(grPulseRecoAll,"Total","p");
 grPulseRecoAll->GetYaxis()->SetTitle("amplitude [GeV]");
 grPulseRecoAll->GetXaxis()->SetTitle("time [ns]");
 grPulseRecoAll->GetXaxis()->SetRangeUser(0, (samples->size()-1) * NFREQ);
 grPulseRecoAll->GetYaxis()->SetRangeUser(-0.5, TMath::MaxElement(grPulseRecoAll->GetN(),grPulseRecoAll->GetY())*1.2);

 grPulseRecoAll->Draw("APC");
 grPulse_noise->Draw("PC");
 for(int iBx=0; iBx<(int)samplesReco->size(); iBx++){
   grPulseReco[iBx]->Draw("PC");
 }

 grPulse->SetMarkerStyle(kFullTriangleUp);
 grPulse->SetMarkerColor(kRed);
 grPulse->Draw("PL");
 pedestalGraph->Draw("PL");
 pedestalGraph->SetLineColor(kBlue);
 pedestalGraph->SetMarkerColor(kBlue);
 leg2->AddEntry(grPulse,"Data","p");
 leg2->AddEntry(grPulse_noise,"Noise","p");
 leg2->AddEntry(pedestalGraph,"Reco pedestal","l");
 
 leg2->Draw();

 ccPulseAndReco->SaveAs("ccpulseAndReco.pdf"); 
 
}


