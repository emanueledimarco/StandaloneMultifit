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

void plotPulse (std::string nameInputFile = "output.root", int nEvent = 10){

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
 
 TTree* tree = (TTree*) file->Get("Samples");
  
 int    nWF;
 std::vector<double>* pulse_signal    = new std::vector<double>;
 std::vector<double>* samples     = new std::vector<double>;
 std::vector<double>* samples_noise = new std::vector<double>;
 
 float NFREQ;
 
 tree->SetBranchAddress("nWF",      &nWF);
 tree->SetBranchAddress("pulse_signal", &pulse_signal);
 tree->SetBranchAddress("samples",   &samples);
 tree->SetBranchAddress("samples_noise",   &samples_noise);
 tree->SetBranchAddress("nFreq",   &NFREQ);
 
 tree->GetEntry(nEvent);

 TCanvas* ccPulse = new TCanvas ("ccPulse","",800,600);
 ccPulse->SetGrid();

 TGraph *grPulse_signal = new TGraph();
 for(int i=0; i<nWF/2; i++){
   grPulse_signal->SetPoint(i, i/4., pulse_signal->at(nWF/2+i));
 }

 grPulse_signal->SetMarkerSize(0.2);
 grPulse_signal->SetMarkerStyle(kFullCircle);
 grPulse_signal->SetMarkerColor(kRed);
 grPulse_signal->SetLineColor(kRed);
 grPulse_signal->SetLineWidth(1);
 
 TGraph *grPulse_noise = new TGraph();
 for(int i=0; i<(int)samples->size(); i++){
   grPulse_noise->SetPoint(i, i * NFREQ, samples_noise->at(i));
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

 grPulse_signal->GetXaxis()->SetTitle("time [ns]");
 grPulse_signal->GetXaxis()->SetRangeUser(0, (samples->size()-1) * NFREQ);
 grPulse_signal->GetYaxis()->SetRangeUser(-0.5, TMath::MaxElement(grPulse->GetN(),grPulse->GetY())*1.2);
 grPulse_signal->Draw("APL");
 grPulse_noise->Draw("PL");
 grPulse->Draw("LP");
 
 
 TLegend* leg = new TLegend(0.91,0.90,0.99,0.99);

 leg->AddEntry(grPulse_signal,"signal input","p");
 leg->AddEntry(grPulse_noise,"noise input","p");
 leg->AddEntry(grPulse,"total","p");
 leg->Draw();
 
 ccPulse->SaveAs("ccpulse.pdf");
 
}


