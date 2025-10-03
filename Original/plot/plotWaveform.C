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

  std::cout << "ccc = " << std::endl;
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
 std::cout << " NFREQ = " << NFREQ << std::endl;
 std::cout << " nWF = " << nWF << std::endl;
 
 TCanvas* ccpulse_signal = new TCanvas ("ccpulse_signal","",800,600);
 TGraph *gr = new TGraph();
 TGraph *gr_sampling = new TGraph();
 for(int i=0; i<nWF; i++){
  gr->SetPoint(i, i/4., pulse_signal->at(i));
 }
 gr->Draw("AL"); 
 gr->SetLineColor(kMagenta);
 gr->SetLineWidth(2);
 gr->GetXaxis()->SetTitle("time [ns]");
 
 gr_sampling->Draw("P");
 gr_sampling->SetMarkerStyle(4);
 gr_sampling->SetMarkerColor(kBlue);
 
 ccpulse_signal->SetGrid();
 
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
 grPulse->GetXaxis()->SetTitle("time [ns]");
 
 grPulse->Draw("ALP");
 
 grPulse->GetXaxis()->SetTitle("time [ns]");
 
 grPulse_noise->Draw("PL");
 
 TLegend* leg = new TLegend(0.91,0.10,0.99,0.90);

 leg->AddEntry(grPulse_noise,"noise input","p");
 
 
 
 leg->Draw();
 
 
 
 
}


