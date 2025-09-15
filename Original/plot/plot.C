//---- plot output of multifit


#include <iostream>
#include <fstream>


void plot (std::string nameInputFile = "output.root",  float valueToWrite = 0){
 
 Color_t* color = new Color_t [200];
 color[0] = kRed ;
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
 
 TCanvas* cc = new TCanvas ("cc","",800,600);
 
 TH1F* histo[200];
 
 TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
 std::vector<double>* samplesReco = new std::vector<double>;
 double amplitudeTruth;
 tree->SetBranchAddress("samplesReco",             &samplesReco);
 tree->SetBranchAddress("amplitudeTruth",          &amplitudeTruth);
 tree->GetEntry(0);
 int totBX = samplesReco->size();
 
 for (int i=0; i<totBX; i++) {
  TString nameHisto = Form ("samplesReco_%d", i);
  TString nameHistoTitle = Form ("BX %d", i-5);
//   TString nameHistoTitle = Form ("BX %d", i);
  histo[i] = new TH1F (nameHisto.Data(), nameHistoTitle.Data(), 3000, 0, 3*amplitudeTruth);
  
  TString toDraw = Form ("samplesReco[%d] >> %s", i, nameHisto.Data());
  tree->Draw(toDraw.Data());
  histo[i] -> SetLineColor(color[i]);
  histo[i] -> SetLineWidth(2);
  leg->AddEntry(histo[i],nameHistoTitle.Data(),"l");
 }

 for (int i=0; i<totBX; i++) {
  if (i==0) {
   histo[i]->Draw();
   histo[i]->GetXaxis()->SetTitle("Energy [GeV]");
  }
  else {
   histo[i]->Draw("same");   
  }
 }
 
 leg->Draw();
 
 TCanvas* ccGood = new TCanvas ("ccGood","",800,600);
 for (int i=0; i<totBX; i++) {
  if (histo[i]->GetMean() > 0.5*amplitudeTruth) {
   histo[i]->Draw();
   histo[i]->GetXaxis()->SetTitle("Energy [GeV]");
  }
 }
 

 
 TCanvas* ccGoodBx5 = new TCanvas ("ccGoodBx5","BX 5",800,600);
 for (int i=0; i<totBX; i++) {
   if ( i == 5)  {
     histo[i]->Draw();
     histo[i]->GetXaxis()->SetTitle("Energy [GeV]");
   }
 }
 
 
 
 
 std::ofstream myfile;
 myfile.open ("results.txt", std::ofstream::out | std::ofstream::app);
 for (int i=0; i<totBX; i++) {
   if (histo[i]->GetMean() > 0.5*amplitudeTruth) {
     myfile << "  " <<  valueToWrite << "  " <<  histo[i]->GetMean() << std::endl ;
   }
  }
 myfile.close(); 
 
 
 
 TCanvas* ccPedestal = new TCanvas ("ccPedestal","",800,600);
 TH1F* histo_pedestal = new TH1F("histo_pedestal", "pedestal", 100, -3, 3);
 
 if (tree->GetListOfBranches()->FindObject("best_pedestal")) {
   TString toDraw = Form ("best_pedestal >> histo_pedestal");
   tree->Draw(toDraw.Data());
   histo_pedestal -> SetLineColor(kRed);
   histo_pedestal -> SetLineWidth(2);
   histo_pedestal->Draw();
   histo_pedestal->GetXaxis()->SetTitle("Pedestal");
 }
 
 
}