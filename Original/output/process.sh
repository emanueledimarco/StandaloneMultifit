
echo "shift,tmean,tmeanerr,trms,trmserr,emean,emeanerr,erms,ermserr" > $1

for i in $(seq -5 6); do
root -l -b  shift_scans/output_-_noisescale1.00_fitPed1_globalshift_${i}.root <<EOF
  RecoAndSim->Draw("timeReco[3]>>ht(2000, -1, 1)");
  RecoAndSim->Draw("samplesReco[3]/signalTruth>>he(1000, 0.95, 1.05)");
  ht->SetName("time_${i}");
  he->SetName("en_${i}");
  ht->SetTitle("shift: ${i} ns");
  he->SetTitle("shift: ${i} ns");

  double f_logn(Double_t* x, Double_t* par) {
          Double_t eta = par[0];
          Double_t sigma = par[1];
          Double_t Epeak = par[2];
          Double_t norm = par[3];
          Double_t pigreco = 3.14159265;
          Double_t Aterm, s0, logterm, expterm, logn;
          Double_t logterms0 = eta * 2.35 / 2 + sqrt(1 + pow((eta * 2.35 / 2), 2));
          s0 = (2 / 2.35) * log(logterms0);
          Aterm = eta / (sqrt(2 * pigreco) * sigma * s0);
          logterm = 1 - (eta / sigma) * (x[0] - Epeak);
          if (logterm < 0) {
              logterm = 0.0001; // Evita valori logaritmici negativi
          }
          expterm = log(logterm) / s0;
          expterm = -0.5 * pow(expterm, 2);
          logn = norm * Aterm * exp(expterm);
          return logn;
  }

  TF1 *logn = new TF1("logn", f_logn, 0.998, 1.005, 4);
  logn->SetParameters(-0.3, 0.001, 1.001, 0.1);
  he->Fit(logn);
  TF1 *gaus = new TF1("f_gaus", "gaus", -1, 1);
  ht->Fit(gaus);
  ht->SaveAs("time_shift_scan_histos/time_histo_${i}.root");
  he->SaveAs("time_shift_scan_histos/energyratio_histo_${i}.root");
  std::ofstream outFile("$1", std::ios::app);
  outFile << "$i" << "," << gaus->GetParameter(1) << "," << gaus->GetParError(1) << "," << gaus->GetParameter(2) << "," << gaus->GetParError(2) << "," << logn->GetParameter(2) << "," << logn->GetParError(2) << "," << logn->GetParameter(1) << "," << logn->GetParError(1) << endl;
  outFile.close();
EOF
done

root rootlogon.C << EOF
  new TTree("t", "t");
  t->ReadFile("$1");
  t->Draw("shift:trms*1000:trmserr*1000", "", "goff");
  TGraphErrors *gt = new TGraphErrors(t->GetSelectedRows(), t->GetV1(), t->GetV2(), t->GetV3(), t->GetV3());
  t->Draw("shift:erms*1e2:ermserr*1e2", "", "goff");
  TGraphErrors *ge = new TGraphErrors(t->GetSelectedRows(), t->GetV1(), t->GetV2(), t->GetV3(), t->GetV3());

  new TCanvas("ct", "ct");
  gt->Draw();
  gt->GetYaxis()->SetTitle("#sigma_{t} [ps]");
  gt->GetXaxis()->SetTitle("Pulse shape shift [ns]");
  gt->GetYaxis()->SetTitleOffset(1.5);
  gt->SetTitle("Time resolution");
  ct->BuildLegend(0.6, 0.70, 0.8, 0.88);
  //ct->SetLogy();
  ct->SaveAs("time_shift_scan_histos/sigma_t.root");
  ct->Print("time_shift_scan_histos/sigma_t.ps");
  new TCanvas("ce", "ce");
  ge->Draw();
  ge->GetYaxis()->SetTitle("#sigma_{E} [%]");
  ge->GetYaxis()->SetTitleOffset(1.5);
  ge->GetXaxis()->SetTitle("Pulse shape shift [ns]");
  ge->SetTitle("Energy resolution");
  ce->BuildLegend(0.45, 0.70, 0.55, 0.88);
  //ce->SetLogy();
  ce->SaveAs("time_shift_scan_histos/sigma_e.root");
  ce->Print("time_shift_scan_histos/sigma_e.ps");
EOF

cd time_shift_scan_histos
ps2pdf sigma_t.ps
ps2pdf sigma_e.ps

cd -
