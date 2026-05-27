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
        logterm = 1 - (eta / sigma) * (Epeak - x[0]);
        if (logterm < 0) {
            logterm = 0.0001; // Evita valori logaritmici negativi
        }
        expterm = log(logterm) / s0;
        expterm = -0.5 * pow(expterm, 2);
        logn = norm * Aterm * exp(expterm);
        return logn;
}

void logn_fit(TH1* h){
  double mean = h->GetMean();
  double sigma = h->GetRMS();
  TF1 *logn = new TF1("logn", f_logn, mean-4*sigma, mean+3*sigma, 4);
  logn->SetParameters(1, sigma, mean, 0.3*h->Integral());
  logn->SetParNames("#eta", "#sigma", "Peak", "Norm.");
  logn->SetParLimits(0, 0.01, 1);
  logn->SetParLimits(1, 0.0001, 10);
  logn->SetParLimits(3, 0, 10000);
  h->Fit(logn, "R");
}
