#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TSpline.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <iostream>
#include <vector>
#include <cmath>

void prova_fixed(int ntoy = 10000)
{
    TSpline3* f = (TSpline3*) gDirectory->Get("Spline3");
    if (!f) {
        std::cerr << "Error: TSpline3 'f' not found!" << std::endl;
        return;
    }

    const double tmin = 100.0, tmax = 180.0, Ts = 6.25, noise_sigma = 0.2, t_ref = 0.0;

    TFile* fout = new TFile("results.root", "RECREATE");
    TTree* tree = new TTree("T", "toy fit results");

    double true_time, true_amp, fit_time, fit_amp;
    tree->Branch("true_time", &true_time, "true_time/D");
    tree->Branch("true_amp",  &true_amp,  "true_amp/D");
    tree->Branch("fit_time",  &fit_time,  "fit_time/D");
    tree->Branch("fit_amp",   &fit_amp,   "fit_amp/D");

    TRandom3 rng(0);

    for (int itoy = 0; itoy < ntoy; ++itoy) {
        // --- generate true parameters ---
        true_time = rng.Uniform(-5.0, 5.0);   // ns
        true_amp  = rng.Uniform(1.0, 100.0); // amplitude

        // --- generate samples ---
        std::vector<double> t_samp, w_samp;
        for (double t = tmin; t <= tmax; t += Ts) {
            t_samp.push_back(t);
            double val = true_amp * f->Eval(t - true_time);
            if (noise_sigma > 0) val += rng.Gaus(0, noise_sigma);
            w_samp.push_back(val);
        }

        const int N = t_samp.size();

        // --- first linear fit: w_j = A * P_j + C * P'_j ---
        TMatrixD M(N, 2);
        TVectorD w(N);
        for (int j = 0; j < N; ++j) {
            double P  = f->Eval(t_samp[j] - t_ref);
            double dP = f->Derivative(t_samp[j] - t_ref);
            M(j, 0) = P;
            M(j, 1) = dP;
            w[j] = w_samp[j];
        }

        // Optional: check correlation
        double meanP = 0, meanPd = 0;
        for (int j = 0; j < N; ++j) { meanP += M(j,0); meanPd += M(j,1); }
        meanP /= N; meanPd /= N;
        double num=0, denomP=0, denomPd=0;
        for (int j = 0; j < N; ++j) {
            double dP  = M(j,0) - meanP;
            double dPd = M(j,1) - meanPd;
            num     += dP * dPd;
            denomP  += dP * dP;
            denomPd += dPd * dPd;
        }
        double corr = num / std::sqrt(denomP * denomPd);
        //std::cout << "[Toy " << itoy << "] correlation(P,P') = " << corr << std::endl;

        // Solve first LS using SVD
        TDecompSVD svd(M);
        Bool_t ok;
        TVectorD x = svd.Solve(w, ok);
        if (!ok) { fit_amp = fit_time = -999; tree->Fill(); continue; }

        double A = x[0];
        double C = x[1];

        if (A <= 0) A = 0;

        // --- compute delta t ---
        double DeltaT = - C / A;

        // --- second fit: shift the spline and refit amplitude only ---
        double sumP2 = 0, sumPw = 0;
        for (int j = 0; j < N; ++j) {
            double P_shifted = f->Eval(t_samp[j] - t_ref - DeltaT);
            sumP2  += P_shifted * P_shifted;
            sumPw  += P_shifted * w_samp[j];
        }
        fit_amp  = sumPw / sumP2;
        fit_time = DeltaT;

        tree->Fill();
    }

    fout->Write();
    fout->Close();

    std::cout << "All toys done. Output written to results.root" << std::endl;
}
