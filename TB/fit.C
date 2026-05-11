#include <TFile.h>
#include <TTree.h>
#include <TSpline.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <iostream>
#include <cmath>

void fit()
{
    // -------------------------
    // Load spline
    // -------------------------
    TFile *fspline = TFile::Open("spline.root");
    if (!fspline || fspline->IsZombie()) {
        std::cerr << "Error opening spline.root\n";
        return;
    }

    TSpline3* f = (TSpline3*) fspline->Get("Spline3");
    if (!f) {
        std::cerr << "Error: TSpline3 not found!\n";
        return;
    }

    // -------------------------
    // Find spline peak
    // -------------------------
    const double Ts = 6.25;
    const double t_ref = -100;

    double t_peak_spline = 0;
    double max_spline = -1e9;

    for (double t = 0; t < 200; t += 0.001) {
        double val = f->Eval(t - t_ref);
        if (val > max_spline) {
            max_spline = val;
            t_peak_spline = t;
        }
    }

    // -------------------------
    // Open waveform file
    // -------------------------
    TFile *fin = TFile::Open("19535_0002_unpacked.root");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error opening data file\n";
        return;
    }

    TTree *tree = (TTree*) fin->Get("h4");
    if (!tree) {
        std::cerr << "Tree not found!\n";
        return;
    }

    Short_t xtal_sample[675][40];
    tree->SetBranchAddress("xtal_sample", xtal_sample);

    // -------------------------
    // OUTPUT SETUP
    // -------------------------
    TFile* fout = new TFile("results.root", "RECREATE");
    TTree* out = new TTree("T", "fit results");

    const int channels[9] = {5,6,13,14,15,16,109,110,119};
    const int NCH = 9;

    int event;

    int out_channels[NCH];
    double fit_time[NCH];
    double fit_amp[NCH];
    double charge[NCH];

    out->Branch("event", &event, "event/I");
    out->Branch("channel_ids", out_channels, "channel_ids[9]/I");
    out->Branch("fit_time", fit_time, "fit_time[9]/D");
    out->Branch("fit_amp",  fit_amp,  "fit_amp[9]/D");
    out->Branch("charge", charge, "charge[9]/D");

    int nentries = tree->GetEntries();

    // =========================
    // EVENT LOOP
    // =========================
    for (event = 0; event < nentries; event++)
    {
        tree->GetEntry(event);
        cout << event << endl;
        const int N = 40;

        // copy channel IDs every event (IMPORTANT)
        for (int i = 0; i < NCH; i++)
            out_channels[i] = channels[i];

        // =========================
        // CHANNEL LOOP
        // =========================
        for (int ich = 0; ich < NCH; ich++)
        {
            int channel = channels[ich];

            double q = 0.0, pedestal = 0.0;

            // -------------------------
            // Peak finding
            // -------------------------
            int jmax_wave = 10;
            double max_wave = -1e9;

            for (int j = 10; j < N; j++) {
                double val = xtal_sample[channel][j] - pedestal;
                if (val > max_wave) {
                    max_wave = val;
                    jmax_wave = j;
                }
            }

            double t_peak_wave = (jmax_wave - 10) * Ts;

            // sub-sample correction
            if (jmax_wave > 10 && jmax_wave < N-1) {
                double y1 = xtal_sample[channel][jmax_wave-1] - pedestal;
                double y2 = xtal_sample[channel][jmax_wave]   - pedestal;
                double y3 = xtal_sample[channel][jmax_wave+1] - pedestal;

                double denom = (y1 - 2*y2 + y3);
                if (fabs(denom) > 1e-6) {
                    double delta = 0.5 * (y1 - y3) / denom;
                    t_peak_wave += delta * Ts;
                }
            }

            double delta_t_align = t_peak_wave - t_peak_spline;

            // -------------------------
            // FIT WINDOW
            // -------------------------
            int j_start = jmax_wave - 3;
            int j_end   = jmax_wave + 7;

            if (j_start < 0) j_start = 0;
            if (j_end >= N) j_end = N - 1;

            int nped = 0;

            int j_ped_start = j_start - 10;
            int j_ped_end   = j_start - 1;

            if (j_ped_start < 0) j_ped_start = 0;

            for (int j = j_ped_start; j <= j_ped_end; j++)
            {
                pedestal += xtal_sample[channel][j];
                nped++;
            }

            if (nped > 0)
                pedestal /= nped;

            for (int j = j_start; j <= j_end; j++)
            {
                q += (xtal_sample[channel][j] - pedestal);
            }

            int Nfit = j_end - j_start + 1;

            double delta_t_align_iter = delta_t_align;

            double A = 0;
            double C = 0;

            // -------------------------
            // ITERATIONS
            // -------------------------
            for (int iter = 0; iter < 3; iter++)
            {
                TMatrixD M(Nfit, 2);
                TVectorD w(Nfit);

                for (int k = 0, j = j_start; j <= j_end; j++, k++)
                {
                    double t = (j - 10) * Ts;
                    double t_shifted = t - delta_t_align_iter;

                    double P  = f->Eval(t_shifted - t_ref);
                    double dP = f->Derivative(t_shifted - t_ref);

                    if (t_shifted - t_ref < 65) { P = 0; dP = 0; }

                    M(k,0) = P;
                    M(k,1) = dP;
                    w(k)   = xtal_sample[channel][j] - pedestal;
                }

                TDecompSVD svd(M);
                Bool_t ok;
                TVectorD x = svd.Solve(w, ok);

                if (!ok) {
                    A = 0;
                    C = 0;
                    break;
                }

                A = x[0];
                C = x[1];

                if (A <= 0) A = 0;

                double dt = -C / A;
                delta_t_align_iter += dt;
            }

            // -------------------------
            // FINAL RESULTS
            // -------------------------
            fit_amp[ich]  = A;
            fit_time[ich] = t_ref + delta_t_align_iter;
            charge[ich] = q;
        }

        out->Fill();
    }

    fout->Write();
    fout->Close();

    std::cout << "Done. Results saved in results.root\n";
}
