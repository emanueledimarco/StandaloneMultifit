TFile *_file0 = TFile::Open("sigma_e.root")
ce->Draw()
//takes sigma_e
int n_e = sigma_e->GetN();
double* y_e = sigma_e->GetY();
for (int i = 0; i < n_e; ++i) {
    y_e[i] = sqrt(y_e[i]*y_e[i] + 0.29*0.29);
}
double ymin_e = y_e[0];
for (int i = 1; i < n_e; ++i) if (y_e[i] < ymin_e) ymin_e = y_e[i];
for (int i = 0; i < n_e; ++i) y_e[i] /= ymin_e;
sigma_e->SaveAs("sigma_e_norm_conv_photostat_term_0.29percent_graph.root")

TFile *_file0 = TFile::Open("sigma_t.root")
ct->Draw()
//takes sigma_t
int n_t = sigma_t->GetN();
double* y_t = sigma_t->GetY();
for (int i = 0; i < n_t; ++i) {
    y_t[i] = sqrt(y_t[i]*y_t[i] + 0.017*0.017);
}
double targetX[2] = {-1.0, 0.0};
double epsilon = 1e-2;
for (int t = 0; t < 2; ++t) {
    int n = sigma_t->GetN();
    double x_val, y_val;
    int removeIndex = -1;
    for (int i = 0; i < n; ++i) {
        sigma_t->GetPoint(i, x_val, y_val);
        if (fabs(x_val - targetX[t]) < epsilon) {
            removeIndex = i;
            break;
        }
    }
    if (removeIndex >= 0) {
        for (int i = removeIndex; i < n-1; ++i) {
            double x_next, y_next;
            sigma_t->GetPoint(i+1, x_next, y_next);
            sigma_t->SetPoint(i, x_next, y_next);
        }
        sigma_t->Set(n-1);
    }
}
n_t = sigma_t->GetN();   // number of points may have changed after removal
double* y_t2 = sigma_t->GetY();
double ymin_t = y_t2[0];
for (int i = 1; i < n_t; ++i) if (y_t2[i] < ymin_t) ymin_t = y_t2[i];
for (int i = 0; i < n_t; ++i) y_t2[i] /= ymin_t;
sigma_t->SaveAs("sigma_t_norm_conv_const_term_17ps_graph.root")
