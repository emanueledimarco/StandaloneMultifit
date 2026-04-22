#include <TGraph.h>
#include <TCanvas.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <TFile.h>
#include <iostream>
#include "PieceWiseCubicSpline.h"


void drawSplineFromFile(const char* file="coeffs_global.txt")
{

    PiecewiseCubicSpline spline(file);

    // Done


    const std::vector<CubicSegment> segments = spline.GetSegments();

    // Determine domain
    double tMin = segments.front().x0;
    double tMax = segments.back().x1;

    // Reconstruct curve
    int N = 2000;
    TGraph* gr = new TGraph(N);

    double dx = (tMax - tMin) / (N - 1);

    for(int i = 0; i < N; i++){
        double x = tMin + i * dx;
        double y = spline.Eval(x);
        gr->SetPoint(i, x, y);
    }

    // Plot
    TCanvas* c = new TCanvas("c", "Spline from file", 900, 600);

    gr->SetLineColor(kBlue);
    gr->SetLineWidth(2);
    gr->Draw("AL");

    gr->SetTitle("Piecewise Cubic Spline from File;X;Y");

    c->BuildLegend();
    c->SaveAs("loaded_spline.png");

    std::cout << "Spline drawn successfully." << std::endl;
}
