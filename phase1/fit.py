import ROOT as rt
import math
import os
from array import array

rt.gInterpreter.AddIncludePath("/usr/include/eigen3")
rt.gROOT.ProcessLine(".L resample.C+")

sampling = 25


def buildFineGraph(
        fcn,
        xmin=0,
        xmax=12,
        step=0.01
):

    xs  = []
    ys  = []
    exs = []
    eys = []

    x = xmin

    while x <= xmax:

        xs.append(x)
        ys.append(fcn.Eval(x))

        exs.append(0.)
        eys.append(1e-6)

        x += step

    gr = rt.TGraphErrors(
        len(xs),
        array('d', xs),
        array('d', ys),
        array('d', exs),
        array('d', eys)
    )

    return gr


def processFile(
        inputFile,
        outdir="coeffs",
        diagnosticPlots=False,
        diagPercentage=0.1
):

    os.makedirs(outdir, exist_ok=True)

    with open(inputFile) as f:

        for iline, line in enumerate(f):


            diag_rand = rt.gRandom.Uniform() < diagPercentage * 1e-2
            isDiag = diagnosticPlots and diag_rand

            toks = line.strip().split()

            if len(toks) < 18:
                continue

            # ------------------------------------------------
            # Geometry
            # ------------------------------------------------

            isEB  = int(toks[0])
            detid = int(toks[4])

            # ------------------------------------------------
            # Samples
            # ------------------------------------------------

            samples = [
                float(x)
                for x in toks[6:]
            ]

            maxv = max(samples)

            ys = [y / maxv for y in samples]

            xs = [
                i * sampling
                for i in range(len(ys))
            ]

            ex = [0.] * len(xs)
            ey = [1e-4] * len(xs)

            gr = rt.TGraphErrors(
                len(xs),
                array('d', xs),
                array('d', ys),
                array('d', ex),
                array('d', ey)
            )

            # ------------------------------------------------
            # Spline interpolation
            # ------------------------------------------------

            spline = rt.TSpline3(
                f"spline_{detid}",
                gr
            )

            if isDiag:

                c1 = rt.TCanvas(
                    f"c_fit_{detid}",
                    "Spline interpolation",
                    800,
                    600
                )

                gr.SetMarkerStyle(20)
                gr.Draw("AP")

                spline.SetLineColor(rt.kRed)
                spline.SetLineWidth(2)
                spline.Draw("same")

                c1.SaveAs(
                    f"{outdir}/fit_spline_{detid}.root"
                )

                c1.SaveAs(
                    f"{outdir}/fit_spline_{detid}.png"
                )

            # ------------------------------------------------
            # Fine sampled graph
            # ------------------------------------------------

            fineGr = buildFineGraph(
                spline,
                xmin=0,
                xmax=12 * sampling,
                step=1
            )

            fineGr.SetName(
                f"fineGraph_{detid}"
            )

            print(fineGr)

            # ------------------------------------------------
            # Output coeff file
            # ------------------------------------------------

            coeffFile = (
                f"{outdir}/"
                f"crystal_{detid}_coeffs"
            )

            # ------------------------------------------------
            # DIRECT C++ CALL
            # ------------------------------------------------

            rt.resample(
                fineGr,
                detid,
                25,
                coeffFile,
                isDiag
            )

            print(
                f"[{iline}] crystal {detid} done"
            )


# ============================================================

if __name__ == "__main__":

    processFile(
        "template_histograms_ECAL_403687.txt",
        outdir="coeffs",
        diagnosticPlots=True,
        diagPercentage=0.1,
    )
