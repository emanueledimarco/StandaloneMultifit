import ROOT as rt
import math
import os
from array import array
from alphaBetaFitter import AlphaBetaFitter

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
        array('d',xs),
        array('d',ys),
        array('d',exs),
        array('d',eys)
    )

    return gr



def processFile(
        inputFile,
        outdir="coeffs",
        diagnosticPlots=False,
        nDiag=10
):

    doEB = True #sbagliato, deve essere per canale

    os.makedirs(outdir,exist_ok=True)

    fitter = AlphaBetaFitter(doEB)

    with open(inputFile) as f:

        for iline, line in enumerate(f):

            isDiag = diagnosticPlots and (iline < nDiag)

            if diagnosticPlots and iline > nDiag: break

            toks = line.strip().split()

            if len(toks) < 18:
                continue

            # ------------------------------------------------
            # Geometry
            # ------------------------------------------------

            ieta  = int(toks[0])
            iphi  = int(toks[1])
            ic    = int(toks[2])
            iz    = int(toks[3])

            detid = int(toks[4])

            # ------------------------------------------------
            # Samples
            # ------------------------------------------------

            samples = [
                float(x)
                for x in toks[6:]
            ]

            maxv = max(samples)

            ys = [y/maxv for y in samples]

            xs = [
                i * sampling
                for i in range(len(ys))
            ]

            ex = [0.] * len(xs)
            ey = [1e-4] * len(xs)

            gr = rt.TGraphErrors(
                len(xs),
                array('d',xs),
                array('d',ys),
                array('d',ex),
                array('d',ey)
            )

            # ------------------------------------------------
            # Fit
            # ------------------------------------------------

            pars, errs = fitter.fit(gr)

            fcn = fitter.getFcn()

            fitFcn = fcn.Clone(
                f"alphabeta_{detid}"
            )

            if isDiag:
                c1 = rt.TCanvas(
                    f"c_fit_{detid}",
                    "AlphaBeta fit",
                    800,
                    600
                )

                gr.SetMarkerStyle(20)
                gr.Draw("AP")

                fcn.Draw("same")

                c1.SaveAs(
                    f"{outdir}/fit_alphabeta_{detid}.root"
                )

            # ------------------------------------------------
            # Fine sampled graph
            # ------------------------------------------------

            fineGr = buildFineGraph(
                fitFcn,
                xmin=0,
                xmax=12*sampling,
                step=1
            )

            fineGr.SetName(
                f"fineGraph_{detid}"
            )

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
        nDiag=10
    )
