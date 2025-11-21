import argparse
import ROOT
import numpy as np
import array

def setStyle():
    ROOT.gStyle.SetTitleFont(132,"xyz");
    ROOT.gStyle.SetTitleFont(132," ");
    ROOT.gStyle.SetTitleSize(0.06,"xyz");
    ROOT.gStyle.SetTitleSize(0.06," ");
    ROOT.gStyle.SetLabelFont(132,"xyz");
    ROOT.gStyle.SetLabelSize(0.05,"xyz");
    ROOT.gStyle.SetTextFont(132);
    ROOT.gStyle.SetTextSize(0.08);
    ROOT.gStyle.SetStatFont(132);

def makeCut(cuts,verbose=False):
    cuts = ["("+c+")" for c in cuts]
    cut="*".join(cuts)
    if verbose:
        print("Selection to be applied = ",cut)
    return cut

def cbFit(h,name,title,xmin=-1,xmax=-1):
    # --- Create RooFit variables ---
    x = ROOT.RooRealVar("x", "E/E_{True}", h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())

    # --- Define the fit subrange ---
    if xmin>=0 and xmax>=0:
        x.setRange("fitRange", xmin, xmax)

    # --- Import histogram into RooDataHist ---
    data = ROOT.RooDataHist("data", "data", ROOT.RooArgList(x), h)

    # ---------------------------
    # Crystal Ball parameters
    # ---------------------------
    meanCB  = ROOT.RooRealVar("meanCB",  "CB mean",  h.GetMean(), h.GetMean()-1, h.GetMean()+1)
    sigmaCB = ROOT.RooRealVar("sigmaCB", "CB sigma", h.GetRMS(), 0.1*h.GetRMS(), 5*h.GetRMS())
    alpha   = ROOT.RooRealVar("alpha",   "alpha",    1.5, 0.1, 5.0)
    n       = ROOT.RooRealVar("n",       "n",        3.0, 0.5, 20.0)
    
    cb = ROOT.RooCBShape("cb", "Crystal Ball", x, meanCB, sigmaCB, alpha, n)

    # ---------------------------
    # Gaussian parameters
    # ---------------------------
    meanG  = ROOT.RooRealVar("meanG",  "Gauss mean",  h.GetMean(), h.GetMean()-1, h.GetMean()+1)
    sigmaG = ROOT.RooRealVar("sigmaG", "Gauss sigma", h.GetRMS()/2, 0.1*h.GetRMS(), 5*h.GetRMS())
    
    gauss = ROOT.RooGaussian("gauss", "Gaussian", x, meanG, sigmaG)
    
    # ---------------------------
    # Fraction and total yield
    # ---------------------------
    fG = ROOT.RooRealVar("fG", "Gaussian fraction", 0.5, 0.0, 1.0)   # fraction for Gaussian
    Ntot = ROOT.RooRealVar("Ntot", "total yield", h.Integral(), 0, 10*h.Integral())

    # ---------------------------
    # Combined PDF: fG * Gauss + (1 - fG) * CB
    # ---------------------------
    shape = ROOT.RooAddPdf("shape", "Gauss + CB (fractions)",
                           ROOT.RooArgList(gauss, cb),
                           ROOT.RooArgList(fG))
    
    # Extended model with a single normalization
    model = ROOT.RooAddPdf("model", "extended model",
                           ROOT.RooArgList(shape),
                           ROOT.RooArgList(Ntot))

    # ---------------------------
    # Fit (extended likelihood)
    # ---------------------------
    result = model.fitTo(data,
                          ROOT.RooFit.Range("fitRange"),
                          ROOT.RooFit.Extended(True),
                          ROOT.RooFit.Save(),
                          ROOT.RooFit.PrintLevel(-1))


    # ---------------------------
    # Plot
    # ---------------------------
    frame = x.frame(ROOT.RooFit.Title(f"{title}"))
    data.plotOn(frame)
    model.plotOn(frame)
    
    # Individual components
    model.plotOn(frame, ROOT.RooFit.Components("gauss"),
                 ROOT.RooFit.LineColor(ROOT.kBlue),
                 ROOT.RooFit.LineStyle(ROOT.kDashed))
    
    model.plotOn(frame, ROOT.RooFit.Components("cb"),
                 ROOT.RooFit.LineColor(ROOT.kRed),
                 ROOT.RooFit.LineStyle(ROOT.kDashed))
    
    c = ROOT.TCanvas("c", "Extended Fit", 800, 600)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    frame.GetYaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetTitleOffset(1.1)

    frame.Draw()
    
    # Add info box with mean and sigma
    pt = ROOT.TPaveText(0.60, 0.65, 0.88, 0.88, "NDC")
    pt.SetFillColor(0)
    pt.SetTextFont(42)
    pt.SetBorderSize(0)
    pt.SetTextSize(0.05)
    
    pt.AddText(f"m_{{core}} = {meanG.getVal():.3f} #pm {meanG.getError():.3f}")
    pt.AddText(f"#sigma_{{core}} = {sigmaG.getVal():.3f} #pm {sigmaG.getError():.3f}")
    
    pt.Draw()
    
    for ext in ["png","pdf"]:
        c.SaveAs(f"{name}.{ext}")

    print("Fit results:")
    result.Print()
    return {"mean":(meanG.getVal(),meanG.getError()),
            "sigma":(sigmaG.getVal(),sigmaG.getError())}
    
def plotSingleResolution(tree,name,title,selection=[],verbose=False):
    sel=makeCut(selection,verbose)

    canvas = ROOT.TCanvas("c1", "resolutions", 800, 600)

    tree.Draw("samplesReco[3]/amplitudeTruth >> resotemp",sel,"goff")
    resotemp = ROOT.gDirectory.Get("resotemp")
    m = resotemp.GetMean()
    s = resotemp.GetRMS()

    reso = ROOT.TH1F("reso","resolution",60,m-5*s,m+5*s)
    tree.Draw("samplesReco[3]/amplitudeTruth >> reso",sel)
    results = cbFit(reso,name,title)
    return results
    
def plotDifferentialResolution(tree,selection=[],verbose=False):
    Ebins=[1,2,5,10,15,20,25,30,50,70,100]
    x,ex,b,eb,s,es = [],[],[],[],[],[]
    print ("Energy bins to be analysed: ",Ebins)
    for ie in range(len(Ebins)-1):
        print ("ie = ",ie)
        addcut = f"signalTruth>{Ebins[ie]} && signalTruth<={Ebins[ie+1]}"
        print ("Processing bin: ",addcut)
        fullsel = selection + [addcut]

        name = f"resolution_E{Ebins[ie]}To{Ebins[ie+1]}"
        title = f"{Ebins[ie]} GeV < E < {Ebins[ie+1]} GeV"
        results = plotSingleResolution(tree,name,title,fullsel)
        b.append(results["mean"][0])
        eb.append(results["mean"][1])
        s.append(results["sigma"][0])
        es.append(results["sigma"][1])
        x.append(np.mean([Ebins[ie],Ebins[ie+1]]))
        ex.append(0)
        

    # Convert Python lists to C-style arrays
    x_arr = array.array('d', x)
    ex_arr = array.array('d', ex)
    b_arr = array.array('d', b)
    eb_arr = array.array('d', eb)
    s_arr = array.array('d', s)
    es_arr = array.array('d', es)
    s2_arr = array.array('d', [si/bi for si,bi in zip(s,b)])
    
    gbias = ROOT.TGraphErrors(len(x), x_arr, b_arr, ex_arr, eb_arr)
    gsigma = ROOT.TGraphErrors(len(x), x_arr, s2_arr, ex_arr, es_arr)

    # Optional: style
    gbias.SetTitle("")
    gbias.SetMarkerStyle(20)
    gbias.GetXaxis().SetTitle("Amplitude [GeV]")
    gbias.GetYaxis().SetTitle("Mean of E/E_{True}")
    
    gsigma.SetTitle("")
    gsigma.SetMarkerStyle(20)
    gsigma.GetXaxis().SetTitle("Amplitude [GeV]")
    gsigma.GetYaxis().SetTitle("#sigma of E/E_{True}")

    # Draw
    c = ROOT.TCanvas("c", "TGraphErrors from arrays", 800, 600)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    gbias.GetYaxis().SetTitleOffset(1.1)
    gbias.GetXaxis().SetTitleOffset(1.1)
    gsigma.GetYaxis().SetTitleOffset(1.1)
    gsigma.GetXaxis().SetTitleOffset(1.1)
    gsigma.GetYaxis().SetRangeUser(0,0.20)

    gbias.Draw("APC")
    for ext in ["png","pdf","root"]:
        c.SaveAs(f"bias_vs_E.{ext}")
    gsigma.Draw("APC")
    for ext in ["png","pdf","root"]:
        c.SaveAs(f"sigma_vs_E.{ext}")
    
def main():
    parser = argparse.ArgumentParser(
        description="Script to run simple analysis from a reconstructed multifit TTree"
    )

    # Required positional argument
    parser.add_argument("input_file", help="Path to the input file")

    # Optional flags
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Enable verbose mode")

    # Option that can appear multiple times
    parser.add_argument("-c", "--cut",
                        action="append",
                        default=["1"],
                        help="cut to be added (can be given multiple times)")

    parser.add_argument("-a", "--analysis",
                        action="append",
                        default=[],
                        help="Analysis that should be run (can be given multiple times)")


    parser.add_argument("-o", "--output",
                        type=str,
                        help="Output file name")

    args = parser.parse_args()

    tfile = ROOT.TFile.Open(args.input_file)
    tree = tfile.Get("RecoAndSim")

    setStyle()

    print (" ==== Analyses to be run: ====\n", args.analysis)
    
    if "single_resolution" in args.analysis:
        plotSingleResolution(tree,"resolution"," AND ".join(args.cut),args.cut,args.verbose)

    if "differential_resolution" in args.analysis:
        plotDifferentialResolution(tree,args.cut,args.verbose)

if __name__ == "__main__":
    main()
