import argparse
import ROOT
import numpy as np

def makeCut(cuts,verbose=False):
    print("single cuts = ",cuts)
    cuts = ["("+c+")" for c in cuts]
    cut="*".join(cuts)
    print("sto qua, verbose = ",verbose)
    if verbose:
        print("Selection to be applied = ",cut)
    return cut

def plotSingleResolution(tree, selection=[],verbose=False):
    sel=makeCut(selection,verbose)

    canvas = ROOT.TCanvas("c1", "resolutions", 800, 600)

    reso = ROOT.TH1F("reso","resolution",100,0.5,1.5)
    tree.Draw("samplesReco[3]/amplitudeTruth >> reso",sel)
    
    pass
    
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
                        default=["single_resolution"],
                        help="Analysis that should be run (can be given multiple times)")


    parser.add_argument("-o", "--output",
                        type=str,
                        help="Output file name")

    args = parser.parse_args()

    tfile = ROOT.TFile.Open(args.input_file)
    tree = tfile.Get("Samples")

    if "single_resolution" in args.analysis:
        plotSingleResolution(tree,args.cut,args.verbose)

if __name__ == "__main__":
    main()
