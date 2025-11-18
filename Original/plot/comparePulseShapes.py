import ROOT
import numpy as np

def style_graph(g, color, marker):
    g.SetLineColor(color)
    g.SetMarkerColor(color)
    g.SetMarkerStyle(marker)
    g.SetLineWidth(2)
    g.SetMarkerSize(1.2)

def comparePulseShapes(pulse_graphs, samplings, offsets, titles, xmin, xmax):

    canvas = ROOT.TCanvas("c1", "Pulse Shapes", 800, 600)

    leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)

    gr_resampled = []
    
    for i,gr in enumerate(pulse_graphs):
        x,y = [],[]
        nbins = int((xmax-xmin)/samplings[i])
        x = np.linspace(xmin,xmax,nbins+1)
        offs = np.full(len(x),offsets[i])
        y = np.array([gr.Eval(xval) for xval in x+offs])
        print ("x = ",x,"\ny = ",y)
        print (f"{gr.GetName()}_resampled_{i}")
        gr_resampled.append(ROOT.TGraph(len(x),x,y))
        gr_resampled[i].SetName(f"{gr.GetName()}_resampled_{i}")
        gr_resampled[i].SetTitle("Pulse Shape")
        gr_resampled[i].SetMarkerStyle(8)
        gr_resampled[i].SetMarkerColor(2+i)
        gr_resampled[i].SetLineColor(2+i)
        
        draw_opt = "APC" if i == 0 else "PC"
        gr_resampled[i].Draw(draw_opt)
        gr_resampled[i].GetXaxis().SetRangeUser(xmin-max(offsets),xmax)
        gr_resampled[i].GetYaxis().SetRangeUser(0,1.2)
        gr_resampled[i].GetXaxis().SetTitle("time (ns)")
        
        leg.AddEntry(gr_resampled[i], titles[i], "lp")
        print ("done")
        
    leg.Draw()
    canvas.SaveAs("pulse_shapes.pdf")
        

if __name__ == "__main__":
    
    files = ["data/EmptyFileCRRC43.root",
             "data/EmptyFileTestBeamPhase2.root"]

    titles = ["Phase-1", "Phase-2"]

    offsets = [12.5,0.0]
    
    graphs = []
    for f in files:
        tf = ROOT.TFile.Open(f)
        gr = tf.Get("PulseShape/grPulseShape")
        graphs.append(gr)

    comparePulseShapes(graphs,[25,6.25],offsets,titles,0,400)
