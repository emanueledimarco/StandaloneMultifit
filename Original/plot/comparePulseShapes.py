import ROOT
import numpy as np

def style_graph(g, color, marker):
    g.SetLineColor(color)
    g.SetMarkerColor(color)
    g.SetMarkerStyle(marker)
    g.SetLineWidth(2)
    g.SetMarkerSize(1.2)


def comparePulseShapes(pulse_graphs, samplings, offsets,
                       titles, xmin, xmax,
                       pedestal_time=0):

    canvas = ROOT.TCanvas("c1", "Pulse Shapes", 800, 600)

    leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)

    gr_resampled = []
    labels = []
    vlines = []
    colors = [ROOT.kRed,ROOT.kBlue]
    
    # Estensione a tempi negativi
    xmin_extended = xmin - pedestal_time

    for i, gr in enumerate(pulse_graphs):

        nbins = int((xmax - xmin_extended) / samplings[i])

        print(f"i = {i}, {samplings[i]}, {xmin_extended}, {xmax}, {(xmax - xmin_extended)}, {nbins}")

        # Campionamento esteso
        x_original = np.linspace(xmin_extended, xmax, nbins + 1)

        # shift asse x per partire da 0
        x = x_original - xmin_extended
        #x = np.linspace(xmin_extended, xmax, nbins + 1)
        
        offs = np.full(len(x), offsets[i])

        # Piedistallo a zero per t<0
        y = np.array([
            0 if xval < 0 else gr.Eval(xval + offsets[i])
            for xval in x_original
        ])

        print("x =", x, "\ny =", y)

        graph = ROOT.TGraph(len(x), x.astype(np.float64),
                            y.astype(np.float64))

        graph.SetName(f"{gr.GetName()}_resampled_{i}")
        graph.SetTitle("Pulse Shape")
        graph.SetMarkerStyle(8)
        graph.SetMarkerColor(colors[i])
        graph.SetLineColor(colors[i])

        gr_resampled.append(graph)

        draw_opt = "APC" if i == 0 else "PC"
        graph.Draw(draw_opt)

        graph.GetXaxis().SetRangeUser(0, xmax - xmin_extended)
        graph.GetYaxis().SetRangeUser(0, 1.2)
        graph.GetXaxis().SetTitle("time (ns)")
        graph.GetYaxis().SetTitle("Amplitude")

        leg.AddEntry(graph, titles[i], "lp")

        # Scrivi numeri e linee verticali sui punti rossi
        if i == 0:
            for j, (xj, yj) in enumerate(zip(x, y)):

                # linea verticale dal punto a y=0
                line = ROOT.TLine(xj, 0, xj, yj)
                line.SetLineColor(ROOT.kRed)
                line.SetLineStyle(2)   # tratteggiata (opzionale)
                line.Draw()
         
                vlines.append(line)  # evita garbage collection
         
                # numero sopra il punto
                text = ROOT.TLatex()
                text.SetTextSize(0.025)
                text.SetTextAlign(22)
         
                text.DrawLatex(xj, yj + 0.03, str(j))
                labels.append(text)
    print("done")

    leg.Draw()
    canvas.Update()
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

    comparePulseShapes(graphs,[25,6.25],offsets,titles,-75,175)
