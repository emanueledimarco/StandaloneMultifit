import ROOT
import math,random
import numpy as np

# dummy phase2 PS params (time units in 1/4 ns)
#                     0 1  2 3 4 5 
ps_params = np.array([1,18,3,5,3,1], dtype=float)

def pyf_tf1_coulomb(x, p):
    return p[1] * x[0] * x[1] / (p[0]**2) * math.exp(-p[2] / p[0])

def pyf_f(t, t0, tau_r):
    tshift = t-t0
    return 1./(np.exp(-1./tau_r * tshift) + 1)

# fast decay time (fluorescence)
def pyf_g(t, t0, tau_f):
    tshift = t-t0
    return 1./(np.exp(tshift/tau_f) + 1)

# slow decay time (phosphorence)
def pyf_h(t, t0, tau_s):
    tshift = t-t0
    return 1./(np.exp(tshift/tau_s) + 1)

def pyf_total(x, p):
    # p[0] = A = pulse amplitude
    # p[1] = t0 = pulse arrival time
    # p[2] = tau_r = rise time (~1 ns)
    # p[3] = tau_f = fast decay time (~5 ns) - fluorescence
    # p[4] = tau_s = fast decay time (~50 ns) - phosphorence
    # p[5] = R = relative weight of the decay times (R=1 means slow decay component is 0. In any case 0 <= R <=1)
    t = x
    return p[0] * pyf_f(t,p[1],p[2]) * (p[5] * pyf_g(t,p[1],p[3]) + (1-p[5]) * pyf_h(t,p[1],p[4])) 

def copy_dummyPUpdfs(sourceFile):
    dirname = "PileupPDFs"
    savdir = ROOT.gDirectory
    targetdir = savdir.mkdir(dirname)

    f = ROOT.TFile.Open(sourceFile)
    f.cd(dirname)
    sourcedir = ROOT.gDirectory
    
    keys = sourcedir.GetListOfKeys()
    for k in keys:
        obj = k.ReadObj()
        print("\t histo ",obj.GetName())
        targetdir.cd()
        obj.Write()
        del obj
    savdir.SaveSelf(ROOT.kTRUE)
    
def makeInputPSGraphFromIdealFunction():

    print("==> create the dummy pulse shape from an ideal parametric function")
    n = 4000 # steps in 1/4 ns
    x = np.linspace(0,1000,n+1)
    y = pyf_total(x,ps_params)

    # for pulse shape need normalized amplitude (1 at maximum)
    y = y/np.max(y)
    
    outFile = ROOT.TFile("data/EmptyFileIdealPSphase2.root","recreate")
    outFile.mkdir("PulseShape","PulseShape")
    outFile.cd("PulseShape")
    
    gr = ROOT.TGraph(n,x,y)
    gr.SetName("grPulseShape")
    gr.SetMarkerStyle(ROOT.kFullCircle)
    gr.SetMarkerSize(0.5)
    # gr.GetXaxis().SetRangeUser(0,50)
    
    # c = ROOT.TCanvas("c","c",600,600)
    # gr.Draw("ALP")
    # c.SaveAs("pstest.pdf")
    gr.Write()
    outFile.cd("../")
    print("pulse shape written.")
    print("==> Now copying the input pileup PDFs...")
    
    # now copy Pileup PDFs from the original file
    copy_dummyPUpdfs("data/EmptyFileCRRC43.root")
    
    outFile.Close()
    print("Done.")

def graph_to_arrays(graph):
    """
    Convert a ROOT TGraphErrors into numpy arrays.
    
    Parameters
    ----------
    graph : ROOT.TGraphErrors
        The input graph.

    Returns
    -------
    x : np.ndarray
        Array of x values.
    y : np.ndarray
        Array of y values.
    ex : np.ndarray
        Array of x errors.
    ey : np.ndarray
        Array of y errors.
    """
    n = graph.GetN()
    
    # Create numpy arrays
    x = np.array([graph.GetPointX(i) for i in range(n)], dtype=float)
    y = np.array([graph.GetPointY(i) for i in range(n)], dtype=float)
    ex = np.array([graph.GetErrorX(i) for i in range(n)], dtype=float)
    ey = np.array([graph.GetErrorY(i) for i in range(n)], dtype=float)
    
    return x, y, ex, ey
    
def arrays_to_graph(x, y, ex=None, ey=None, name="graph"):
    """
    Create a ROOT.TGraphErrors from numpy arrays.
    
    Parameters
    ----------
    x, y : array-like
        Arrays of x and y values.
    ex, ey : array-like or None
        Arrays of x and y errors (optional). If None, errors are set to zero.
    name : str
        Name of the resulting graph.
    
    Returns
    -------
    ROOT.TGraphErrors
        A TGraphErrors object containing the input data.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    
    if ex is None:
        ex = np.zeros_like(x)
    if ey is None:
        ey = np.zeros_like(y)
    
    n = len(x)
    graph = ROOT.TGraphErrors(n, x, y, ex, ey)
    graph.SetName(name)
    
    return graph



def makeInputPSGraphFromTestBeam(tbfile="pulse_tgraph_fromTB.root",graph="pulse_shape"):

    print("taking Pulse shape ",graph," from test-beam file ",tbfile)
    f = ROOT.TFile.Open(tbfile)
    gr = f.Get(graph)

    x, y, ex, ey = graph_to_arrays(gr)
    
    # shift the graph to the left
    x_sel = x[y>1e-3]
    y_sel = y[y>1e-3]
    ex_sel = ex[y>1e-3]
    ey_sel = ey[y>1e-3]
    
    xmin=min(x_sel)
    print ("Shifting xmin by ",xmin)
    x_sel -= xmin
    
    outFile = ROOT.TFile("data/EmptyFileTestBeamPhase2.root","recreate")
    outFile.mkdir("PulseShape","PulseShape")
    outFile.cd("PulseShape")

    grPs = arrays_to_graph(x_sel,y_sel,ex_sel,ey_sel,"PulseShape")
    
    grPs.SetName("grPulseShape")
    grPs.SetMarkerStyle(ROOT.kFullCircle)
    grPs.SetMarkerSize(0.5)
    grPs.Write()
    outFile.cd("../")
    print("pulse shape written.")
    print("==> Now copying the input pileup PDFs...")
    
    # now copy Pileup PDFs from the original file
    copy_dummyPUpdfs("data/EmptyFileCRRC43.root")
    
    outFile.Close()
    print("Done.")


def two_closest_points(value, data):
    data = np.asarray(data, dtype=float)
    # Compute absolute differences
    diff = np.abs(data - value)
    # Get indices of the two smallest differences
    indices = np.argsort(diff)[:2]
    # Return both values and indices
    return data[sorted(indices)], sorted(indices)

def correlation_from_covariance(covariance):
    v = np.sqrt(np.diag(covariance))
    outer_v = np.outer(v, v)
    correlation = covariance / outer_v
    correlation[covariance == 0] = 0
    return correlation


def numpy_to_TH2F(matrix, name="PulseCovariance", title="Template Covariance Matrix"):
    if matrix.ndim != 2:
        raise ValueError("Input must be a 2D array.")

    ny, nx = matrix.shape  # rows (y), cols (x)
    # Axis ranges: 0..nx and 0..ny
    h2 = ROOT.TH2F(name, title, nx, 0, nx, ny, 0, ny)

    for iy in range(ny):
        for ix in range(nx):
            h2.SetBinContent(ix + 1, iy + 1, float(matrix[iy, ix]))

    return h2

def makePulseCovariance(rfile,pulse_shape,time_bias=0,time_spread=100*1e-3,n_events=1e4):

    NFREQ = 6.25
    NSAMPLES = 16;
    NPRESAMPLES = 6
    maxSample = 9

    fout = ROOT.TFile.Open("data/PulseCovarianceTestBeamPhase2.root","recreate")
    
    f = ROOT.TFile.Open(rfile)
    ps = f.Get(pulse_shape)

    x, y, ex, ey = graph_to_arrays(ps)
    ymax = max(y)
    
    covmat = np.matrix(np.zeros((16,16), dtype=np.float64))
    t0s = []

    template_atMax = ps.Eval((maxSample-NPRESAMPLES) * NFREQ)
    print ("template val at max sample = ",template_atMax)
    
    run_ps_params = ps_params.copy()
    for i in range(n_events):
        tshift = random.uniform(-1*time_spread,time_spread)
        t0s.append(18.+tshift)

        run_ps_params[1] = 18. + tshift
        pulse_atMax = ps.Eval((maxSample-NPRESAMPLES) * NFREQ + tshift)
        
        # now take the digitized points
        for i in range(NSAMPLES):
            for j in range(NSAMPLES):
                sample_i = ps.Eval(i * NFREQ + tshift)/pulse_atMax;
                template_i = ps.Eval(i * NFREQ)/template_atMax;
                sample_j = ps.Eval(j * NFREQ + tshift)/pulse_atMax;
                template_j = ps.Eval(j * NFREQ)/template_atMax;
                cov_timespread = (sample_i - template_i) * (sample_j - template_j);
                closest_pts_i = two_closest_points(i * NFREQ + tshift, x)[1]
                closest_pts_j = two_closest_points(j * NFREQ + tshift, x)[1]
                sampleerr_i = np.average(ey[closest_pts_i])
                sampleerr_j = np.average(ey[closest_pts_j])
                cov_template = sampleerr_i * sampleerr_j
                covmat[i,j] += math.hypot(cov_timespread,cov_template)
                
    covmat /= n_events
    print("covmat = ",covmat)

    corrmat = correlation_from_covariance(covmat)
    print("corrmat = ",corrmat)

    fout.cd()
    covmat_th2 = numpy_to_TH2F(covmat)
    corrmat_th2 = numpy_to_TH2F(corrmat,"PulseCorrelation","Template Correlation Matrix")
    covmat_th2.Write()
    corrmat_th2.Write()
    fout.Close()
    
    
    # plot the largest positive and negative variations, to have an idea
    # maxs = sorted(t0s)[-6:-1]
    # mins = sorted(t0s)[0:6]
    # print("max spreads = ",mins+maxs)
    
    # canvas = ROOT.TCanvas("c1", "Pulse Variations", 800, 600)
    # legend = ROOT.TLegend(0.65, 0.5, 0.9, 0.9)
    # legend.SetHeader("Pulse variations", "C")
    
    # colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kMagenta, ROOT.kOrange + 1]
    # graphs = []

    # for i,s in enumerate(mins+maxs):
    #     ps_params[1] = s
    #     pulse = ROOT.TF1(f"pulse{i}", pyf_total, 0, 60., ps_params.size)
    #     pulse.SetParameters(*ps_params)
        
    #     pulse.SetLineColor(colors[i % len(colors)])
    #     pulse.SetLineWidth(2)
    #     pulse.SetName("")
    #     pulse.GetXaxis().SetTitle("time [ns]")
    #     pulse.GetYaxis().SetTitle("p.d.f.")
    #     if i == 0:
    #         pulse.Draw()  # first one defines axes
    #     else:
    #         pulse.Draw("SAME")
    #     legend.AddEntry(pulse, f"Event {i} (t0 shift={s*4:.3f} ns)", "l")
    #     graphs.append(pulse)

    # legend.Draw()
    # canvas.Update()
    # canvas.SaveAs("pulses_spread.root")
    # canvas.SaveAs("pulses_spread.pdf")
    
if __name__ == "__main__":
    #makeInputPSGraphFromIdealFunction()
    #makeInputPSGraphFromTestBeam("data/pulse_tgraph_fromTB.root")    
    makePulseCovariance("data/EmptyFileTestBeamPhase2.root","PulseShape/grPulseShape",n_events=100)
