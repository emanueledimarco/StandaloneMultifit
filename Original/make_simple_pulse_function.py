import ROOT
import math
import numpy as np

# dummy phase2 PS params (time units in 1/4 ns)
#                     0 1  2 3 4 5 
ps_params = np.array([1,18,3,5,3,1])

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

def pyf_total(t, p):
    # p[0] = A = pulse amplitude
    # p[1] = t0 = pulse arrival time
    # p[2] = tau_r = rise time (~1 ns)
    # p[3] = tau_f = fast decay time (~5 ns) - fluorescence
    # p[4] = tau_s = fast decay time (~50 ns) - phosphorence
    # p[5] = R = relative weight of the decay times (R=1 means slow decay component is 0. In any case 0 <= R <=1)
    return p[0] * pyf_f(t,p[1],p[2]) * (p[5] * pyf_g(t,p[1],p[3]) + (1-p[5]) * pyf_h(t,p[1],p[4])) 

def makeInputPSGraph():
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
    outFile.Close()
    

if __name__ == "__main__":
    makeInputPSGraph()

