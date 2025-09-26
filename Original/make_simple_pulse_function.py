import ROOT
import math
import numpy as np

def pyf_tf1_coulomb(x, p):
    return p[1] * x[0] * x[1] / (p[0]**2) * math.exp(-p[2] / p[0])


def pyf_f(t, t0, tau_r):
    return 1./(math.exp(-1.*(t-t0)/tau_r) + 1)

# fast decay time (fluorescence)
def pyf_g(t, t0, tau_f):
    return 1./(math.exp((t-t0)/tau_f) + 1)

# slow decay time (phosphorence)
def pyf_h(t, t0, tau_s):
    return 1./(math.exp((t-t0)/tau_s) + 1)

def pyf_total(t, p):
    # p[0] = A = pulse amplitude
    # p[1] = t0 = pulse arrival time
    # p[2] = tau_r = rise time (~1 ns)
    # p[3] = tau_f = fast decay time (~5 ns) - fluorescence
    # p[4] = tau_s = fast decay time (~50 ns) - phosphorence
    # p[5] = R = relative weight of the decay times (R=1 means slow decay component is 0. In any case 0 <= R <=1)
    return p[0] * pyf_f(t,p[1],p[2]) * (p[5] * pyf_g(t,p[1],p[3]) + (1-p[5]) * pyf_h(t,p[1],p[4])) 


