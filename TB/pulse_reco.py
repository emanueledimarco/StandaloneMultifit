
import numpy as np
import ROOT
import sys

# ============================================================
# 1. Piecewise Cubic Spline (C++-equivalent, vectorized eval)
# ============================================================

class PiecewiseCubicSplineNP:
    """
    Exact translation of your C++ PiecewiseCubicSpline:
    - segment-based cubic polynomials
    - vectorized evaluation
    """

    def __init__(self, segments, Ts):
        self.x0 = np.array([s["x0"] for s in segments])
        self.x1 = np.array([s["x1"] for s in segments])
        self.xc = np.array([s["xc"] for s in segments])

        self.a = np.array([s["a"] for s in segments])
        self.b = np.array([s["b"] for s in segments])
        self.c = np.array([s["c"] for s in segments])
        self.d = np.array([s["d"] for s in segments])
        self.Ts = Ts

    def _select(self, x):

        x = np.asarray(x)[None, ...]  # (1, E, C, N)

        x0 = self.x0[:, None, None, None]
        x1 = self.x1[:, None, None, None]

        inside = (x >= x0) & (x <= x1)

        idx = np.argmax(inside, axis=0)  # (E, C, N)

        xc = self.xc[idx]
        dx = x[0] - xc
        #print("x", x, "dx", dx)
        mask = np.abs(dx) <= (self.Ts / 2.0)

        return idx, dx, mask


    def eval(self, x):

        idx, dx, mask = self._select(x)

        a = self.a[idx]
        b = self.b[idx]
        c = self.c[idx]
        d = self.d[idx]

        f = a + b*dx + c*dx*dx + d*dx*dx*dx

        return np.where(mask, f, 0.0)



    def derivative(self, x):

        idx, dx, mask = self._select(x)

        b = self.b[idx]
        c = self.c[idx]
        d = self.d[idx]

        deriv = b + 2*c*dx + 3*d*dx*dx

        return np.where(mask, deriv, 0.0)


# ============================================================
# 2. Iterative pulse fit (ROOT-equivalent Gauss-Newton)
# ============================================================

def fit_pulse_iterative(waveforms, pulse, t, t_data_peak, t_template_peak, n_iter=4):

    E, C, N = waveforms.shape

    # initial alignment from DATA only
    dt = np.full((E, C), t_data_peak - t_template_peak)   # (E, C)
    A  = np.ones((E, C))

    #print("t", t, "t_template_peak,", t_template_peak, "t_centered", t_centered)

    # DEBUG (optional)
    #p = pulse.eval(t_centered[None, None, :])
    #print(p)
    #sys.exit(0)

    for _ in range(n_iter):

        # -------------------------
        # build shifted time grid
        # -------------------------
        # (E, C, N)
        t_shift = t[None, None, :] - dt[:, :, None]

        # -------------------------
        # evaluate pulse
        # -------------------------
        P  = pulse.eval(t_shift)        # (E, C, N)
        dP = pulse.derivative(t_shift) # (E, C, N)

        #print("p", ",".join(map(str, P[1, 3, :])))
        #print("wf", ",".join(map(str, waveforms[1, 3, :])))
        #sys.exit(0)

        # -------------------------
        # projections (sum over samples)
        # -------------------------
        Ap = np.sum(waveforms * P, axis=2)   # (E, C)
        Ad = np.sum(waveforms * dP, axis=2)

        PP   = np.sum(P * P, axis=2)
        PdP  = np.sum(P * dP, axis=2)
        dPdP = np.sum(dP * dP, axis=2)

        # -------------------------
        # solve 2x2 system
        # -------------------------
        denom = PP * dPdP - PdP * PdP
        denom = np.clip(denom, 1e-12, None)

        A_new = (Ap * dPdP - Ad * PdP) / denom
        Ccorr = (Ad * PP - Ap * PdP) / denom

        # -------------------------
        # update
        # -------------------------
        dt += (-Ccorr / np.clip(A_new, 1e-12, None))
        A = A_new

    return A, dt


def fit_channels(waveforms, pulse, t, t_data_peak, t_pulse_peak):
    """
    waveforms: (E, C, N)
    peak_idx: (E, C)

    returns:
        amp (E, C)
        time (E, C)
    """

    E, C, N = waveforms.shape

    # flatten channels into events for vectorized processing
    A, dt = fit_pulse_iterative(waveforms, pulse, t, t_data_peak, t_pulse_peak, n_iter=4)

    A = A.reshape(E, C)
    dt = dt.reshape(E, C)

    return A, dt
