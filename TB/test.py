import numpy as np
import uproot
import awkward as ak

from pulse_reco import PiecewiseCubicSplineNP, fit_channels


# ============================================================
# 1. Load spline segments (C++ header style txt file)
# ============================================================

def load_segments_txt(filename):
    """
    Format expected (same idea as your C++ loader):

    xc x0 x1 a b c d
    """

    segs = []

    with open(filename, "r") as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue

            vals = line.split()

            segs.append({
                "xc": float(vals[0]),
                "x0": float(vals[1]),
                "x1": float(vals[2]),
                "a":  float(vals[3]),
                "b":  float(vals[4]),
                "c":  float(vals[5]),
                "d":  float(vals[6]),
            })

    return segs


# ============================================================
# 2. Peak finding (vectorized per event/channel)
# ============================================================

def find_peaks(waves):
    """
    waves: (E, C, S)
    """

    return np.argmax(waves, axis=2)


# ============================================================
# 3. Main reconstruction (ROOT-equivalent logic)
# ============================================================

def run_reco(input_file, spline_file, output_file):

    Ts = 6.25
    t_ref = -100

    # -------------------------
    # Load spline
    # -------------------------
    segs = load_segments_txt(spline_file)
    pulse = PiecewiseCubicSplineNP(segs, Ts)

    # -------------------------
    # Load ROOT file
    # -------------------------
    f = uproot.open(input_file)
    tree = f["h4"]

    waves = tree["xtal_sample"].array(library="np")  # (E, C, S)

    E, C, S = waves.shape

    channels = np.array([5,6,13,14,15,16,109,110,119])
    waves = waves[:, channels, :]   # (E, Csel, S)

    peak_idx = find_peaks(waves)   # (E, C)

    t_peak = peak_idx * Ts         # (E, C)

    j_start = np.clip(peak_idx - 3, 0, S-1)
    j_end   = np.clip(peak_idx + 7, 0, S-1)

    ped_start = j_start - 11   # 10 samples before signal window
    ped_end   = j_start - 1

    ped_start = np.clip(ped_start, 0, S-1)
    ped_end   = np.clip(ped_end,   0, S-1)

    idx = np.arange(S)[None, None, :]

    ped_mask = (idx >= ped_start[:, :, None]) & (idx <= ped_end[:, :, None])
    sig_mask = (idx >= j_start[:, :, None]) & (idx <= j_end[:, :, None])

    pedestal = (
        np.sum(waves * ped_mask, axis=2) /
        np.maximum(ped_mask.sum(axis=2), 1)
    )

    waves_corr = waves - pedestal[:, :, None]


    charge = np.sum(waves_corr * sig_mask, axis=2)

    t_grid = np.arange(S) * Ts


    pulse_values = pulse.eval(t_grid)
    t_pulse_peak = t_grid[np.argmax(pulse_values)]


    print("about to enter fit_channels")

    amp, dt = fit_channels(
        waves_corr,
        pulse,
        t_grid,
        t_peak,          # data peak (E, C)
        t_pulse_peak     # model peak (scalar)
    )

    fit_time = t_peak + dt

    print("fit_time: ", fit_time)
    print("fit_amp: ", amp)
    # -------------------------
    # Output tree
    # -------------------------
    with uproot.recreate(output_file) as fout:

        fout["T"] = {
            "event": np.arange(E),
            "fit_amp": amp,
            "fit_time": fit_time,
            "charge": charge,
            "peak_idx": peak_idx,
        }


# ============================================================
# 4. Run
# ============================================================

if __name__ == "__main__":
    run_reco(
        input_file="19535_0002_unpacked.root",
        spline_file="coeffs_global.txt",
        output_file="results.root"
    )
