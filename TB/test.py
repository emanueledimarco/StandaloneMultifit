import numpy as np
import uproot

from pulse_reco import PiecewiseCubicSplineNP, fit_channels



def split(waveforms, threshold=None, pre=5, post=10, baseline_samples=10):

    # Assume waveforms is shape (E, C, S)
    E, C, S = waveforms.shape

    if threshold is not None:
      argmax_idx = np.argmax(waveforms > threshold, axis=2)  # shape (E, C)
    else:
      argmax_idx = np.argmax(waveforms, axis=2)  # shape (E, C)

    # Step 2: Build offsets
    window_offsets = np.arange(-int(pre), int(post)).reshape(1, 1, -1)         # shape (1,1,15)
    baseline_offsets = np.arange(-int(pre)-int(baseline_samples), -int(pre)).reshape(1, 1, -1)      # shape (1,1,10)

    # Expand argmax index for broadcasting
    argmax_exp = argmax_idx[:, :, np.newaxis]  # shape (E, C, 1)

    # Add offsets and wrap with modulo S to stay in bounds
    window_indices   = (argmax_exp + window_offsets) % S        # shape (E, C, 15)
    baseline_indices = (argmax_exp + baseline_offsets) % S      # shape (E, C, 10)

    # Build broadcasted event/channel indices
    event_idx = np.arange(E)[:, None, None]
    chan_idx  = np.arange(C)[None, :, None]

    # Extract waveform windows and baseline windows
    window_waveforms   = waveforms[event_idx, chan_idx, window_indices]      # (E, C, 15)
    baseline_waveforms = waveforms[event_idx, chan_idx, baseline_indices]    # (E, C, 10)

    # Step 3: Compute baseline mean
    baseline = np.mean(baseline_waveforms, axis=2)       # shape (E, C)
    baseline_std = np.std(baseline_waveforms, axis=2)    # shape (E, C)
    baseline_integral = np.sum(baseline_waveforms, axis=2)  # shape (E, C)

    return argmax_idx, baseline, baseline_std, baseline_integral, window_waveforms


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

    signal_samples_pre_peak = 3
    signal_samples_post_peak = 7

    charge_zerosup_peak_threshold = 10

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

    channels = np.array([5,6,13,14,15,16,109,110,119])
    waves = waves[:, channels, :]   # (E, Csel, S)

    print(waves.shape)

    max_idx, baselines, baselines_std, baseline_integral, signal_window = split(waves, pre=signal_samples_pre_peak, post=signal_samples_post_peak)
    # print(baseline_integral.shape)

    waves = waves - np.repeat(baselines[:, :, np.newaxis], waves.shape[2], axis=2)  # baseline subtraction
    signal_window = signal_window - np.repeat(baselines[:, :, np.newaxis], signal_window.shape[2], axis=2) # baseline subtraction in the signal window


    charge = np.sum(signal_window, axis=2)

    t_grid = np.arange(signal_window.shape[2]) * Ts

    pulse_values = pulse.eval(t_grid)
    t_pulse_peak = t_grid[np.argmax(pulse_values)]


    print("about to enter fit_channels")

    amp, dt = fit_channels(
        signal_window,
        pulse,
        t_grid,
        signal_samples_pre_peak,  # data peak (scalar)
        t_pulse_peak     # model peak (scalar)
    )

    fit_time = dt + np.ones(dt.shape)*signal_samples_pre_peak

    print("fit_time: ", fit_time)
    print("fit_amp: ", amp)

    # -------------------------
    # Output tree
    # -------------------------
    with uproot.recreate(output_file) as fout:

        fout["T"] = {
            "event": np.arange(fit_time.shape[0]),
            "fit_amp": amp,
            "fit_time": fit_time,
            "charge": charge,
            "peak_idx": max_idx,
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
