#!/usr/bin/env python3
"""Frequency-response identification of the compliant joint.

Input:  results/chirp_plant.csv  (log chirp added to the motor torque command, 5 kHz log)
Output: results/identified.txt   (key=value, read by the C++ controllers for notch placement)
        results/bode_plant.png

The FRF tau_cmd -> w_m_hat uses the drive's own speed estimate (the true speed is never used).
Because the soft position hold closes a loop around the plant, the plain H1 estimate
S_uy / S_uu is biased toward -1/C wherever sensor noise drives the controller harder than the
chirp does (it shows up as a flat line at 1/Kp above ~450 Hz). Using the known injection r as
an instrument, H = S_ry / S_ru, removes that bias. A two-mass model with a pure delay is fitted in the complex-log domain:

    H(s) = (Jl s^2 + d s + k) / ( s [ Jm Jl s^2 + (Jm + Jl) d s + (Jm + Jl) k ] ) * exp(-s Td)

with all quantities reflected to the motor side.
"""
import argparse
import pathlib

import numpy as np
import pandas as pd
from scipy import optimize, signal

import plotstyle

N_GEAR = 50.0


def two_mass(p, w):
    Jm, Jl, k, d, Td = p
    s = 1j * w
    num = Jl * s**2 + d * s + k
    den = s * (Jm * Jl * s**2 + (Jm + Jl) * d * s + (Jm + Jl) * k)
    return num / den * np.exp(-s * Td)


def estimate_frf(df, fs=5000.0, t0=1.0, t1=25.0, nperseg=2**14):
    m = (df.t > t0) & (df.t < t1)
    r = df.tau_inj[m].to_numpy()
    u = df.tau_cmd[m].to_numpy()
    y = df.w_m_hat[m].to_numpy()
    f, Prr = signal.welch(r, fs, nperseg=nperseg)
    _, Pyy = signal.welch(y, fs, nperseg=nperseg)
    _, Pry = signal.csd(r, y, fs, nperseg=nperseg)
    _, Pru = signal.csd(r, u, fs, nperseg=nperseg)
    H = Pry / Pru
    coh = np.abs(Pry) ** 2 / (Prr * Pyy)  # how much of y the injection explains
    return f, H, coh


def fit(f, H, coh, fmin=3.0, fmax=300.0, npts=240):
    # Resample on a log grid so the fit is not dominated by the dense high-frequency bins.
    fg = np.geomspace(fmin, fmax, npts)
    Hg = np.interp(fg, f, H.real) + 1j * np.interp(fg, f, H.imag)
    cg = np.interp(fg, f, coh)
    wgt = np.sqrt(np.clip(cg, 0.0, 1.0))
    w = 2 * np.pi * fg

    def resid(x):
        Hm = two_mass(np.exp(x), w)
        r = np.log(Hm / Hg)  # real part = log-magnitude error, imag part = phase error
        return np.concatenate([wgt * r.real, wgt * np.angle(np.exp(1j * r.imag))])

    # Deliberately rough initial guess (order of magnitude only).
    x0 = np.log([1e-4, 1e-4, 1.0, 1e-3, 1e-4])
    sol = optimize.least_squares(resid, x0, method="lm", max_nfev=20000)
    return np.exp(sol.x), sol


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    a = ap.parse_args()
    out = pathlib.Path(a.results)
    df = pd.read_csv(out / "chirp_plant.csv")
    f, H, coh = estimate_frf(df)
    p, sol = fit(f, H, coh)
    Jm, Jl, k, d, Td = p
    f_ar = np.sqrt(k / Jl) / (2 * np.pi)
    f_r = np.sqrt(k * (1 / Jm + 1 / Jl)) / (2 * np.pi)

    truth = {"Jm": 6e-5, "Jl": 0.6 / N_GEAR**2, "k": 4000 / N_GEAR**2, "d": 2.0 / N_GEAR**2}
    print("parameter      identified       true (simulator)   error")
    for name, val in zip(["Jm", "Jl", "k", "d"], [Jm, Jl, k, d]):
        print(f"  {name:4s}  {val:14.4e}  {truth[name]:14.4e}   {100 * (val / truth[name] - 1):+6.1f} %")
    f_ar_true = np.sqrt(truth["k"] / truth["Jl"]) / (2 * np.pi)
    f_r_true = np.sqrt(truth["k"] * (1 / truth["Jm"] + 1 / truth["Jl"])) / (2 * np.pi)
    print(f"  anti-resonance {f_ar:.2f} Hz (true {f_ar_true:.2f}),  resonance {f_r:.2f} Hz (true {f_r_true:.2f})")
    print(f"  loop delay {1e6 * Td:.0f} us,  fit cost {sol.cost:.3e}")

    with open(out / "identified.txt", "w") as fh:
        fh.write("# written by scripts/identify.py from results/chirp_plant.csv\n")
        fh.write(f"f_antires_hz={f_ar:.4f}\nf_res_hz={f_r:.4f}\n")
        fh.write(f"Jm={Jm:.6e}\nJl_out={Jl * N_GEAR**2:.6f}\nk_out={k * N_GEAR**2:.3f}\n")
        fh.write(f"d_out={d * N_GEAR**2:.4f}\ndelay_us={1e6 * Td:.2f}\n")

    plotstyle.apply()
    import matplotlib.pyplot as plt

    # Only show bins the injection actually explains. Above ~350 Hz the 400 Hz speed observer
    # rolls the response off into quantization noise and the estimate is meaningless.
    sel = (f > 1.5) & (f < 600) & (coh > 0.8)
    Hfit = two_mass(p, 2 * np.pi * f[sel])
    fig, ax = plt.subplots(3, 1, figsize=(7.2, 7.0), sharex=True,
                           gridspec_kw={"height_ratios": [3, 2, 1.2]})
    ax[0].semilogx(f[sel], 20 * np.log10(np.abs(H[sel])), color=plotstyle.SERIES[0], label="measured (chirp as instrument)")
    ax[0].semilogx(f[sel], 20 * np.log10(np.abs(Hfit)), color=plotstyle.SERIES[1], ls="--", label="two-mass fit")
    for fx, lab, dx, ha in [(f_ar, f"anti-resonance {f_ar:.1f} Hz", -4, "right"),
                            (f_r, f"resonance {f_r:.1f} Hz", 4, "left")]:
        ax[0].axvline(fx, color=plotstyle.MUTED, lw=0.8, ls=":")
        ax[0].annotate(lab, (fx, ax[0].get_ylim()[1]), xytext=(dx, -12), textcoords="offset points",
                       fontsize=8, color=plotstyle.MUTED, ha=ha)
    ax[0].set_ylabel("|w_m / tau_m|  [dB rad/s/Nm]")
    ax[0].legend(loc="lower left")
    ax[0].set_title("Motor torque to motor speed, identified from a log chirp")
    ax[1].semilogx(f[sel], np.degrees(np.angle(H[sel])), color=plotstyle.SERIES[0])
    ax[1].semilogx(f[sel], np.degrees(np.angle(Hfit)), color=plotstyle.SERIES[1], ls="--")
    ax[1].set_ylabel("phase [deg]")
    allf = (f > 1.5) & (f < 600)
    ax[2].semilogx(f[allf], coh[allf], color=plotstyle.SERIES[0])
    ax[2].axhline(0.8, color=plotstyle.MUTED, lw=0.8, ls=":")
    ax[2].set_ylabel("coherence")
    ax[2].set_ylim(0, 1.05)
    ax[2].set_xlabel("frequency [Hz]")
    fig.tight_layout()
    fig.savefig(out / "bode_plant.png")
    print(f"wrote {out / 'identified.txt'} and {out / 'bode_plant.png'}")


if __name__ == "__main__":
    main()
