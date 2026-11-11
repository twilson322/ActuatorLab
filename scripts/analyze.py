#!/usr/bin/env python3
"""Turns the CSV logs in results/ into figures and markdown tables.

    python3 scripts/analyze.py --results results
"""
import argparse
import pathlib

import numpy as np
import pandas as pd
from scipy import signal

import plotstyle
from plotstyle import CONTROLLER_COLOR as COL

plotstyle.apply()
import matplotlib.pyplot as plt  # noqa: E402

FS = 5000.0
LABEL = {
    "motor_only": "motor encoder only",
    "motor_defl": "motor encoder + deflection FF",
    "motor_defl_velnotch": "+ notch in speed loop",
    "motor_defl_dob": "+ DOB in speed loop",
    "dual_pi_notch": "dual encoder PI + notch",
    "dual_pi_nonotch": "dual encoder PI, no notch",
    "dual_rls": "dual encoder + RLS",
    "adrc_dual": "ADRC",
    "stsmc_dual": "super-twisting SMC",
    "lqg": "LQG (dual encoder)",
    "lqg_nodist": "LQG without load state",
    "lqg_rls": "LQG + RLS re-synthesis",
}


def md_table(df, floatfmt=".2f"):
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            cells.append(f"{v:{floatfmt}}" if isinstance(v, (float, np.floating)) else str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


# --------------------------------------------------------------------------- tracking
def tracking(res):
    show = ["motor_only", "motor_defl", "dual_pi_notch", "stsmc_dual", "lqg"]
    fig, ax = plt.subplots(2, 1, figsize=(7.4, 6.2))
    for name in show:
        d = pd.read_csv(res / f"track_{name}.csv")
        ax[0].plot(d.t, 1e3 * d.err, color=COL[name], lw=1.2, label=LABEL[name])
        z = (d.t > 1.55) & (d.t < 2.25)
        ax[1].plot(d.t[z], 1e3 * d.err[z], color=COL[name], lw=1.4, label=LABEL[name])
    for a in ax:
        a.axvspan(1.6, 2.2, color="#e4e3dd", alpha=0.6, lw=0)
        a.axhline(0, color=plotstyle.MUTED, lw=0.6)
        a.set_ylabel("link error [mrad]")
    ax[0].set_title("Link tracking error: 90 deg and -45 deg quintic moves, then a 25 Nm load step")
    ax[0].set_ylim(-14, 14)
    ax[0].legend(ncol=2, fontsize=8, loc="lower left")
    ax[1].set_title("Load step at t = 1.6 s (shaded while applied)")
    ax[1].set_xlabel("time [s]")
    fig.tight_layout()
    fig.savefig(res / "track_error.png")
    plt.close(fig)

    m = pd.read_csv(res / "track_metrics.csv")
    m.insert(1, "description", [LABEL[c] for c in m.controller])
    m["dist_recovery_ms"] = [f"{v:.0f}" if v < 599 else ">600" for v in m.dist_recovery_ms]
    keep = m[["controller", "description", "rms_err_motion_mrad", "max_err_motion_mrad",
              "dist_peak_mrad", "dist_recovery_ms", "peak_iq_A", "chatter_Nm"]]
    keep.columns = ["config", "what it is", "RMS err [mrad]", "max err [mrad]",
                    "load peak [mrad]", "recovery to 2 mrad [ms]", "peak iq [A]", "torque chatter [Nm]"]
    return md_table(keep, ".3f")


# --------------------------------------------------------------------------- loop gains
def loop_gain(res, name):
    d = pd.read_csv(res / f"loopgain_{name}.csv")
    m = (d.t > 1.2) & (d.t < 30.8)
    r = d.tau_inj[m].to_numpy()
    u = d.tau_cmd[m].to_numpy()
    c = u - r
    f, Prr = signal.welch(r, FS, nperseg=2**14)
    _, Puu = signal.welch(u, FS, nperseg=2**14)
    _, Pru = signal.csd(r, u, FS, nperseg=2**14)
    _, Prc = signal.csd(r, c, FS, nperseg=2**14)
    L = -Prc / Pru  # c = -L u, with the injection as instrument
    coh = np.abs(Pru) ** 2 / (Prr * Puu)
    ok = (f > 1.0) & (f < 800.0) & (coh > 0.8)
    return f[ok], L[ok]


def margins(f, L):
    mag, ph = np.abs(L), np.angle(L, deg=True)
    pms, gms, fx = [], [], []
    for i in range(len(f) - 1):
        if mag[i] >= 1.0 > mag[i + 1]:
            pms.append(180.0 + ph[i])
            fx.append(f[i])
        # phase crossing of +/-180 deg (angle() wraps there)
        if abs(ph[i] - ph[i + 1]) > 300.0 and mag[i] < 1.0:
            gms.append(-20 * np.log10(mag[i]))
    pms = [((p + 180) % 360) - 180 for p in pms]
    Ms = np.max(np.abs(1.0 / (1.0 + L)))
    return (min(pms) if pms else np.nan, min(gms) if gms else np.inf,
            fx[int(np.argmin(pms))] if pms else np.nan, Ms)


def loopgains(res):
    names = ["motor_defl", "motor_defl_velnotch", "dual_pi_notch", "adrc_dual", "stsmc_dual", "lqg"]
    rows = []
    fig, ax = plt.subplots(2, 1, figsize=(7.4, 6.0), sharex=True)
    for n in names:
        f, L = loop_gain(res, n)
        pm, gm, fc, Ms = margins(f, L)
        rows.append({"config": n, "crossover [Hz]": fc, "phase margin [deg]": pm,
                     "gain margin [dB]": gm, "peak sensitivity Ms": Ms,
                     "Ms [dB]": 20 * np.log10(Ms)})
        if n in ("motor_defl", "motor_defl_velnotch", "dual_pi_notch", "lqg"):
            ax[0].semilogx(f, 20 * np.log10(np.abs(L)), color=COL[n], lw=1.2, label=LABEL[n])
            ax[1].semilogx(f, np.degrees(np.unwrap(np.angle(L))), color=COL[n], lw=1.2)
    ax[0].axhline(0, color=plotstyle.MUTED, lw=0.8)
    ax[0].set_ylabel("|L| [dB]")
    ax[0].set_title("Measured loop gain, broken at the motor torque input")
    ax[0].legend(fontsize=8)
    ax[1].set_ylabel("phase [deg]")
    ax[1].set_xlabel("frequency [Hz]")
    fig.tight_layout()
    fig.savefig(res / "loop_gain.png")
    plt.close(fig)
    df = pd.DataFrame(rows)
    df.to_csv(res / "margins.csv", index=False)
    return md_table(df, ".1f")


# --------------------------------------------------------------------------- tuning grid
def tuning(res):
    t = pd.read_csv(res / "tuning.csv")
    fig, ax = plt.subplots(1, 2, figsize=(8.0, 3.4), sharey=True)
    grids = {}
    for a, fam, title in [(ax[0], "dual_pi_nonotch", "link PI, no notch"),
                          (ax[1], "dual_pi_notch", f"link PI + notch")]:
        x = t[t.family == fam].copy()
        x["kp"] = x.params.str.extract(r"link.kp=([\d.]+)").astype(float)
        x["ki"] = x.params.str.extract(r"link.ki=([\d.]+)").astype(float)
        g = x.pivot(index="kp", columns="ki", values="rms_mrad")
        grids[fam] = g
        im = a.imshow(np.log10(g.values), origin="lower", cmap="Blues", vmin=-0.4, vmax=2.0, aspect="auto")
        a.set_xticks(range(len(g.columns)), [f"{v:g}" for v in g.columns])
        a.set_yticks(range(len(g.index)), [f"{v:g}" for v in g.index])
        a.set_xlabel("link ki [1/s]")
        a.set_title(title)
        a.grid(False)
        for i in range(g.shape[0]):
            for j in range(g.shape[1]):
                v = g.values[i, j]
                a.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=7,
                       color="white" if np.log10(v) > 1.0 else plotstyle.INK)
    ax[0].set_ylabel("link kp")
    fig.suptitle("RMS tracking error [mrad] over the link-loop gain grid", fontsize=11)
    fig.tight_layout()
    fig.savefig(res / "link_gain_grid.png")
    plt.close(fig)
    nn, wn = grids["dual_pi_nonotch"], grids["dual_pi_notch"]
    ok_nn = int((nn.values < 1.0).sum())
    ok_wn = int((wn.values < 1.0).sum())
    return ok_nn, ok_wn, nn.size


# --------------------------------------------------------------------------- Monte Carlo
def montecarlo(res):
    mc = pd.read_csv(res / "montecarlo.csv")
    order = ["motor_only", "motor_defl", "motor_defl_velnotch", "motor_defl_dob", "dual_pi_notch",
             "adrc_dual", "stsmc_dual", "lqg", "lqg_nodist"]
    rows = []
    for n in order:
        x = mc[mc.controller == n]
        s = x[x.stable == 1]
        # "Well behaved": stable, below 2 mrad RMS during motion and no sustained oscillation.
        good = s[(s.rms_err_motion_mrad < 2.0) & (s.chatter_Nm < 0.08)]
        rows.append({
            "config": n,
            "stable [%]": 100.0 * len(s) / len(x),
            "well behaved [%]": 100.0 * len(good) / len(x),
            "RMS median [mrad]": s.rms_err_motion_mrad.median(),
            "RMS p90 [mrad]": s.rms_err_motion_mrad.quantile(0.9),
            "load peak median [mrad]": s.dist_peak_mrad.median(),
            "recovered within 600 ms [%]": 100.0 * (s.dist_recovery_ms < 599).mean(),
        })
    df = pd.DataFrame(rows)

    fig, ax = plt.subplots(1, 2, figsize=(8.6, 4.0))
    show = ["motor_defl", "dual_pi_notch", "adrc_dual", "stsmc_dual", "lqg"]
    for a, col, title in [(ax[0], "rms_err_motion_mrad", "RMS tracking error [mrad]"),
                          (ax[1], "dist_peak_mrad", "peak deviation under 25 Nm load [mrad]")]:
        data = [mc[(mc.controller == n) & (mc.stable == 1)][col].to_numpy() for n in show]
        bp = a.boxplot(data, orientation="horizontal", widths=0.55, patch_artist=True, showfliers=True,
                       flierprops={"markersize": 3, "markeredgewidth": 0.6})
        for patch, n in zip(bp["boxes"], show):
            patch.set_facecolor(COL[n])
            patch.set_alpha(0.75)
            patch.set_edgecolor(plotstyle.INK)
        for med in bp["medians"]:
            med.set_color(plotstyle.INK)
        a.set_xscale("log")
        a.set_title(title)
        a.set_yticks(range(1, len(show) + 1), [LABEL[n] for n in show] if a is ax[0] else [""] * len(show))
    fig.suptitle("200 randomized plants (link inertia 0.5x-3x, stiffness +/-30%, friction 0.5x-2x, gravity 0.6x-1.6x)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(res / "montecarlo.png")
    plt.close(fig)
    df.to_csv(res / "montecarlo_summary.csv", index=False)

    # Where does LQG lose? Correlate failure with the inertia ratio.
    lq = mc[mc.controller == "lqg"]
    worst = lq.sort_values("rms_err_motion_mrad", ascending=False).head(5)[["Jl", "k", "mgl", "rms_err_motion_mrad", "stable"]]
    return md_table(df, ".1f"), worst


# --------------------------------------------------------------------------- adaptation
def adaptation(res):
    fig, ax = plt.subplots(3, 1, figsize=(7.4, 6.6), sharex=True)
    for n in ["lqg", "lqg_rls"]:
        d = pd.read_csv(res / f"adapt_{n}.csv")
        ax[0].plot(d.t, 1e3 * d.err, color=COL[n], lw=1.0, label=LABEL[n])
    d = pd.read_csv(res / "adapt_lqg_rls.csv")
    ax[1].plot(d.t, d.d1, color=COL["lqg_rls"], label="estimated")
    ax[1].step([0, 4, 8], [0.6, 1.5, 1.5], where="post", color=plotstyle.INK, lw=1.0, ls="--", label="true")
    ax[2].plot(d.t, d.d2, color=COL["lqg_rls"], label="estimated")
    ax[2].step([0, 4, 8], [11.8, 24.0, 24.0], where="post", color=plotstyle.INK, lw=1.0, ls="--", label="true")
    for a in ax:
        a.axvline(4.0, color=plotstyle.MUTED, lw=0.8, ls=":")
    ax[0].set_ylabel("link error [mrad]")
    ax[0].set_title("Payload change at t = 4 s: link inertia 0.6 -> 1.5 kg m^2, gravity 11.8 -> 24 Nm")
    ax[0].legend(fontsize=8, loc="lower left")
    ax[1].set_ylabel("J_link [kg m^2]")
    ax[1].legend(fontsize=8, loc="upper left")
    ax[2].set_ylabel("m g l [Nm]")
    ax[2].set_xlabel("time [s]")
    fig.tight_layout()
    fig.savefig(res / "adaptation.png")
    plt.close(fig)
    a = pd.read_csv(res / "adapt_metrics.csv")
    a["rms_err_mrad"] = a.rms_err_mrad.map(lambda v: "" if pd.isna(v) else f"{v:.3f}")
    a["diverged_at_s"] = a.diverged_at_s.map(lambda v: "" if pd.isna(v) else f"{v:.2f}")
    return md_table(a)


# --------------------------------------------------------------------------- speed-loop sweep
def sweep(res):
    s = pd.read_csv(res / "sweep.csv")
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    for notch, name in [(0, "motor_defl"), (1, "motor_defl_velnotch")]:
        x = s[(s.notch == notch) & (s.stable == 1)]
        ax.plot(x.bw_hz, x.rms_err_motion_mrad, marker="o", ms=4, color=COL[name], label=LABEL[name])
    ax.set_yscale("log")
    ax.set_xlabel("speed-loop design bandwidth [Hz] (rigid-body)")
    ax.set_ylabel("RMS link error [mrad]")
    ax.set_title("Collocated speed loop: bandwidth sweep")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(res / "speed_loop_sweep.png")
    plt.close(fig)
    unstable = s[s.stable == 0]
    return unstable


# --------------------------------------------------------------------------- RT bench
def rt(res):
    out = []
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    for fn, lab, color in [("rt_bench.csv", "idle machine", plotstyle.SERIES[0]),
                           ("rt_bench_stress.csv", "2 busy threads", plotstyle.SERIES[1])]:
        p = res / fn
        if not p.exists():
            continue
        d = pd.read_csv(p)
        bins = np.geomspace(1, 2e4, 60)
        ax.hist(d.wake_us.clip(lower=1), bins=bins, histtype="step", lw=1.5, color=color, label=lab)
        out.append((lab, d.wake_us.quantile(0.5), d.wake_us.quantile(0.99), d.wake_us.quantile(0.999),
                    d.wake_us.max(), int((d.wake_us > 1000).sum()), len(d)))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axvline(1000, color=plotstyle.MUTED, lw=0.8, ls=":")
    ax.set_xlabel("wake-up latency of the 1 kHz host thread [us]")
    ax.set_ylabel("cycles")
    ax.set_title("SCHED_FIFO on a stock (non-PREEMPT_RT) VM kernel")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(res / "rt_latency.png")
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    res = pathlib.Path(ap.parse_args().results)
    parts = ["# Results (generated by scripts/analyze.py)\n"]
    parts.append("## Nominal plant: tracking and load rejection\n\n" + tracking(res))
    parts.append("## Measured loop margins\n\n" + loopgains(res))
    ok_nn, ok_wn, tot = tuning(res)
    parts.append(f"## Link-loop gain grid\n\nGrid points with RMS error below 1 mrad: {ok_nn}/{tot} "
                 f"without the notch, {ok_wn}/{tot} with it.")
    mc_table, worst = montecarlo(res)
    parts.append("## Monte Carlo over 200 plants\n\n" + mc_table)
    parts.append("Worst five LQG runs:\n\n" + md_table(worst, ".3f"))
    parts.append("## Payload change\n\n" + adaptation(res))
    uns = sweep(res)
    parts.append(f"## Speed-loop sweep\n\nUnstable points: {len(uns)}\n\n" + (md_table(uns) if len(uns) else ""))
    rows = rt(res)
    if rows:
        df = pd.DataFrame(rows, columns=["condition", "p50 [us]", "p99 [us]", "p99.9 [us]", "max [us]",
                                         "late > 1 ms", "cycles"])
        parts.append("## 1 kHz loop timing\n\n" + md_table(df, ".1f"))
    (res / "summary.md").write_text("\n\n".join(parts) + "\n")
    print((res / "summary.md").read_text())


if __name__ == "__main__":
    main()
