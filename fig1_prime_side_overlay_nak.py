#!/usr/bin/env python3
   

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import mpmath as mp
from scipy.special import erfc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GAMMA_1 = "14.134725141734693790457251983562470270784257115699"


def configure_matplotlib() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "font.size": 10,
        "axes.labelsize": 10,
        "legend.fontsize": 9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
        "axes.unicode_minus": False,
    })


def prime_powers(N: int):
                                                                                 
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    primes = np.nonzero(sieve)[0]
    del sieve
    logn, lam = [], []
    k = 1
    while True:
        base = primes[primes.astype(np.float64) ** k <= N]
        if base.size == 0:
            break
        logp = np.log(base.astype(np.float64))
        logn.append(k * logp)
        lam.append(logp)
        k += 1
    return np.concatenate(logn), np.concatenate(lam)


def prime_field(s: complex, logn, lam_w) -> complex:
                                                                                
    return complex(-np.sum(lam_w * np.exp(-s * logn)))


def main() -> None:
    ap = argparse.ArgumentParser(description="Regenerate Figure 1 (prime-side check).")
    ap.add_argument("--X", type=float, default=1e4)
    ap.add_argument("--alpha", type=float, default=0.2)
    ap.add_argument("--N", type=int, default=30_000_000)
    ap.add_argument("--lam-min", type=float, default=0.25)
    ap.add_argument("--lam-max", type=float, default=3.5)
    ap.add_argument("--check-step", type=float, default=0.025,
                    help="lambda step for the max-discrepancy check")
    ap.add_argument("--plot-step", type=float, default=0.125,
                    help="lambda step of the plotted prime-side markers")
    ap.add_argument("--output-dir", type=str, default="figures")
    args = ap.parse_args()

    configure_matplotlib()
    mp.mp.dps = 30
    X, alpha, N = float(args.X), float(args.alpha), int(args.N)
    Y = math.log(X)
    cutoff_arg = (math.log(N) - Y) / (2 * alpha * math.sqrt(Y))
    if cutoff_arg < 6:
        print(f"WARNING: cutoff argument {cutoff_arg:.2f} < 6; truncation may "
              "dominate the discrepancy (manuscript rule, Eq. cutoff-rule).")

    rho0 = mp.mpc(mp.mpf("0.5"), mp.mpf(GAMMA_1))
    logn, lam = prime_powers(N)
    weight = 0.5 * erfc((logn - Y) / (2 * alpha * math.sqrt(Y)))
    lam_w = lam * weight

    def defect(l: float) -> float:
        s = rho0 + mp.mpf(l) / Y
        zl = complex(mp.zeta(s, derivative=1) / mp.zeta(s))
        P = prime_field(complex(s), logn, lam_w)
        return (l / Y) * (P - zl).real                                    

    def L_finite(l):
        return -np.exp(-l + alpha**2 * l**2 / Y)

                                                                    
    n_chk = int(round((args.lam_max - args.lam_min) / args.check_step)) + 1
    lam_chk = np.linspace(args.lam_min, args.lam_max, n_chk)
    D_chk = np.array([defect(l) for l in lam_chk])
    disc = np.abs(D_chk - L_finite(lam_chk))
    i_max = int(np.argmax(disc))

                
    n_pl = int(round((args.lam_max - args.lam_min) / args.plot_step)) + 1
    lam_pl = np.linspace(args.lam_min, args.lam_max, n_pl)
    D_pl = np.array([D_chk[np.argmin(np.abs(lam_chk - l))]
                     if np.min(np.abs(lam_chk - l)) < 1e-12 else defect(l)
                     for l in lam_pl])
    lam_fine = np.linspace(args.lam_min, args.lam_max, 400)

    fig, ax = plt.subplots(figsize=(5.6, 3.6))
                                                                             
    ax.plot(lam_fine, -np.exp(-lam_fine), color="0.35", ls="--", lw=1.0,
            zorder=3, label=r"limiting profile $-e^{-\lambda}$")
    ax.plot(lam_fine, L_finite(lam_fine), color="#1f4e8c", lw=1.4,
            label=r"selected-zero prediction $\mathcal{L}_{\rho_0,X,\alpha}(s_\lambda)$")
    ax.plot(lam_pl, D_pl, ls="none", marker="o", ms=3.6, mfc="none",
            mec="#b2182b", mew=1.0,
            label=r"prime-side defect $\Delta_{X,\alpha}(s_\lambda)$")
    ax.set_xlabel(r"$\lambda$  ($s_\lambda=\rho_0+\lambda/\log X$)")
    ax.set_ylabel("defect")
    ax.set_xlim(args.lam_min - 0.05, args.lam_max + 0.05)
    ax.legend(loc="lower right", frameon=False)
    k = math.log10(X)
    Xs = rf"10^{ {int(round(k))}} " if abs(k - round(k)) < 1e-12 else f"{X:g}"
    Ns = f"{N:,}".replace(",", "{,}")
    ax.text(0.03, 0.95, rf"$X={Xs}$, $\alpha={alpha:g}$, $N={Ns}$",
            transform=ax.transAxes, va="top", fontsize=9)

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    pdf = out / "prime_side_first_zero_overlay.pdf"
    png = out / "prime_side_first_zero_overlay.png"
    fig.savefig(pdf)
    fig.savefig(png)
    plt.close(fig)

    meta = {
        "script": "fig1_prime_side_overlay.py",
        "rho0": f"0.5 + {GAMMA_1} i",
        "X": X, "alpha": alpha, "N": N, "Y": Y,
        "cutoff_argument": cutoff_arg,
        "prime_powers_used": int(logn.size),
        "lambda_range": [args.lam_min, args.lam_max],
        "check_step": args.check_step,
        "max_discrepancy": float(disc[i_max]),
        "lambda_at_max": float(lam_chk[i_max]),
        "plotted": {"lambda": lam_pl.tolist(), "Delta": D_pl.tolist(),
                    "L": L_finite(lam_pl).tolist()},
    }
    (out / "prime_side_first_zero_overlay.json").write_text(json.dumps(meta, indent=2))

    print(f"prime powers <= N: {logn.size}; cutoff argument = {cutoff_arg:.2f}")
    print(f"max |Delta - L| over lambda in [{args.lam_min}, {args.lam_max}] "
          f"(step {args.check_step}) = {disc[i_max]:.3e} at lambda = {lam_chk[i_max]:.3f}")
    print(pdf)
    print(png)


if __name__ == "__main__":
    main()
