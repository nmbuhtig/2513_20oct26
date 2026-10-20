#!/usr/bin/env python3
   

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import cm
from matplotlib.colors import TwoSlopeNorm


                                                                            
       
                                                                            
def configure_matplotlib(usetex: bool = False) -> None:
                                                
    plt.rcParams.update({
        "text.usetex": usetex,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.labelsize": 10,
        "font.size": 10,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
        "axes.unicode_minus": False,
    })


def format_X(X: float) -> str:
                                                                             
    k = math.log10(X)
    if abs(k - round(k)) < 1e-12:
        return rf"10^{ {int(round(k))}} "
    return f"{X:g}"


def offset_tag(beta_offset: float) -> str:
                                                                               
    return f"{beta_offset:.2f}".replace("-", "m").replace(".", "p")


                                                                            
                  
                                                                            
def selected_zero_profile(
    a: np.ndarray,
    b: np.ndarray,
    Y: float,
    alpha: float,
    beta_offset: float,
    eps: float,
) -> np.ndarray:
           
    lam = a + 1j * b
    radius = np.hypot(a, b)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        Z = -(beta_offset + a / Y) * Y * np.real(
            np.exp(-lam + (alpha**2 * lam**2) / Y) / lam
        )
    return np.where(radius < eps, np.nan, Z)


def make_grid(
    a_min: float, a_max: float, b_min: float, b_max: float, grid_size: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    a_vals = np.linspace(a_min, a_max, grid_size)
    b_vals = np.linspace(b_min, b_max, grid_size)
    A, B = np.meshgrid(a_vals, b_vals)
    return a_vals, b_vals, A, B


                                                                            
                              
                                                                            
def finite_stats(Z: np.ndarray) -> Dict[str, float]:
    finite = Z[np.isfinite(Z)]
    if finite.size == 0:
        nan = float("nan")
        return {"min": nan, "max": nan, "mean": nan,
                "median": nan, "std": nan, "max_abs": nan}
    return {
        "min": float(np.min(finite)),
        "max": float(np.max(finite)),
        "mean": float(np.mean(finite)),
        "median": float(np.median(finite)),
        "std": float(np.std(finite)),
        "max_abs": float(np.max(np.abs(finite))),
    }


def robust_vmax(Z: np.ndarray, percentile: float) -> float:
                                                                           
    finite = Z[np.isfinite(Z)]
    if finite.size == 0:
        return 1.0
    vmax = float(np.percentile(np.abs(finite), percentile))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = float(np.max(np.abs(finite))) or 1.0
    return vmax


def make_norm(vmax: float) -> TwoSlopeNorm:
    return TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)


def clip_info(Z: np.ndarray, vmax: float) -> Dict[str, Any]:
                                                                                  
    finite = Z[np.isfinite(Z)]
    below = bool(np.any(finite < -vmax))
    above = bool(np.any(finite > vmax))
    extend = {(False, False): "neither", (True, False): "min",
              (False, True): "max", (True, True): "both"}[(below, above)]
    frac = float(np.mean(np.abs(finite) > vmax)) if finite.size else 0.0
    return {"color_vmin": -vmax, "color_vmax": vmax,
            "clipped_fraction": frac, "colorbar_extend": extend}


                                                                            
               
                                                                            
def downsample_grid(
    a_vals: np.ndarray, b_vals: np.ndarray, Z: np.ndarray, target: int = 51,
) -> Dict[str, Any]:
                                                               
    ia = np.unique(np.linspace(0, len(a_vals) - 1, min(target, len(a_vals))).astype(int))
    ib = np.unique(np.linspace(0, len(b_vals) - 1, min(target, len(b_vals))).astype(int))
    Z_small = Z[np.ix_(ib, ia)]
    return {
        "a_values": [float(a_vals[i]) for i in ia],
        "b_values": [float(b_vals[i]) for i in ib],
        "z_values": [[None if not np.isfinite(x) else float(x) for x in row]
                     for row in Z_small],
    }


def make_metadata(
    *,
    X: float,
    Y: float,
    alpha: float,
    beta_offset: float | None,
    domain: Dict[str, Any],
    rendering: Dict[str, Any],
    stats: Dict[str, Any],
    color: Dict[str, Any],
    outputs: Dict[str, str],
    downsampled_grid: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    return {
        "figure_type": "selected_zero_profile_surface",
        "script": "SM_5.py",
        "manuscript_caption_note": (
            "This figure visualizes the selected-zero local theorem profile, "
            "not the full prime-force defect. It omits the nonlocal remainder."
        ),
        "formula_latex": (
            r"\mathcal L_{\rho_0,X,\alpha}\!\left(\rho_0+\frac{\lambda}{Y}\right)"
            r"= -\left(\beta_0-\frac12+\frac{a}{Y}\right)Y"
            r"\operatorname{Re}\left("
            r"\frac{\exp(-\lambda+\alpha^2\lambda^2/Y)}{\lambda}"
            r"\right),\quad \lambda=a+ib."
        ),
        "parameters": {
            "X": float(X),
            "Y": float(Y),
            "alpha": float(alpha),
            "beta_offset": None if beta_offset is None else float(beta_offset),
            "domain": domain,
            "rendering": rendering,
        },
        "statistics": stats,
        "color_scale": color,
        "outputs": outputs,
        "downsampled_grid": downsampled_grid,
    }


def save_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


                                                                            
         
                                                                            
def set_axes_style(ax, title: str, elev: float, azim: float) -> None:
    ax.set_xlabel(r"$a=\operatorname{Re}\lambda$", labelpad=6)
    ax.set_ylabel(r"$b=\operatorname{Im}\lambda$", labelpad=6)
    ax.set_zlabel(r"$\mathcal{L}_{\rho_0,X,\alpha}$", labelpad=6)
    ax.set_title(title, pad=4)
    ax.view_init(elev=elev, azim=azim)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.zaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.xaxis.pane.set_alpha(0.03)
    ax.yaxis.pane.set_alpha(0.03)
    ax.zaxis.pane.set_alpha(0.03)
    ax.grid(True)


def draw_surface(ax, A, B, Z, *, norm: TwoSlopeNorm, mesh_stride: int):
                                                                                
    surf = ax.plot_surface(
        A, B, Z,
        cmap=cm.RdBu_r, norm=norm,
        linewidth=0.0, edgecolors="none", antialiased=True, alpha=0.88,
        rcount=min(A.shape[0], 160), ccount=min(A.shape[1], 160),
    )
    surf.set_rasterized(True)

    wire = ax.plot_wireframe(
        A[::mesh_stride, ::mesh_stride],
        B[::mesh_stride, ::mesh_stride],
        Z[::mesh_stride, ::mesh_stride],
        color="black", linewidth=0.33, alpha=0.55,
    )
    wire.set_rasterized(True)

                                                                    
    ax.contour(A, B, Z, levels=[0], colors="black", linewidths=1.0,
               offset=float(np.nanmin(Z)))
    return surf


def add_colorbar(fig, surf, cax, extend: str):
    cbar = fig.colorbar(surf, cax=cax, extend=extend)
                                                                             
    cbar.ax.set_title(r"$\mathcal{L}_{\rho_0,X,\alpha}$", fontsize=8, pad=8)
    cbar.ax.tick_params(labelsize=7)
    return cbar


def save_individual_figure(
    *, A, B, Z, X, Y, alpha, beta_offset, png_path, pdf_path,
    title, mesh_stride, elev, azim, percentile,
) -> Dict[str, Any]:
    vmax = robust_vmax(Z, percentile)
    info = clip_info(Z, vmax)

    fig = plt.figure(figsize=(6.75, 5.2))
                                                                        
                                
    ax = fig.add_axes([0.00, 0.02, 0.80, 0.90], projection="3d")
    cax = fig.add_axes([0.87, 0.22, 0.022, 0.52])

    surf = draw_surface(ax, A, B, Z, norm=make_norm(vmax), mesh_stride=mesh_stride)

    if abs(beta_offset) < 1e-14:
        profile = "Critical-line selected-zero profile"
    else:
        profile = rf"Off-line selected-zero profile, $\beta_0-1/2={beta_offset:g}$"
    figure_title = (
        profile + "\n" + rf"$\alpha={alpha:g}$, $X={format_X(X)}$, $Y\approx {Y:.3f}$"
        if title else ""
    )
    set_axes_style(ax, figure_title, elev, azim)
    add_colorbar(fig, surf, cax, info["colorbar_extend"])

    fig.savefig(png_path, format="png")
    fig.savefig(pdf_path, format="pdf")
    plt.close(fig)
    return info


def save_combined_figure(
    *, A, B, Z_critical, Z_offline, X, Y, alpha, beta_offset_offline,
    png_path, pdf_path, title, mesh_stride, elev, azim, percentile, color_scale,
) -> Dict[str, Any]:
    if color_scale == "shared":
        v = robust_vmax(np.concatenate([Z_critical.ravel(), Z_offline.ravel()]),
                        percentile)
        v1 = v2 = v
    else:
        v1 = robust_vmax(Z_critical, percentile)
        v2 = robust_vmax(Z_offline, percentile)
    info1, info2 = clip_info(Z_critical, v1), clip_info(Z_offline, v2)

    fig = plt.figure(figsize=(12.5, 5.2))
    top = 0.86 if title else 0.94
                                                                               
    ax1 = fig.add_axes([0.00, 0.02, 0.43, top], projection="3d")
    ax2 = fig.add_axes([0.50, 0.02, 0.43, top], projection="3d")
    cax2 = fig.add_axes([0.955, 0.22, 0.011, 0.52])
    cax1 = fig.add_axes([0.455, 0.22, 0.011, 0.52]) if color_scale != "shared" else None

    surf1 = draw_surface(ax1, A, B, Z_critical, norm=make_norm(v1), mesh_stride=mesh_stride)
    surf2 = draw_surface(ax2, A, B, Z_offline, norm=make_norm(v2), mesh_stride=mesh_stride)

    title1 = "Critical-line selected-zero profile" if title else ""
    title2 = (rf"Off-line selected-zero profile, $\beta_0-1/2={beta_offset_offline:g}$"
              if title else "")
    set_axes_style(ax1, title1, elev, azim)
    set_axes_style(ax2, title2, elev, azim)

    if title:
        fig.suptitle(
            rf"Selected-zero profiles for $\alpha={alpha:g}$, $X={format_X(X)}$, "
            rf"$Y\approx {Y:.3f}$",
            y=0.99,
        )

    if color_scale == "shared":
        both = np.concatenate([Z_critical.ravel(), Z_offline.ravel()])
        add_colorbar(fig, surf2, cax2, clip_info(both, v2)["colorbar_extend"])
    else:
        add_colorbar(fig, surf1, cax1, info1["colorbar_extend"])
        add_colorbar(fig, surf2, cax2, info2["colorbar_extend"])

    fig.savefig(png_path, format="png")
    fig.savefig(pdf_path, format="pdf")
    plt.close(fig)
    return {"mode": color_scale, "critical_line": info1, "offline": info2}


                                                                            
              
                                                                            
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate all selected-zero profile graphics: individual and combined."
    )
    p.add_argument("--X", type=float, default=1e6)
    p.add_argument("--alpha", type=float, default=1.0)
    p.add_argument("--offline-beta-offset", type=float, default=0.10)
    p.add_argument("--a-min", type=float, default=-2.75)
    p.add_argument("--a-max", type=float, default=2.75)
    p.add_argument("--b-min", type=float, default=-2.75)
    p.add_argument("--b-max", type=float, default=2.75)
    p.add_argument("--grid-size", type=int, default=190)
    p.add_argument("--eps", type=float, default=0.055)
    p.add_argument("--mesh-stride", type=int, default=10)
    p.add_argument("--elev", type=float, default=32.0)
    p.add_argument("--azim", type=float, default=-55.0)
    p.add_argument("--percentile", type=float, default=99.0)
    p.add_argument("--color-scale", choices=["per-panel", "shared"], default="per-panel",
                   help="per-panel: one colourbar per panel (default); "
                        "shared: one common normalization and colourbar.")
    p.add_argument("--output-dir", type=str, default="prime_force_figures_final")
    p.add_argument("--usetex", action="store_true")
    p.add_argument("--no-title", action="store_true")
    return p


def main() -> None:
    args = build_parser().parse_args()
    if args.X <= 1:
        raise ValueError("X must be greater than 1.")
    if args.alpha <= 0:
        raise ValueError("alpha must be positive.")
    if args.grid_size < 30:
        raise ValueError("grid-size should be at least 30.")
    if args.mesh_stride < 2:
        raise ValueError("mesh-stride should be at least 2.")
    if not 50.0 <= args.percentile <= 100.0:
        raise ValueError("percentile should lie in [50, 100].")

    configure_matplotlib(usetex=args.usetex)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    X = float(args.X)
    Y = float(np.log(X))
    alpha = float(args.alpha)
    beta_off = float(args.offline_beta_offset)
    title = not args.no_title

    a_vals, b_vals, A, B = make_grid(args.a_min, args.a_max, args.b_min, args.b_max,
                                     args.grid_size)
    Z_critical = selected_zero_profile(A, B, Y=Y, alpha=alpha, beta_offset=0.0, eps=args.eps)
    Z_offline = selected_zero_profile(A, B, Y=Y, alpha=alpha, beta_offset=beta_off, eps=args.eps)

    domain = {
        "a_min": float(args.a_min), "a_max": float(args.a_max),
        "b_min": float(args.b_min), "b_max": float(args.b_max),
        "grid_size": int(args.grid_size), "masked_radius_eps": float(args.eps),
    }
    rendering = {
        "mesh_stride": int(args.mesh_stride), "elev": float(args.elev),
        "azim": float(args.azim), "color_norm_percentile": float(args.percentile),
        "usetex": bool(args.usetex), "title": title,
        "surface_rasterized_in_pdf": True,
        "note": ("Colour ranges are zero-centred and clipped at the stated percentile "
                 "of |L|; each colourbar matches the normalization of its panel."),
    }
    common = dict(mesh_stride=args.mesh_stride, elev=args.elev, azim=args.azim,
                  percentile=args.percentile, title=title)

                                       
    crit = {ext: out_dir / f"critical_line_selected_zero_profile.{ext}"
            for ext in ("png", "pdf", "json")}
    info = save_individual_figure(A=A, B=B, Z=Z_critical, X=X, Y=Y, alpha=alpha,
                                  beta_offset=0.0, png_path=crit["png"],
                                  pdf_path=crit["pdf"], **common)
    save_json(crit["json"], make_metadata(
        X=X, Y=Y, alpha=alpha, beta_offset=0.0, domain=domain, rendering=rendering,
        stats=finite_stats(Z_critical), color=info,
        outputs={"png": crit["png"].name, "pdf": crit["pdf"].name},
        downsampled_grid=downsample_grid(a_vals, b_vals, Z_critical)))

                                                                  
    tag = offset_tag(beta_off)
    off = {ext: out_dir / f"offline_selected_zero_profile_beta_offset_{tag}.{ext}"
           for ext in ("png", "pdf", "json")}
    info = save_individual_figure(A=A, B=B, Z=Z_offline, X=X, Y=Y, alpha=alpha,
                                  beta_offset=beta_off, png_path=off["png"],
                                  pdf_path=off["pdf"], **common)
    save_json(off["json"], make_metadata(
        X=X, Y=Y, alpha=alpha, beta_offset=beta_off, domain=domain, rendering=rendering,
        stats=finite_stats(Z_offline), color=info,
        outputs={"png": off["png"].name, "pdf": off["pdf"].name},
        downsampled_grid=downsample_grid(a_vals, b_vals, Z_offline)))

                                              
    combo = {ext: out_dir / f"selected_zero_profiles_comparison.{ext}"
             for ext in ("png", "pdf", "json")}
    info = save_combined_figure(A=A, B=B, Z_critical=Z_critical, Z_offline=Z_offline,
                                X=X, Y=Y, alpha=alpha, beta_offset_offline=beta_off,
                                png_path=combo["png"], pdf_path=combo["pdf"],
                                color_scale=args.color_scale, **common)
    save_json(combo["json"], make_metadata(
        X=X, Y=Y, alpha=alpha, beta_offset=None, domain=domain, rendering=rendering,
        stats={"critical_line": finite_stats(Z_critical), "offline": finite_stats(Z_offline)},
        color=info,
        outputs={"png": combo["png"].name, "pdf": combo["pdf"].name},
        downsampled_grid=None))

    for group in (crit, off, combo):
        for p in group.values():
            print(p)


if __name__ == "__main__":
    main()
