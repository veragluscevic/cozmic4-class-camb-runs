#!/usr/bin/env python
"""
Plot (d_dmeff / d_cdm)^2 for a given (n, m) parameter set.

Reads sigma values from sim-table.dat for halfmode, midpoint, and envelope
types, loads the corresponding transfer function files from output/, and
plots the squared ratio against the CDM reference, alongside a WDM
transfer function for the same particle mass.

Usage:
    python plot_transfer.py -n 2 -m 1e-1
"""

import argparse
import os
import math
import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SIM_TABLE = os.path.join(SCRIPT_DIR, "sim-table.dat")
CDM_FILE = os.path.join(SCRIPT_DIR, "CDM_class_sync_tk.dat")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")


def format_sci(val):
    """Format a number in scientific notation matching the filename convention.

    Strips trailing zeros from the mantissa, e.g. 4.2000e-28 -> '4.2e-28',
    9.0e-23 -> '9e-23', 2.6145e-24 -> '2.6145e-24'.
    """
    exp = int(math.floor(math.log10(abs(val))))
    mantissa = round(val / 10**exp, 4)
    return f"{mantissa:g}e{exp}"


def parse_sim_table(filepath, n_target, m_target):
    """Return list of (type, sigma) entries from sim-table.dat for given (n, m)."""
    entries = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            n_val = int(parts[0])
            m_val = float(parts[1])
            sigma_str = parts[2]
            stype = parts[3]
            if n_val == n_target and math.isclose(m_val, m_target, rel_tol=1e-6):
                if sigma_str.lower() == "nan":
                    continue
                entries.append({"type": stype, "sigma": float(sigma_str)})
    return entries


def build_filename(n, m, sigma):
    return f"n{n}_{format_sci(m)}GeV_{format_sci(sigma)}_sync_tk.dat"


def load_columns(filepath):
    """Load k (col 0), d_cdm (col 3), d_dmeff (col 4) from a transfer file."""
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 3], data[:, 4]


H = 0.7
OMEGA_MH2 = 0.11711
WDM_A = 0.0437
WDM_B = -1.188
WDM_NU = 1.049
WDM_THETA = 2.012
WDM_ETA = 0.2463


def transfer_wdm(k, mwdm):
    alpha = WDM_A * (mwdm**WDM_B) * ((OMEGA_MH2 / 0.12)**WDM_ETA) * ((H / 0.6736)**WDM_THETA) / H
    return (1 + (alpha * k) ** (2 * WDM_NU)) ** (-5.0 / WDM_NU)


def main():
    parser = argparse.ArgumentParser(
        description="Plot (d_dmeff / d_cdm)^2 for a given (n, m) parameter set."
    )
    parser.add_argument("-n", type=int, default=2, help="power-law index n (default: 2)")
    parser.add_argument("-m", type=float, default=1e-2, help="DM mass in GeV (default: 1e-2)")
    parser.add_argument("-s", "--save", action="store_true", help="save figure to plots/")
    args = parser.parse_args()

    n = args.n
    m = args.m

    k_cdm, d_cdm, _ = load_columns(CDM_FILE)

    entries = parse_sim_table(SIM_TABLE, n, m)
    if not entries:
        print(f"No entries found in sim-table.dat for n={n}, m={m}")
        raise SystemExit(1)

    LW = 2.2

    styles = {
        "halfmode":  {"color": "C0", "ls": "-"},
        "midpoint":  {"color": "C1", "ls": "--"},
        "envelope":  {"color": "C2", "ls": "-."},
    }

    fig, ax = plt.subplots(figsize=(8, 5))
    plotted = False

    for stype in ["halfmode", "midpoint", "envelope"]:
        entry = next((e for e in entries if e["type"] == stype), None)
        if entry is None:
            continue

        fname = build_filename(n, m, entry["sigma"])
        fpath = os.path.join(OUTPUT_DIR, fname)
        if not os.path.isfile(fpath):
            print(f"  File not found, skipping: {fname}")
            continue

        k, _, d_dmeff = load_columns(fpath)
        ratio_sq = (d_dmeff / d_cdm) ** 2

        st = styles[stype]
        label = rf"{stype} ($\sigma = {format_sci(entry['sigma'])}$)"
        ax.semilogx(k, ratio_sq, color=st["color"], ls=st["ls"], lw=LW, label=label)
        plotted = True

    if not plotted:
        print("No transfer function files found for the given parameters.")
        raise SystemExit(1)

    MWDM_KEV = 5.9
    k_plot = np.logspace(np.log10(1), np.log10(k_cdm.max()), 500)
    ax.semilogx(k_plot, transfer_wdm(k_plot, MWDM_KEV)**2, color="black", ls="-", lw=LW,
                label=rf"WDM ($m_\mathrm{{wdm}} = {MWDM_KEV}$ keV)")

    m_label = format_sci(m)
    ax.set_xlabel(r"$k$ [h/Mpc]")
    ax.set_ylabel(r"$T^2 / T_{\mathrm{CDM}}^2$")
    ax.set_title(rf"$n = {n}$, $m = {m_label}$ GeV")
    ax.legend()
    ax.set_xlim(1, k_cdm.max())
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if args.save:
        plotdir = os.path.join(SCRIPT_DIR, "plots")
        os.makedirs(plotdir, exist_ok=True)
        outfig = os.path.join(plotdir, f"transfer_ratio_n{n}_{m_label}GeV.png")
        fig.savefig(outfig, dpi=150)
        print(f"Saved: {outfig}")
    plt.show()


if __name__ == "__main__":
    main()
