#!/usr/bin/env python
"""
Analyze transfer functions: compute k_hm for each (n, m, sigma) entry.

k_hm is the wavenumber where d_dmeff/d_cdm first drops to 0.5
(equivalently, where T^2/T_CDM^2 = 0.25).

Shares file conventions with plot_transfer.py.

Usage:
    python analyze_pk.py -n 2 -m 1e-2
    python analyze_pk.py -m all
    python analyze_pk.py -m all --save-khm
    python analyze_pk.py -n 2 -m 1e-2 --plot
"""

import argparse
import subprocess
import sys
import os
import re
import math
import numpy as np

from plot_transfer import (
    SCRIPT_DIR, SIM_TABLE, CDM_FILE, OUTPUT_DIR, H,
    format_sci, parse_sim_table, unique_masses,
    build_filename, load_columns,
)


def compute_khm(k, d_dmeff, d_cdm, threshold=0.5):
    """Find the smallest k where d_dmeff/d_cdm drops to `threshold`.

    Uses linear interpolation between adjacent grid points.
    Returns None if the ratio never drops below the threshold.
    """
    ratio = d_dmeff / d_cdm
    for i in range(len(ratio) - 1):
        if ratio[i] >= threshold and ratio[i + 1] < threshold:
            k_hm = k[i] + (threshold - ratio[i]) / (ratio[i + 1] - ratio[i]) * (k[i + 1] - k[i])
            return float(k_hm)
    return None


def save_khm_to_table(results):
    """Append or update k_hm column in sim-table.dat.

    Reads the file, adds/updates a k_hm column for matching rows,
    and writes it back.
    """
    with open(SIM_TABLE) as f:
        lines = f.readlines()

    def find_khm(n_val, m_val, sigma_str, stype):
        """Look up k_hm from results, matching sigma by float value."""
        sigma_file = sigma_str.lower()
        if sigma_file == "nan":
            return None
        sigma_fval = float(sigma_file)
        for r in results:
            if (r["n"] == int(n_val)
                    and math.isclose(r["m"], m_val, rel_tol=1e-6)
                    and math.isclose(r["sigma"], sigma_fval, rel_tol=1e-6)
                    and r["type"] == stype):
                return r["khm"]
        return None

    header_re = re.compile(r"^#\s*n\s+m\s+sigma\s+type\s+status")
    new_lines = []
    for line in lines:
        stripped = line.rstrip("\n")
        if header_re.match(stripped):
            if "khm" not in stripped:
                stripped = stripped + "  khm[1/Mpc]"
            new_lines.append(stripped + "\n")
            continue

        if not stripped.strip() or stripped.strip().startswith("#"):
            new_lines.append(stripped + "\n")
            continue

        parts = stripped.split()
        if len(parts) < 5:
            new_lines.append(stripped + "\n")
            continue

        n_str, m_str, sigma_str, stype, status = parts[0], parts[1], parts[2], parts[3], parts[4]
        matched = find_khm(n_str, float(m_str), sigma_str, stype)
        if matched is not None:
            khm_str = f"{matched * H:.2f}"
        elif sigma_str.lower() == "nan":
            khm_str = "nan"
        elif len(parts) >= 6:
            khm_str = parts[5]
        else:
            khm_str = "nan"

        new_lines.append(f"{n_str:>3}  {m_str:>4}  {sigma_str:>12}  {stype:<8}  {status:<4}  {khm_str}\n")

    with open(SIM_TABLE, "w") as f:
        f.writelines(new_lines)
    print(f"Updated {SIM_TABLE}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze transfer functions: compute k_hm."
    )
    parser.add_argument("-n", type=int, default=2, help="power-law index n (default: 2)")
    parser.add_argument("-m", type=str, default="1e-2",
                        help="DM mass in GeV, or 'all' (default: 1e-2)")
    parser.add_argument("--mwdm", type=float, default=5.9, help="WDM mass in keV (default: 5.9)")
    parser.add_argument("--no-midpoint", action="store_true", help="skip midpoint entries")
    parser.add_argument("-s", "--save", action="store_true", help="save plot figure to plots/")
    parser.add_argument("--save-khm", action="store_true", help="write k_hm values into sim-table.dat")
    parser.add_argument("--plot", action="store_true", help="call plot_transfer.py with the same arguments")
    args = parser.parse_args()

    n = args.n
    show_all = args.m.lower() == "all"

    k_cdm, d_cdm, _ = load_columns(CDM_FILE)

    if show_all:
        all_entries = parse_sim_table(SIM_TABLE, n)
    else:
        m_val = float(args.m)
        all_entries = parse_sim_table(SIM_TABLE, n, m_val)
    masses = unique_masses(all_entries)

    if not all_entries:
        print(f"No entries found in sim-table.dat for n={n}")
        raise SystemExit(1)

    results = []

    print(f"{'n':>3}  {'m':>6}  {'type':<9}  {'sigma':>14}  {'k_hm [1/Mpc]':>14}")
    print("-" * 58)

    for m in masses:
        m_entries = [e for e in all_entries if math.isclose(e["m"], m, rel_tol=1e-6)]
        for stype in ["halfmode", "midpoint", "envelope"]:
            if stype == "midpoint" and args.no_midpoint:
                continue
            entry = next((e for e in m_entries if e["type"] == stype), None)
            if entry is None:
                continue

            fname = build_filename(n, m, entry["sigma"])
            fpath = os.path.join(OUTPUT_DIR, fname)
            if not os.path.isfile(fpath):
                print(f"{n:>3}  {format_sci(m):>6}  {stype:<9}  {format_sci(entry['sigma']):>14}  {'(no file)':>14}")
                continue

            k, _, d_dmeff = load_columns(fpath)
            khm = compute_khm(k, d_dmeff, d_cdm)

            khm_str = f"{khm * H:.2f}" if khm is not None else "n/a"
            print(f"{n:>3}  {format_sci(m):>6}  {stype:<9}  {format_sci(entry['sigma']):>14}  {khm_str:>14}")

            results.append({
                "n": n,
                "m": m,
                "sigma": entry["sigma"],
                "type": stype,
                "khm": khm,
            })

    if args.save_khm:
        save_khm_to_table(results)

    if args.plot:
        cmd = [sys.executable, os.path.join(SCRIPT_DIR, "plot_transfer.py"),
               "-n", str(n), "-m", args.m,
               "--mwdm", str(args.mwdm)]
        if args.no_midpoint:
            cmd.append("--no-midpoint")
        if args.save:
            cmd.append("-s")
        subprocess.run(cmd)


if __name__ == "__main__":
    main()
