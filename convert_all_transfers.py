#!/usr/bin/env python
"""Convert all CLASS transfer outputs to CAMB 13-column format.

For each *_sync_tk.dat file in the input directory, looks for the
corresponding *_newt_tk.dat and *_sync_background.dat, then calls
class_to_camb.py to produce a CAMB-format file in the output directory.
"""

import argparse
import glob
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Convert all CLASS transfer outputs to CAMB format.")
    parser.add_argument("--class-tk-dir", type=str, default="output",
                        help="directory with CLASS output files (default: output/)")
    parser.add_argument("--output-dir", type=str, default="transfers",
                        help="directory for CAMB-format files (default: transfers/)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    converter = os.path.join(script_dir, "class_to_camb.py")

    sync_files = sorted(glob.glob(os.path.join(args.class_tk_dir, "*_sync_tk.dat")))

    if not sync_files:
        print(f"No *_sync_tk.dat files found in {args.class_tk_dir}/")
        sys.exit(1)

    converted = 0
    skipped = 0

    for sync_tk in sync_files:
        base = os.path.basename(sync_tk).replace("_sync_tk.dat", "")
        newt_tk = os.path.join(args.class_tk_dir, f"{base}_newt_tk.dat")
        bg = os.path.join(args.class_tk_dir, f"{base}_sync_background.dat")
        out = os.path.join(args.output_dir, f"camb_{base}_tk.dat")

        if not os.path.isfile(newt_tk):
            print(f"SKIP {base}: missing {newt_tk}")
            skipped += 1
            continue
        if not os.path.isfile(bg):
            print(f"SKIP {base}: missing {bg}")
            skipped += 1
            continue

        print(f"Converting {base} ...")
        result = subprocess.run(
            [sys.executable, converter, sync_tk, newt_tk, bg, "-o", out])
        if result.returncode == 0:
            converted += 1
        else:
            print(f"  *** FAILED for {base}")
            skipped += 1

    print("===========================================")
    print(f"Done: {converted} converted, {skipped} skipped")


if __name__ == "__main__":
    main()
