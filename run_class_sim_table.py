#!/usr/bin/env python
"""Run CLASS on all ini files in a directory.

Equivalent to run_class.sh but with configurable paths for the CLASS
executable, ini file directory, and output directory.
"""

import argparse
import glob
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Run CLASS on all ini files in a directory.")
    parser.add_argument("--class-exe", type=str,
                        default=os.environ.get("CLASS_EXE",
                                               "../class_dmeff_rui_used/class"),
                        help="path to CLASS executable (default: CLASS_EXE env "
                             "var or ../class_dmeff_rui_used/class)")
    parser.add_argument("--ini-dir", type=str, default="inis",
                        help="directory containing ini files (default: inis/)")
    parser.add_argument("--output-dir", type=str, default="output",
                        help="directory for CLASS output and logs (default: output/)")
    args = parser.parse_args()

    if not os.path.isfile(args.class_exe):
        print(f"ERROR: CLASS executable not found at {args.class_exe}")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    ini_files = sorted(glob.glob(os.path.join(args.ini_dir, "*.ini")))
    total = len(ini_files)

    if total == 0:
        print(f"No ini files found in {args.ini_dir}/")
        sys.exit(1)

    print(f"Running CLASS for {total} ini files...")
    print(f"Executable: {args.class_exe}")
    print(f"Ini dir:    {args.ini_dir}/")
    print(f"Output dir: {args.output_dir}/")
    print("===========================================")

    succeeded = 0
    failed = 0

    for i, ini in enumerate(ini_files):
        name = os.path.basename(ini)
        log_path = os.path.join(args.output_dir,
                                os.path.splitext(name)[0] + ".log")
        print(f"[{i+1}/{total}] {args.class_exe} {ini}")

        with open(log_path, "w") as log_file:
            result = subprocess.run([args.class_exe, ini],
                                    stdout=log_file, stderr=subprocess.STDOUT)
        if result.returncode == 0:
            succeeded += 1
        else:
            failed += 1
            print(f"  *** FAILED (see {log_path})")

    print("===========================================")
    print(f"Done: {succeeded} succeeded, {failed} failed out of {total}")


if __name__ == "__main__":
    main()
