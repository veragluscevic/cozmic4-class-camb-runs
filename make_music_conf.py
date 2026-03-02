#!/usr/bin/env python3
"""Generate a MUSIC configuration file for a zoom-in halo.

Reads halo-specific parameters (ref_offset, ref_extent, seeds) from a
zoom-key file and substitutes them into a template config, along with the
specified transfer-function and white-noise filenames.

Usage
-----
    python make_music_conf.py <transfer_file> <wnoise_file> <halo_name> \
        [--template TEMPLATE] [--keyfile KEYFILE]

Output is written to  configs/music_<halo_name>_<transfer_stem>.conf
"""

import argparse
import os
import re
import sys


def parse_keyfile(path):
    """Return a dict mapping halo name -> dict of zoom parameters."""
    halos = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 10:
                print(f"Warning: skipping malformed line in {path}: {line}",
                      file=sys.stderr)
                continue
            name = parts[0]
            halos[name] = {
                "ref_offset": f"{parts[1]}, {parts[2]}, {parts[3]}",
                "ref_extent": f"{parts[4]}, {parts[5]}, {parts[6]}",
                "seed11": parts[7],
                "seed12": parts[8],
                "seed13": parts[9],
            }
    return halos


KEY_SECTIONS = {
    "ref_offset": "setup",
    "ref_extent": "setup",
    "transfer_file": "cosmology",
    "seed[10]": "random",
    "seed[11]": "random",
    "seed[12]": "random",
    "seed[13]": "random",
    "filename": "output",
}


def replace_value(lines, key, new_value):
    """Replace the value for 'key'; insert into the correct section if missing."""
    pattern = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*).*$")
    for i, line in enumerate(lines):
        m = pattern.match(line)
        if m:
            lines[i] = m.group(1) + new_value
            return

    new_line = f"{key:16s}= {new_value}"
    section = KEY_SECTIONS.get(key)
    if section is None:
        print(f"Warning: key '{key}' not found and has no known section; "
              "appending to end of file", file=sys.stderr)
        lines.append(new_line)
        return

    section_header = re.compile(rf"^\[{re.escape(section)}\]\s*$")
    next_section = re.compile(r"^\[.+\]\s*$")
    in_section = False
    insert_at = None
    for i, line in enumerate(lines):
        if section_header.match(line):
            in_section = True
            continue
        if in_section:
            if next_section.match(line):
                insert_at = i
                break
            insert_at = i + 1

    if insert_at is None:
        lines.append(f"[{section}]")
        lines.append(new_line)
    else:
        lines.insert(insert_at, new_line)

    print(f"Note: key '{key}' not in template; inserted into [{section}]",
          file=sys.stderr)


def main():
    default_template = "template-music.conf"
    default_keyfile = "zoom-key.dat"

    parser = argparse.ArgumentParser(
        description="Generate a MUSIC config for a zoom-in halo.")
    parser.add_argument("transfer_file", help="Transfer-function filename")
    parser.add_argument("wnoise_file", help="White-noise filename (seed[10])")
    parser.add_argument("halo_name", help="Halo name as listed in the key file")
    parser.add_argument("--template", default=default_template,
                        help=f"Template config file (default: {default_template})")
    parser.add_argument("--keyfile", default=default_keyfile,
                        help=f"Zoom-key file (default: {default_keyfile})")
    parser.add_argument("--icdir", default="ic",
                        help="Output directory for ICs (default: ic)")
    args = parser.parse_args()

    halos = parse_keyfile(args.keyfile)
    if args.halo_name not in halos:
        sys.exit(f"Error: halo '{args.halo_name}' not found in {args.keyfile}. "
                 f"Available: {', '.join(halos.keys())}")

    halo = halos[args.halo_name]

    with open(args.template) as f:
        lines = [line.rstrip("\n") for line in f]

    transfer_stem = os.path.splitext(os.path.basename(args.transfer_file))[0]

    replace_value(lines, "ref_offset", halo["ref_offset"])
    replace_value(lines, "ref_extent", halo["ref_extent"])
    replace_value(lines, "transfer_file", args.transfer_file)
    replace_value(lines, "seed[10]", args.wnoise_file)
    replace_value(lines, "seed[11]", halo["seed11"])
    replace_value(lines, "seed[12]", halo["seed12"])
    replace_value(lines, "seed[13]", halo["seed13"])
    ic_path = os.path.join(args.icdir, f"ic_gadget_{args.halo_name}_{transfer_stem}")
    replace_value(lines, "filename", ic_path)
    out_name = f"music_{args.halo_name}_{transfer_stem}.conf"
    out_path = os.path.join("configs", out_name)
    os.makedirs("configs", exist_ok=True)

    with open(out_path, "w") as f:
        for line in lines:
            f.write(line + "\n")

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
