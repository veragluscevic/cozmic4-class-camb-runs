#!/usr/bin/env python3
"""Generate CLASS ini files for all valid rows in sim-table.dat."""

import os
import argparse


def clean_sigma(sigma_str):
    """Strip trailing zeros from mantissa of scientific notation."""
    if 'e' in sigma_str or 'E' in sigma_str:
        sep = 'e' if 'e' in sigma_str else 'E'
        mantissa, exponent = sigma_str.split(sep)
        mantissa = mantissa.rstrip('0').rstrip('.')
        return f"{mantissa}e{exponent}"
    return sigma_str


def format_mass(mass_str):
    mass_str = mass_str.strip()
    if mass_str == '1':
        return '1e0'
    return mass_str


def generate_ini(template, n, mass_fmt, sigma_fmt, gauge_short):
    root = f"output/n{n}_{mass_fmt}GeV_{sigma_fmt}_{gauge_short}_"

    lines = template.split('\n')
    new_lines = []
    for line in lines:
        if line.startswith('root = '):
            new_lines.append(f'root = {root}')
        elif line.startswith('m_dmeff = '):
            new_lines.append(f'm_dmeff = {mass_fmt}')
        elif line.startswith('sigma_dmeff = '):
            new_lines.append(f'sigma_dmeff = {sigma_fmt}')
        elif line.startswith('npow_dmeff = '):
            new_lines.append(f'npow_dmeff = {n}')
        else:
            new_lines.append(line)

    return '\n'.join(new_lines)


def main():
    parser = argparse.ArgumentParser(
        description="Generate CLASS ini files from sim-table.dat")
    parser.add_argument(
        "--skip-done", action="store_true",
        help="skip rows whose status is 'done'")
    args = parser.parse_args()

    with open('minimal_newtonian.ini', 'r') as f:
        newt_template = f.read()
    with open('minimal_syncronous.ini', 'r') as f:
        sync_template = f.read()

    os.makedirs('inis', exist_ok=True)

    count = 0
    skipped = 0
    with open('sim-table.dat', 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue

            n, m, sigma, stype, status = (
                parts[0], parts[1], parts[2], parts[3], parts[4]
            )

            if sigma.lower() == 'nan':
                continue

            if args.skip_done and status == 'done':
                skipped += 1
                continue

            mass_fmt = format_mass(m)
            sigma_fmt = clean_sigma(sigma)

            for gauge_short, template in [('newt', newt_template),
                                          ('sync', sync_template)]:
                content = generate_ini(template, n, mass_fmt, sigma_fmt,
                                       gauge_short)
                filename = f"n{n}_{mass_fmt}GeV_{sigma_fmt}_{gauge_short}.ini"
                with open(os.path.join('inis', filename), 'w') as out:
                    out.write(content)
                count += 1

    print(f"Created {count} ini files in inis/")
    if skipped:
        print(f"Skipped {skipped} rows with status 'done'")


if __name__ == '__main__':
    main()
