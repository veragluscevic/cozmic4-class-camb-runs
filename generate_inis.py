#!/usr/bin/env python3
"""Generate CLASS ini files for all valid rows in sim-table.csv."""

import os
import csv


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
    with open('minimal_newtonian.ini', 'r') as f:
        newt_template = f.read()
    with open('minimal_syncronous.ini', 'r') as f:
        sync_template = f.read()

    os.makedirs('inis', exist_ok=True)

    count = 0
    with open('sim-table.csv', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            n = row['n'].strip()
            m = row['m'].strip()
            sigma = row['sigma'].strip()

            if sigma.lower() == 'nan':
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


if __name__ == '__main__':
    main()
