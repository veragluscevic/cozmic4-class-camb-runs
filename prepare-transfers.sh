#!/usr/bin/env bash
set -euo pipefail

TABLE="sim-table.dat"
INI_DIR="inis"
NEWT_TEMPLATE="minimal_newtonian.ini"
SYNC_TEMPLATE="minimal_syncronous.ini"

mkdir -p "$INI_DIR"

# --- Generate missing ini files ---
clean_sigma() {
    # Strip trailing zeros from mantissa: 4.2000e-28 -> 4.2e-28, 8.0000e-22 -> 8e-22
    echo "$1" | sed -E 's/(\.[0-9]*[1-9])0+e/\1e/g; s/\.0+e/e/g'
}

format_mass() {
    if [ "$1" = "1" ]; then echo "1e0"; else echo "$1"; fi
}

generate_ini() {
    local template_file="$1" n="$2" mass_fmt="$3" sigma_fmt="$4" gauge="$5"
    local root="output/n${n}_${mass_fmt}GeV_${sigma_fmt}_${gauge}_"
    sed \
        -e "s|^root = .*|root = ${root}|" \
        -e "s|^m_dmeff = .*|m_dmeff = ${mass_fmt}|" \
        -e "s|^sigma_dmeff = .*|sigma_dmeff = ${sigma_fmt}|" \
        -e "s|^npow_dmeff = .*|npow_dmeff = ${n}|" \
        "$template_file"
}

ini_count=0
while read -r line; do
    [[ -z "$line" || "$line" == \#* ]] && continue
    read -r n m sigma stype status rest <<< "$line"
    [[ "$(echo "$sigma" | tr '[:upper:]' '[:lower:]')" == "nan" ]] && continue

    mass_fmt=$(format_mass "$m")
    sigma_fmt=$(clean_sigma "$sigma")

    for gauge in newt sync; do
        fname="n${n}_${mass_fmt}GeV_${sigma_fmt}_${gauge}.ini"
        if [ ! -f "$INI_DIR/$fname" ]; then
            template="$NEWT_TEMPLATE"
            [ "$gauge" = "sync" ] && template="$SYNC_TEMPLATE"
            generate_ini "$template" "$n" "$mass_fmt" "$sigma_fmt" "$gauge" \
                > "$INI_DIR/$fname"
            ini_count=$((ini_count + 1))
        fi
    done
done < "$TABLE"

if [ "$ini_count" -gt 0 ]; then
    echo "Generated $ini_count missing ini files in $INI_DIR/"
fi

# --- Run analyze for (n, m) pairs that are not all "done" ---
pairs=$(awk '!/^#/ && NF>=5 && $5 != "done" {print $1, $2}' "$TABLE" | sort -u)

if [ -z "$pairs" ]; then
    echo "All entries are done, nothing to recalculate."
    exit 0
fi

echo "$pairs" | while read -r n m; do
    echo "=== n=$n  m=$m ==="
    python analyze_pk.py -n "$n" -m "$m" --recalculate-sigma --plot
done
