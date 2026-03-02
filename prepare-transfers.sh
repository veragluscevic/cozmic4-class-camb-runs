#!/usr/bin/env bash
set -euo pipefail

# --- Parse arguments ---
INI_DIR="inis"
CLASS_TK_DIR="output"
CONV_DIR="transfers"
DO_ANALYZE=false

usage() {
    echo "Usage: $0 [--ini-dir DIR] [--class-tk-dir DIR] [--output-dir DIR] [--analyze]"
    echo "  --ini-dir      directory for ini files (default: inis/)"
    echo "  --class-tk-dir directory for CLASS transfer output (default: output/)"
    echo "  --output-dir   directory for CAMB-converted transfers (default: transfers/)"
    echo "  --analyze      also run analyze_pk.py for each (n, m) pair"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ini-dir)      INI_DIR="$2"; shift 2 ;;
        --class-tk-dir) CLASS_TK_DIR="$2"; shift 2 ;;
        --output-dir)   CONV_DIR="$2"; shift 2 ;;
        --analyze)      DO_ANALYZE=true; shift ;;
        -h|--help)      usage ;;
        *)              echo "Unknown option: $1"; usage ;;
    esac
done

TABLE="sim-table.dat"
NEWT_TEMPLATE="COZMIC1-files/minimal_newtonian.ini"
SYNC_TEMPLATE="COZMIC1-files/minimal_syncronous.ini"

mkdir -p "$INI_DIR"

# --- Helper functions ---
clean_sigma() {
    # Strip trailing zeros from mantissa: 4.2000e-28 -> 4.2e-28, 8.0000e-22 -> 8e-22
    echo "$1" | sed -E 's/(\.[0-9]*[1-9])0+e/\1e/g; s/\.0+e/e/g'
}

format_mass() {
    if [ "$1" = "1" ]; then echo "1e0"; else echo "$1"; fi
}

generate_ini() {
    local template_file="$1" n="$2" mass_fmt="$3" sigma_fmt="$4" gauge="$5"
    local root="${CLASS_TK_DIR}/n${n}_${mass_fmt}GeV_${sigma_fmt}_${gauge}_"
    sed \
        -e "s|^root = .*|root = ${root}|" \
        -e "s|^m_dmeff = .*|m_dmeff = ${mass_fmt}|" \
        -e "s|^sigma_dmeff = .*|sigma_dmeff = ${sigma_fmt}|" \
        -e "s|^npow_dmeff = .*|npow_dmeff = ${n}|" \
        "$template_file"
}

# --- Step 1: Generate missing ini files ---
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

# --- Step 2: Run CLASS for non-done entries with missing output ---
CLASS_EXE="${CLASS_EXE:-../class_dmeff_rui_used/class}"
mkdir -p "$CLASS_TK_DIR"

if [[ ! -x "$CLASS_EXE" ]]; then
    echo "ERROR: CLASS executable not found at $CLASS_EXE"
    exit 1
fi

class_ran=0
class_skip=0

while read -r line; do
    [[ -z "$line" || "$line" == \#* ]] && continue
    read -r n m sigma stype status rest <<< "$line"
    [[ "$status" == "done" ]] && continue
    [[ "$(echo "$sigma" | tr '[:upper:]' '[:lower:]')" == "nan" ]] && continue

    mass_fmt=$(format_mass "$m")
    sigma_fmt=$(clean_sigma "$sigma")
    base="n${n}_${mass_fmt}GeV_${sigma_fmt}"

    for gauge in sync newt; do
        ini="$INI_DIR/${base}_${gauge}.ini"
        if [[ ! -f "$ini" ]]; then
            echo "SKIP CLASS $base ($gauge): missing $ini"
            class_skip=$((class_skip + 1))
            continue
        fi

        if [ "$gauge" = "sync" ]; then
            out_check="$CLASS_TK_DIR/${base}_sync_tk.dat"
        else
            out_check="$CLASS_TK_DIR/${base}_newt_tk.dat"
        fi

        if [[ -f "$out_check" ]]; then
            continue
        fi

        echo "Running CLASS: $base ($gauge) ..."
        if $CLASS_EXE "$ini" > "$CLASS_TK_DIR/${base}_${gauge}.log" 2>&1; then
            class_ran=$((class_ran + 1))
        else
            class_skip=$((class_skip + 1))
            echo "  *** FAILED (see $CLASS_TK_DIR/${base}_${gauge}.log)"
        fi
    done
done < "$TABLE"

if [ "$class_ran" -gt 0 ] || [ "$class_skip" -gt 0 ]; then
    echo "==========================================="
    echo "CLASS: $class_ran ran, $class_skip skipped/failed"
fi

# --- Step 3 (optional): Run analyze for (n, m) pairs that are not all "done" ---
if [ "$DO_ANALYZE" = true ]; then
    pairs=$(awk '!/^#/ && NF>=5 && $5 != "done" {print $1, $2}' "$TABLE" | sort -u)

    if [ -z "$pairs" ]; then
        echo "All entries are done, nothing to analyze."
    else
        echo "$pairs" | while read -r n m; do
            echo "=== n=$n  m=$m ==="
            python analyze_pk.py -n "$n" -m "$m" --recalculate-sigma --plot
        done
    fi
fi

# --- Step 4: Convert non-done transfers to CAMB format ---
mkdir -p "$CONV_DIR"

conv_count=0
conv_skip=0

while read -r line; do
    [[ -z "$line" || "$line" == \#* ]] && continue
    read -r n m sigma stype status rest <<< "$line"
    [[ "$status" == "done" ]] && continue
    [[ "$(echo "$sigma" | tr '[:upper:]' '[:lower:]')" == "nan" ]] && continue

    mass_fmt=$(format_mass "$m")
    sigma_fmt=$(clean_sigma "$sigma")
    base="n${n}_${mass_fmt}GeV_${sigma_fmt}"

    sync_tk="$CLASS_TK_DIR/${base}_sync_tk.dat"
    newt_tk="$CLASS_TK_DIR/${base}_newt_tk.dat"
    bg="$CLASS_TK_DIR/${base}_sync_background.dat"
    out="$CONV_DIR/camb_${base}_tk.dat"

    if [[ ! -f "$sync_tk" ]]; then
        echo "SKIP convert $base: missing $sync_tk"
        conv_skip=$((conv_skip + 1))
        continue
    fi
    if [[ ! -f "$newt_tk" ]]; then
        echo "SKIP convert $base: missing $newt_tk"
        conv_skip=$((conv_skip + 1))
        continue
    fi
    if [[ ! -f "$bg" ]]; then
        echo "SKIP convert $base: missing $bg"
        conv_skip=$((conv_skip + 1))
        continue
    fi

    echo "Converting $base ..."
    python class_to_camb.py "$sync_tk" "$newt_tk" "$bg" -o "$out"
    conv_count=$((conv_count + 1))
done < "$TABLE"

echo "==========================================="
echo "Conversion: $conv_count converted, $conv_skip skipped"
