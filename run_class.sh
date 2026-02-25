#!/usr/bin/env bash
set -euo pipefail

CLASS_EXE="../class_dmeff_rui_used/class"
INI_DIR="inis"
OUT_DIR="output"

if [[ ! -x "$CLASS_EXE" ]]; then
    echo "ERROR: CLASS executable not found at $CLASS_EXE"
    exit 1
fi

mkdir -p "$OUT_DIR"

ini_files=("$INI_DIR"/*.ini)
total=${#ini_files[@]}

if [[ $total -eq 0 ]]; then
    echo "No ini files found in $INI_DIR/"
    exit 1
fi

echo "Running CLASS for $total ini files..."
echo "Executable: $CLASS_EXE"
echo "Output dir: $OUT_DIR/"
echo "==========================================="

failed=0
succeeded=0

for i in "${!ini_files[@]}"; do
    ini="${ini_files[$i]}"
    name=$(basename "$ini")
    count=$((i + 1))

    echo "[$count/$total] $CLASS_EXE $ini"

    if $CLASS_EXE "$ini" > "${OUT_DIR}/${name%.ini}.log" 2>&1; then
        succeeded=$((succeeded + 1))
    else
        failed=$((failed + 1))
        echo "  *** FAILED (see ${OUT_DIR}/${name%.ini}.log)"
    fi
done

echo "==========================================="
echo "Done: $succeeded succeeded, $failed failed out of $total"
