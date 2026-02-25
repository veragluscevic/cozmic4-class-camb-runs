#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="output"
CONV_DIR="converted"

mkdir -p "$CONV_DIR"

count=0
failed=0

for sync_tk in "$OUT_DIR"/*_sync_tk.dat; do
    base=$(basename "$sync_tk" _sync_tk.dat)
    newt_tk="$OUT_DIR/${base}_newt_tk.dat"
    bg="$OUT_DIR/${base}_sync_background.dat"
    out="$CONV_DIR/camb_${base}_tk.dat"

    if [[ ! -f "$newt_tk" ]]; then
        echo "SKIP $base: missing $newt_tk"
        failed=$((failed + 1))
        continue
    fi
    if [[ ! -f "$bg" ]]; then
        echo "SKIP $base: missing $bg"
        failed=$((failed + 1))
        continue
    fi

    echo "Converting $base ..."
    echo "  python class_to_camb.py $sync_tk $newt_tk $bg -o $out"
    python class_to_camb.py "$sync_tk" "$newt_tk" "$bg" -o "$out"
    count=$((count + 1))
done

echo "==========================================="
echo "Done: $count converted, $failed skipped"
