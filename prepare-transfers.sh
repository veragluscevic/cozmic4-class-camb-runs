#!/usr/bin/env bash
set -euo pipefail

TABLE="sim-table.dat"

# Collect unique (n, m) pairs where at least one row has status != "done"
pairs=$(awk '!/^#/ && NF>=5 && $5 != "done" {print $1, $2}' "$TABLE" | sort -u)

if [ -z "$pairs" ]; then
    echo "All entries are done, nothing to run."
    exit 0
fi

echo "$pairs" | while read -r n m; do
    echo "=== n=$n  m=$m ==="
    python analyze_pk.py -n "$n" -m "$m" --recalculate-sigma --plot
done
