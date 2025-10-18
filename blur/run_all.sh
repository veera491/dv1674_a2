#!/usr/bin/env bash
set -euo pipefail
R=${1:-15}
T=${2:-8}
mkdir -p data_o
for f in data/im{1,2,3,4}.ppm; do
  base=$(basename "$f" .ppm)
  ./blur_par "$R" "$f" "data_o/output_${base}_r${R}.ppm" "$T"
done
