#!/usr/bin/env bash
# Copy the Prokka protein files (from step 00) into a folder you can upload from.
# Usage:  bash scripts/04_export_for_signalp.sh <destination_folder> [max_sequences_per_file]
#   e.g.  bash scripts/04_export_for_signalp.sh "/mnt/c/Users/Sarah Elufisan/Desktop/signalp_upload"
#         bash scripts/04_export_for_signalp.sh "/mnt/c/.../signalp_upload" 1000   # split into parts
# Upload each file to SignalP-6.0 (organism: Other). Name each downloaded
# prediction_results.txt after its input file, e.g. <genome_id>.txt or <genome_id>__part2.txt.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
dest=${1:?give a destination folder}; chunk=${2:-0}
mkdir -p "$dest"
while read -r id; do
  f="$OUT/faa/$id.faa"
  [[ -s "$f" ]] || { echo "missing $f — run step 00 first"; exit 1; }
  if (( chunk > 0 )); then
    awk -v n="$chunk" -v pre="$dest/${id}__part" '
      /^>/ { if (c % n == 0) { p++; out = pre p ".faa" } c++ } { print > out }' "$f"
  else
    cp "$f" "$dest/$id.faa"
  fi
done < "$OUT/genome_ids.txt"
echo "Wrote $(ls "$dest"/*.faa | wc -l) protein files to $dest"
