#!/usr/bin/env bash
# Step 04 — SignalP-6.0 signal peptide prediction + secretome benchmarking.
#
# Two ways to get SignalP results:
#  A) local SignalP-6 (DTU licence): if `signalp6` is installed, missing genomes are run here.
#  B) run SignalP-6 elsewhere (DTU web server or another computer) on the Prokka protein files
#     made by step 00 ($OUT/faa/<genome_id>.faa — use THESE files so protein IDs match), then put
#     each result at:  $OUT/signalp/<genome_id>/prediction_results.txt
#     `bash scripts/04_export_for_signalp.sh <folder>` copies the protein files out for you, and
#     `bash scripts/04_import_signalp.sh <folder>` puts downloaded results in the right place.
#
# NOTE: SignalP-6 offers only "eukarya" or "other"; Gram-negative bacteria are run as "other".
set -euo pipefail
source "$(dirname "$0")/../config.sh"
mkdir -p "$OUT/signalp"

missing=()
while read -r id; do
  [[ -s "$OUT/signalp/$id/prediction_results.txt" ]] || missing+=("$id")
done < "$OUT/genome_ids.txt"

if (( ${#missing[@]} > 0 )); then
  if command -v signalp6 >/dev/null; then
    for id in "${missing[@]}"; do
      echo "  [signalp6] $id"
      signalp6 --fastafile "$OUT/faa/$id.faa" --organism other --output_dir "$OUT/signalp/$id" \
               --format none --mode "$SIGNALP_MODE" --write_procs "$THREADS" > "$OUT/logs/signalp_$id.log" 2>&1
    done
  else
    echo "SignalP results are missing for ${#missing[@]} genome(s), and signalp6 is not installed here:"
    printf '   %s\n' "${missing[@]}"
    echo "Run SignalP-6 elsewhere (organism: Other) on $OUT/faa/<genome_id>.faa and place each result at"
    echo "   $OUT/signalp/<genome_id>/prediction_results.txt"
    echo "Helpers: scripts/04_export_for_signalp.sh and scripts/04_import_signalp.sh"
    exit 1
  fi
fi

python3 "$SCRIPTS/04_secretome.py" --out "$OUT" --foreground "$FOREGROUND" \
        --benchmark "$SECRETOME_BENCHMARK"
