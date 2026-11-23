#!/usr/bin/env bash
# Step 07 — synteny (foreground vs HD100) and genomic-island staging/summary.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
python3 "$SCRIPTS/07_synteny.py" --out "$OUT" --query "$FOREGROUND" --ref "$REFERENCE_ID" --threads "$THREADS"

mkdir -p "$OUT/islands/upload" "$OUT/islands/results"
for s in $FOREGROUND $REFERENCE_ID; do
  cp "$OUT/gbk/$s.gbk" "$OUT/islands/upload/$s.gbk"
  n=$(grep -c '^LOCUS' "$OUT/gbk/$s.gbk")
  [[ "$n" -gt 1 ]] && echo "NOTE: $s has $n contigs — IslandViewer needs a complete genome or reference-guided contig ordering."
done
echo "Upload $OUT/islands/upload/*.gbk to IslandViewer 4, save results as $OUT/islands/results/<id>.tsv"
python3 "$SCRIPTS/07_islands.py" --out "$OUT"
