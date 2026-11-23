#!/usr/bin/env bash
# Step 01 — pairwise ANI (skani) and genomospecies assignment at 95% ANI.
# If METADATA (config.sh) is set, that curated table becomes $OUT/metadata.tsv and the skani
# clusters are only reported as a cross-check (ARI). Otherwise skani clusters are used.
# The 'group' column is used for the Mann–Whitney comparisons in step 03; by default it
# equals the genomospecies — edit it (e.g. habitat, host range) if you compare other groups.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
mkdir -p "$OUT/ani"
ls "$OUT"/fna/*.fna > "$OUT/ani/genome_list.txt"
skani triangle -l "$OUT/ani/genome_list.txt" -t "$THREADS" -E > "$OUT/ani/skani_edges.tsv"
python3 "$SCRIPTS/01_ani_groups.py" --edges "$OUT/ani/skani_edges.tsv" \
        --ids "$OUT/genome_ids.txt" --cutoff "$ANI_CUTOFF" \
        --matrix "$OUT/ani/ani_matrix.tsv" --metadata "$OUT/metadata.tsv" --curated "${METADATA:-}"
