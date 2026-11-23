#!/usr/bin/env bash
# Step 00 — collect genomes, annotate with Prokka (uniform gene calls for all 41 genomes),
# and stage .faa / .ffn / .gff / .gbk / .fna files under $OUT/annot.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
mkdir -p "$OUT"/{annot,faa,ffn,gff,gbk,fna,logs}

shopt -s nullglob
files=("$GENOME_DIR"/*.fna "$GENOME_DIR"/*.fa "$GENOME_DIR"/*.fasta "$GENOME_DIR"/*.fas \
       "$GENOME_DIR"/*.fna.gz "$GENOME_DIR"/*.fa.gz "$GENOME_DIR"/*.fasta.gz)
[[ ${#files[@]} -eq 0 ]] && { echo "No genome FASTA files found in $GENOME_DIR"; exit 1; }
echo "Found ${#files[@]} genome files."

: > "$OUT/genome_ids.txt"
: > "$OUT/locus_tags.tsv"
declare -A seen_id seen_tag
for f in "${files[@]}"; do
  base=$(basename "$f"); id=${base%.gz}; id=${id%.*}
  if [[ -n "${seen_id[$id]:-}" ]]; then
    echo "ERROR: two files give the same genome ID '$id' (e.g. .fna and .fasta). Keep only one."; exit 1
  fi
  seen_id[$id]=1
  echo "$id" >> "$OUT/genome_ids.txt"
  # unique, Prokka-compliant locus tag from the strain name (last '_' field), e.g. kdesi, HD100, HCB162
  tag=$(echo "${id##*_}" | tr -cd 'A-Za-z0-9' | cut -c1-10); [[ -z "$tag" ]] && tag=G
  t=$tag; k=2; while [[ -n "${seen_tag[$t]:-}" ]]; do t="${tag:0:8}$k"; k=$((k+1)); done
  seen_tag[$t]=1; printf "%s\t%s\n" "$id" "$t" >> "$OUT/locus_tags.tsv"
  if [[ -s "$OUT/faa/$id.faa" ]]; then echo "  [skip] $id already annotated"; continue; fi
  tmp="$OUT/annot/$id.input.fna"
  if [[ "$f" == *.gz ]]; then gunzip -c "$f" > "$tmp"; else cp "$f" "$tmp"; fi
  locus=$t
  echo "  [prokka] $id  (locus tag $locus)"
  $PROKKA_CMD --kingdom Bacteria --genus Bdellovibrio --usegenus --gcode 11 \
         --locustag "$locus" --prefix "$id" --outdir "$OUT/annot/$id" \
         --cpus "$THREADS" --compliant --centre X --force --quiet "$tmp" \
         > "$OUT/logs/prokka_$id.log" 2>&1
  for ext in faa ffn gff gbk fna; do cp "$OUT/annot/$id/$id.$ext" "$OUT/$ext/$id.$ext"; done
  rm -f "$tmp"
done

n=$(wc -l < "$OUT/genome_ids.txt")
echo "Genome IDs ($n):"; tr '\n' ' ' < "$OUT/genome_ids.txt"; echo
for s in $FOREGROUND $BACKGROUND $REFERENCE_ID $SECRETOME_BENCHMARK; do
  grep -qx "$s" "$OUT/genome_ids.txt" || echo "WARNING: strain ID '$s' in config.sh not found among genome IDs"
done
