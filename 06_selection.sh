#!/usr/bin/env bash
# Step 06 — branch-site selection analysis on the foreground lineage
# (representatives -> OrthoFinder -> filters -> MAFFT -> codon back-translation -> tree
#  -> codeml M0 + saturation filter -> branch-site test -> enrichment).
set -euo pipefail
source "$(dirname "$0")/../config.sh"
mkdir -p "$OUT/selection" "$OUT/logs"

if [[ "$BACKGROUND" == "auto" ]]; then
  echo "Choosing one representative genome per genomospecies:"
  STRAINS=$(python3 "$SCRIPTS/06_pick_representatives.py" --out "$OUT" --foreground "$FOREGROUND" --metadata "${METADATA:-}" --exclude "${CODEML_EXCLUDE:-}")
else
  STRAINS="$FOREGROUND $BACKGROUND"
fi
echo "$STRAINS" | tr ' ' '\n' > "$OUT/selection/codeml_strains.txt"
n=$(wc -l < "$OUT/selection/codeml_strains.txt")
echo "codeml strain set: $n genomes (foreground $FOREGROUND)"

# A different strain set invalidates earlier results -> keep them aside, start clean
OF="$OUT/selection/orthofinder"
if [[ -f "$OUT/selection/.strainset" ]] && ! cmp -s "$OUT/selection/.strainset" "$OUT/selection/codeml_strains.txt"; then
  old="$OUT/selection_previous_$(date +%Y%m%d_%H%M)"
  echo "Strain set changed — moving earlier selection results to $old"
  mkdir -p "$old"
  for d in orthofinder alignments codeml_runs orthogroup_candidates.tsv codeml_branchsite_results.tsv \
           saturation_filter.tsv positive_selection_enrichment.tsv tree_M0.nwk tree_branchsite.nwk \
           codeml_aliases.tsv concat_sco.faa concat_sco.treefile; do
    if [[ -e "$OUT/selection/$d" ]]; then mv "$OUT/selection/$d" "$old/"; fi
  done
fi
cp "$OUT/selection/codeml_strains.txt" "$OUT/selection/.strainset"
mkdir -p "$OF/proteomes"

if ! ls "$OF"/proteomes/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv >/dev/null 2>&1; then
  for s in $STRAINS; do cp "$OUT/faa/$s.faa" "$OF/proteomes/$s.faa"; done
  echo "Running OrthoFinder on $n proteomes (log: $OUT/logs/orthofinder.log)"
  orthofinder -f "$OF/proteomes" -t "$THREADS" -a "$THREADS" > "$OUT/logs/orthofinder.log" 2>&1
fi

if [[ ! -s "$OUT/selection/orthogroup_candidates.tsv" ]]; then
  python3 "$SCRIPTS/06a_select_orthogroups.py" --out "$OUT" --strains "$STRAINS" \
     --orthofinder-dir "$OF/proteomes" --mode "$OG_MODE" --per-class "$OG_PER_CLASS" \
     --min-codons "$OG_MIN_CODONS" --max-gap "$OG_MAX_GAP" --min-seg "$OG_MIN_SEGSITES" \
     --threads "$THREADS" --seed "$RANDOM_SEED"
fi

python3 "$SCRIPTS/06b_tree.py" --out "$OUT" --strains "$STRAINS" --species-tree "$SPECIES_TREE" --threads "$THREADS"
python3 "$SCRIPTS/06c_codeml.py" --out "$OUT" --strains "$STRAINS" --threads "$THREADS" \
        --alpha "$FDR_ALPHA" --ds-max "$DS_MAX"
