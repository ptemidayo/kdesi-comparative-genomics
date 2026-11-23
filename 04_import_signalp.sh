#!/usr/bin/env bash
# Put SignalP-6 results produced elsewhere into $OUT/signalp/<genome_id>/prediction_results.txt
# Usage:  bash scripts/04_import_signalp.sh <folder_with_results>
# The folder holds one prediction_results.txt per input file, renamed to match it:
#   <genome_id>.txt                          (whole proteome in one job), or
#   <genome_id>__part1.txt, __part2.txt ...  (if the proteome was split) — parts are merged.
# Every protein in $OUT/faa/<genome_id>.faa must be covered; the count is checked.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
src=${1:?give the folder with SignalP result files}
shopt -s nullglob
ok=0; bad=0
while read -r id; do
  parts=()
  for p in "$src/$id.txt" "$src/${id}__part"*.txt; do [[ -s "$p" ]] && parts+=("$p"); done
  if (( ${#parts[@]} == 0 )); then echo "  [missing] $id"; bad=$((bad+1)); continue; fi
  mkdir -p "$OUT/signalp/$id"
  cat "${parts[@]}" | tr -d '\r' > "$OUT/signalp/$id/prediction_results.txt"
  n_pred=$(grep -vc '^#' "$OUT/signalp/$id/prediction_results.txt" || true)
  n_prot=$(grep -c '^>' "$OUT/faa/$id.faa")
  if [[ "$n_pred" -eq "$n_prot" ]]; then echo "  [ok] $id  ($n_pred proteins)"; ok=$((ok+1))
  else
    mv "$OUT/signalp/$id/prediction_results.txt" "$OUT/signalp/$id/prediction_results.INCOMPLETE.txt"
    echo "  [CHECK] $id: $n_pred predictions but $n_prot proteins — a part is missing or the wrong .faa was used (not imported)"
    bad=$((bad+1))
  fi
done < "$OUT/genome_ids.txt"
echo "Imported $ok genome(s); $bad need attention."
