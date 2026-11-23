#!/usr/bin/env bash
# Step 01 — Pfam-A (release 36.0) domain annotation with hmmscan (HMMER v3.3.2).
set -euo pipefail
source "$(dirname "$0")/../config.sh"
mkdir -p "$PFAM_DIR" "$OUT/pfam/domtbl"

hmmscan -h | grep -q "HMMER 3.3.2" || echo "NOTE: Methods specify HMMER 3.3.2; installed: $(hmmscan -h | grep '^# HMMER')"

if [[ ! -s "$PFAM_DIR/Pfam-A.hmm.h3m" ]]; then
  echo "Downloading Pfam-A $PFAM_RELEASE ..."
  wget -q -O "$PFAM_DIR/Pfam-A.hmm.gz" \
    "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam${PFAM_RELEASE}/Pfam-A.hmm.gz"
  gunzip -f "$PFAM_DIR/Pfam-A.hmm.gz"
  hmmpress -f "$PFAM_DIR/Pfam-A.hmm"
fi

# Check the Pfam database is complete (e.g. after a crash during hmmpress)
for ext in h3m h3i h3f h3p; do
  [[ -s "$PFAM_DIR/Pfam-A.hmm.$ext" ]] || { echo "Pfam database incomplete ($ext missing) — rebuilding"; \
    rm -f "$PFAM_DIR"/Pfam-A.hmm.h3?; hmmpress -f "$PFAM_DIR/Pfam-A.hmm"; break; }
done

# Run genomes in parallel: JOBS concurrent hmmscan processes x CPU_PER_JOB threads each
CPU_PER_JOB=$(( THREADS < 4 ? THREADS : 4 )); (( CPU_PER_JOB < 1 )) && CPU_PER_JOB=1
JOBS=$(( THREADS / CPU_PER_JOB )); (( JOBS < 1 )) && JOBS=1
rm -f "$OUT"/pfam/domtbl/*.tmp                      # half-written results from an interrupted run
done_n=$(ls "$OUT"/pfam/domtbl/*.domtbl 2>/dev/null | wc -l)
echo "hmmscan: $done_n of $(wc -l < "$OUT/genome_ids.txt") genomes already done; $JOBS job(s) x $CPU_PER_JOB thread(s)"
export PFAM_DIR OUT PFAM_EVALUE CPU_PER_JOB
xargs -a "$OUT/genome_ids.txt" -P "$JOBS" -I{} bash -c '
  id={}; out="$OUT/pfam/domtbl/$id.domtbl"
  [[ -s "$out" ]] && { echo "  [skip] $id"; exit 0; }
  echo "  [hmmscan] $id"
  hmmscan --cpu "$CPU_PER_JOB" --noali -E "$PFAM_EVALUE" --domE "$PFAM_EVALUE" \
          --domtblout "$out.tmp" "$PFAM_DIR/Pfam-A.hmm" "$OUT/faa/$id.faa" > /dev/null
  mv "$out.tmp" "$out"'
echo "hmmscan finished: $(ls "$OUT"/pfam/domtbl/*.domtbl | wc -l) domain tables."
