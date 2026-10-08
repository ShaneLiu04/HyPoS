#!/bin/bash
# AR004 baseline capture (BEFORE any code change) — Release build, OMP=1.
# Usage: bash bench_ar004_baseline.sh <outdir>
set -u
OUT="${1:-/root/ar004-baseline}"
BIN=/root/build-release/hypos
mkdir -p "$OUT"

run_case() {
  local tag="$1"; shift
  local np="$1"; shift
  local runs=3
  for r in 1 2 3; do
    local dir="$OUT/${tag}_r${r}"
    mkdir -p "$dir"
    OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
      mpirun --oversubscribe -np "$np" "$BIN" \
      --nx 256 --ny 256 --omp-threads 1 --output-format json \
      --output-dir "$dir" "$@" > "$dir/stdout.log" 2>&1
  done
}

# Per-iteration timing baseline (tol 0 -> fixed 500 iterations)
run_case jacobi_np1 1 --solver jacobi --max-iter 500 --tol 0.0
run_case jacobi_np4 4 --solver jacobi --max-iter 500 --tol 0.0
run_case rbgs_np1_500 1 --solver red_black_gs --max-iter 500 --tol 0.0
run_case rbgs_np4_500 4 --solver red_black_gs --max-iter 500 --tol 0.0
run_case cg_np1_500 1 --solver cg --max-iter 500 --tol 0.0
run_case cg_np4_500 4 --solver cg --max-iter 500 --tol 0.0

# Convergence iteration-count baseline (64^2)
run_case jacobi64_conv 1 --nx 64 --ny 64 --solver jacobi --max-iter 200000 --tol 1e-6
run_case rbgs64_conv 1 --nx 64 --ny 64 --solver red_black_gs --max-iter 100000 --tol 1e-6
run_case cg64_conv 1 --nx 64 --ny 64 --solver cg --max-iter 30000 --tol 1e-7

# Collect medians (no python3 in this WSL — awk-based)
echo "case|iter_ms_med|iters_med|final_res_med"
for d in "$OUT"/*/; do
  tag=$(basename "$d")
  for r in 1 2 3; do
    f="$d/r$r/performance_report.json"
    [ -f "$f" ] || continue
    grep -o '"iter_time_ms"[^,}]*' "$f" | head -1
    grep -o '"iterations"[^,}]*' "$f" | head -1
    grep -o '"final_residual"[^,}]*' "$f" | head -1
  done | awk -v tag="$tag" '
    /iter_time_ms/ { gsub(/[^0-9.eE-]/,"",$2); ms[NR]=$2 }
    /"iterations"/ { it[NR]=$2 }
    /final_residual/ { gsub(/[^0-9.eE-]/,"",$2); fr[NR]=$2 }
    END {
      # simple median over the collected values
      n=0; for (k in ms) {n++; v[n]=ms[k]} asort(v); m=(n>0)?v[int((n+1)/2)]:0
      ni=0; for (k in it) {ni++; w[ni]=it[k]} asort(w); mi=(ni>0)?w[int((ni+1)/2)]:0
      nf=0; for (k in fr) {nf++; x[nf]=fr[k]} asort(x); mf=(nf>0)?x[int((nf+1)/2)]:0
      printf "%s|%.4f|%d|%.6e\n", tag, m, mi, mf
    }'
done
