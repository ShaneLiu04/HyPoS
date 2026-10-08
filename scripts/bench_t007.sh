#!/bin/bash
# T007 benchmark: 256^2, np=1/4, OMP=1, 3 runs each, mpibin vs binary.
# Metric: io_write total_sec from --enable-profiling JSON report (rank 0).
set -u
HYPOS=/root/build-release/hypos
OUT=/root/bench
mkdir -p "$OUT"
export OMP_NUM_THREADS=1
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

for fmt in mpibin binary; do
  for np in 1 4; do
    for run in 1 2 3; do
      d="$OUT/${fmt}_np${np}_r${run}"
      rm -rf "$d"; mkdir -p "$d"
      mpirun.openmpi -n "$np" --oversubscribe "$HYPOS" \
        --nx 256 --ny 256 --max-iter 100 --save-interval 10 \
        --enable-profiling --output-format "$fmt" --output-dir "$d" \
        > "$d/stdout.log" 2>&1
      t=$(python3 -c "import re,sys; s=open('$d/stdout.log').read(); m=re.search(r'\"io_write\": \{\s*\"total_sec\": ([0-9.]+)', s); print(m.group(1) if m else 'MISSING')" 2>/dev/null || \
          grep -A1 '"io_write"' "$d/stdout.log" | grep total_sec | head -1 | sed 's/[^0-9.]//g')
      echo "$fmt np=$np run=$run io_write_total_sec=$t"
    done
  done
done
