#!/bin/bash
# AR006 halo communication benchmark: p2p (pack) vs datatype (derived
# datatypes, direct placement) vs collective (dist graph + alltoallw)
# (design 4.6 FP6; decision D6: record-style verdict, no pass/fail gate).
#
# Workload: Jacobi 256^2, --max-iter 200, tol 0.0, residual check pushed
# out of the loop (--residual-check-interval 10000). np=1 is the no-face
# anchor; np=4 is the 2x2 partition with real faces. OMP_NUM_THREADS=1
# per rank. Timing comes from the "halo_exchange" profiler region
# (default non-overlap path: begin+end inside the region; "halo_wait"
# only exists under --overlap-comm and is 0 here).
#
# Environment: WSL, Release build at /root/build-release.
set -u
EXE=/root/build-release/hypos
OUT=/root/bench-ar006
REPS=3
NX=256

export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

mkdir -p "$OUT"

# median of up to 3 values (mawk-safe: no asort)
med3() {
    echo "$1 $2 $3" | awk '{
        v1=$1; v2=$2; v3=$3;
        if (v1<=v2) { if (v2<=v3) print v2; else if (v1<=v3) print v3; else print v1 }
        else        { if (v1<=v3) print v1; else if (v2<=v3) print v3; else print v2 }
    }'
}

# run_case <np> <mode> -> echoes "np mode median_sec"
run_case() {
    local np=$1 mode=$2 i d sec secs
    secs=""
    for i in $(seq 1 $REPS); do
        d="$OUT/np${np}-${mode}-$i"
        rm -rf "$d"; mkdir -p "$d"
        OMP_NUM_THREADS=1 mpirun --oversubscribe -np "$np" \
            "$EXE" --nx "$NX" --ny "$NX" --max-iter 200 --tol 0.0 \
            --residual-check-interval 10000 --comm-mode "$mode" \
            --enable-profiling --output-dir "$d" > "$d/run.log" 2>&1
        sec=$(grep -A2 '"halo_exchange"' "$d/run.log" | grep -o '"total_sec": [0-9.eE+-]*' \
              | head -1 | sed 's/.*: //')
        if [ -z "$sec" ]; then echo "ERROR: no halo_exchange region in $d/run.log" >&2; exit 1; fi
        secs="$secs $sec"
    done
    local med
    med=$(med3 $secs)
    echo "$np $mode $med"
}

echo "== AR006 halo comm benchmark (Jacobi ${NX}^2, iter 200, OMP=1, median of $REPS) =="
echo "np comm_mode halo_exchange_sec"
for np in 1 4; do
    for mode in p2p datatype collective; do
        run_case "$np" "$mode"
    done
done
echo "== done =="
