#!/bin/bash
# AR009 pipelined CG benchmark: pcg vs cg on the manufactured-solution
# workload (design FP7 / B1 acceptance — record-style verdict; iteration
# counts and wall times enter PERFORMANCE §16, incl. the R1① 512^2
# half-item compensation).
#
# Workload: 256^2 and 512^2, tol 1e-6 (absolute true-residual L2, same
# convention as §14/§15), np = 1 / 4 (2x2), OMP_NUM_THREADS=1 per rank,
# median of 3 runs, WSL Release build at /root/build-release.
#
# Log-line formats parsed:
#   pcg: "PipelinedCG finished in N iterations, time = T s"
#   cg:  "CG finished in N iterations, time = T s"
set -u
EXE=/root/build-release/hypos
OUT=/root/bench-ar009
REPS=3

export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

mkdir -p "$OUT"

med3() {
    echo "$1 $2 $3" | awk '{
        v1=$1; v2=$2; v3=$3;
        if (v1<=v2) { if (v2<=v3) print v2; else if (v1<=v3) print v3; else print v1 }
        else        { if (v1<=v3) print v1; else if (v2<=v3) print v3; else print v2 }
    }'
}

# run_case <np> <nx> <solver> -> echoes "np nx solver iters med_sec"
run_case() {
    local np=$1 nx=$2 solver=$3 i d line it sec secs
    local maxiter=5000
    secs=""
    it=""
    for i in $(seq 1 $REPS); do
        d="$OUT/np${np}-n${nx}-${solver}-$i"
        rm -rf "$d"; mkdir -p "$d"
        OMP_NUM_THREADS=1 mpirun --oversubscribe -np "$np" \
            "$EXE" --nx "$nx" --ny "$nx" --max-iter "$maxiter" --tol 1e-6 \
            --solver "$solver" --output-dir "$d" > "$d/run.log" 2>&1
        line=$(grep -E 'finished in' "$d/run.log" | tail -1)
        if [ -z "$line" ]; then echo "ERROR: no finish line in $d/run.log" >&2; exit 1; fi
        sec=$(echo "$line" | grep -oE 'time = [0-9.eE+-]+' | sed 's/time = //')
        if [ -z "$sec" ]; then echo "ERROR: no time in: $line" >&2; exit 1; fi
        it=$(echo "$line" | grep -oE 'finished in [0-9]+ iterations' | grep -oE '[0-9]+')
        if [ -z "$it" ]; then echo "ERROR: parse: $line" >&2; exit 1; fi
        secs="$secs $sec"
    done
    local med
    med=$(med3 $secs)
    echo "$np $nx $solver $it $med"
}

echo "== AR009 pcg vs cg benchmark (manufactured solution, tol 1e-6, OMP=1, median of $REPS) =="
echo "np nx solver iterations wall_sec"
for nx in 256 512; do
    for np in 1 4; do
        for solver in cg pcg; do
            run_case "$np" "$nx" "$solver"
        done
    done
done
echo "== done =="
