#!/bin/bash
# AR007 two-level multigrid benchmark: mg2 vs red_black_gs vs jacobi on the
# manufactured-solution workload (design 4.6 FP6 / srs R1 acceptance 2:
# record-style verdict — mg2 must use strictly fewer total smoothing
# sweeps than plain RBGS iterations; actual speedup is recorded as-is).
#
# Workload: 256^2, tol 1e-2 (both baselines converge within max-iter at this
# tolerance; at tol 1e-3 red_black_gs exceeds 100k iterations — recorded as
# "did not converge" in PERFORMANCE §14), np = 1 / 4 (2x2), OMP_NUM_THREADS=1
# per rank, median of 3 runs. Iteration counts and wall times come from the
# solver "finished in ..." log line; mg2 smoothing sweeps = 4 per cycle
# (nu1 = nu2 = 2).
#
# Environment: WSL, Release build at /root/build-release.
set -u
EXE=/root/build-release/hypos
OUT=/root/bench-ar007
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

# run_case <np> <solver> -> echoes "np solver iters sweeps med_sec"
# CLI solver names: jacobi / red_black_gs / mg2 (see --help).
run_case() {
    local np=$1 solver=$2 i d line it sweeps sec secs
    local maxiter=100000
    case "$solver" in
        mg2)          maxiter=2000 ;;
        red_black_gs) maxiter=200000 ;;
        jacobi)       maxiter=400000 ;;
    esac
    secs=""
    for i in $(seq 1 $REPS); do
        d="$OUT/np${np}-${solver}-$i"
        rm -rf "$d"; mkdir -p "$d"
        OMP_NUM_THREADS=1 mpirun --oversubscribe -np "$np" \
            "$EXE" --nx "$NX" --ny "$NX" --max-iter "$maxiter" --tol 1e-2 \
            --solver "$solver" --output-dir "$d" > "$d/run.log" 2>&1
        line=$(grep -E 'finished in' "$d/run.log" | tail -1)
        if [ -z "$line" ]; then echo "ERROR: no finish line in $d/run.log" >&2; exit 1; fi
        sec=$(echo "$line" | grep -oE 'time = [0-9.eE+-]+ s' | sed 's/time = //; s/ s//')
        if [ -z "$sec" ]; then
            sec=$(echo "$line" | grep -oE 'time = [0-9.eE+-]+' | sed 's/time = //')
        fi
        if [ -z "$sec" ]; then echo "ERROR: no time in finish line: $line" >&2; exit 1; fi
        if [ "$solver" = mg2 ]; then
            it=$(echo "$line" | grep -oE 'finished in [0-9]+ cycles' | grep -oE '[0-9]+')
            sweeps=$(echo "$line" | grep -oE '\([0-9]+ smoothing sweeps' | grep -oE '[0-9]+')
        else
            it=$(echo "$line" | grep -oE 'finished in [0-9]+ iterations' | grep -oE '[0-9]+')
            sweeps=$it
        fi
        if [ -z "$it" ] || [ -z "$sec" ]; then echo "ERROR: parse: $line" >&2; exit 1; fi
        secs="$secs $sec"
    done
    local med
    med=$(med3 $secs)
    echo "$np $solver $it $sweeps $med"
}

echo "== AR007 mg2 benchmark (${NX}^2 manufactured solution, tol 1e-2, OMP=1, median of $REPS) =="
echo "np solver iterations smoothing_sweeps wall_sec"
for np in 1 4; do
    for solver in red_black_gs jacobi mg2; do
        run_case "$np" "$solver"
    done
done
echo "== done =="
