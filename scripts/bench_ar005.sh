#!/bin/bash
# AR005 I/O benchmark: ASCII vs appended-binary .vti + elementwise vs
# block-write .bin (design §4.1 FP6, §6.5-3; decision D6 gate:
# median(after) <= median(baseline)).
#
# Usage: bench_ar005.sh baseline|after
#   baseline - run on the PRE-change tree (ASCII vtk, elementwise bin)
#   after    - run on the POST-change tree (appended vtk, block bin)
#
# Environment: WSL, Release build at /root/build-release, OMP_NUM_THREADS=1.
# io timing comes from the "io_write" profiler region (loop-path write via
# --save-interval 10; the final write happens after the profile report is
# printed and therefore does not enter it -- see main.cpp ordering).
set -u
MODE=${1:?usage: bench_ar005.sh baseline|after}
EXE=/root/build-release/hypos
OUT=/root/bench-ar005-$MODE
REPS=3

mkdir -p "$OUT"

# median of up to 3 values (mawk-safe: no asort)
med3() {
    echo "$1 $2 $3" | awk '{
        v1=$1; v2=$2; v3=$3;
        if (v1<=v2) { if (v2<=v3) print v2; else if (v1<=v3) print v3; else print v1 }
        else        { if (v1<=v3) print v1; else if (v2<=v3) print v3; else print v2 }
    }'
}

# run_case <nx> <fmt> -> echoes "nx fmt median_sec bytes"
run_case() {
    local nx=$1 fmt=$2 i d f size secs sec
    size=0
    secs=""
    for i in $(seq 1 $REPS); do
        d="$OUT/${nx}-${fmt}-$i"
        rm -rf "$d"; mkdir -p "$d"
        OMP_NUM_THREADS=1 "$EXE" --nx "$nx" --ny "$nx" --max-iter 10 \
            --save-interval 10 --output-format "$fmt" --enable-profiling \
            --output-dir "$d" > "$d/run.log" 2>&1
        if [ "$fmt" = "vtk" ]; then
            f=$(ls "$d"/solution_10_r0.vti 2>/dev/null)
        else
            f=$(ls "$d"/solution_10_r0.bin 2>/dev/null)
        fi
        if [ -z "$f" ]; then echo "ERROR: no output file in $d" >&2; exit 1; fi
        size=$(stat -c %s "$f")
        sec=$(grep -A2 '"io_write"' "$d/run.log" | grep -o '"total_sec": [0-9.eE+-]*' \
              | head -1 | sed 's/.*: //')
        if [ -z "$sec" ]; then echo "ERROR: no io_write region in $d/run.log" >&2; exit 1; fi
        secs="$secs $sec"
    done
    local med
    med=$(med3 $secs)
    echo "$nx $fmt $med $size"
}

echo "== AR005 I/O benchmark ($MODE, np=1, OMP=1, median of $REPS) =="
echo "nx format io_write_sec bytes"
for nx in 256 512; do
    run_case "$nx" vtk
    run_case "$nx" binary
done
echo "== done ($MODE) =="
