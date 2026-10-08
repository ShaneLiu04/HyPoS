#!/bin/bash
# AR003 ST manual executions: ST-005 (save-interval 50/120) and ST-003
# (mpibin vs binary data-area equivalence, np=1). Automated ST cases are
# covered by the 16 ctest entries; this script only covers the scenarios
# that need srs-literal parameters or cross-format comparison.
set -u
export OMP_NUM_THREADS=4
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
H=/root/build-release/hypos
R=/root/st_ar003
rm -rf "$R"
mkdir -p "$R"

echo "== ST-005: mpibin --save-interval 50 --max-iter 120, np=4 =="
mkdir -p "$R/s50"
mpirun.openmpi -n 4 --oversubscribe "$H" --output-format mpibin --nx 32 --ny 32 \
    --max-iter 120 --save-interval 50 --output-dir "$R/s50" > "$R/s50/run.log" 2>&1
echo "exit=$?"
ls "$R/s50" | sort

echo "== ST-003: mpibin vs binary data area, np=1 =="
mkdir -p "$R/mpibin" "$R/binary"
mpirun.openmpi -n 1 "$H" --output-format mpibin --nx 32 --ny 32 --max-iter 5 \
    --output-dir "$R/mpibin" > "$R/mpibin/run.log" 2>&1
mpirun.openmpi -n 1 "$H" --output-format binary --nx 32 --ny 32 --max-iter 5 \
    --output-dir "$R/binary" > "$R/binary/run.log" 2>&1
ls "$R/mpibin" "$R/binary"
# mpibin: 72-byte header; binary shard: 56-byte header (7 x Index)
dd if="$R/mpibin/solution_5.bin" of="$R/mpibin_data.bin" bs=8 skip=9 status=none
dd if="$R/binary/solution_5_r0.bin" of="$R/binary_data.bin" bs=8 skip=7 status=none
if cmp -s "$R/mpibin_data.bin" "$R/binary_data.bin"; then
    echo "ST-003: data areas byte-identical (32*32*8 = $((32*32*8)) bytes)"
else
    echo "ST-003: MISMATCH"
    cmp "$R/mpibin_data.bin" "$R/binary/solution_5_r0.bin" | head -2
fi
