#!/bin/bash
set -e

# HyPoS Benchmark Script
# Runs a quick benchmark with different configurations

echo "========================================"
echo "HyPoS Benchmark Suite"
echo "========================================"

BUILD_DIR="build"
if [ ! -d "$BUILD_DIR" ]; then
    echo "Error: Build directory not found. Run: cmake -B build && cmake --build build"
    exit 1
fi

HYPO="$BUILD_DIR/hypos"
if [ ! -f "$HYPO" ]; then
    echo "Error: hypos executable not found. Build first."
    exit 1
fi

mkdir -p benchmark_results

# Test 1: Small serial test
echo ""
echo "Test 1: Small serial (128x128)"
mpirun -np 1 "$HYPO" --nx 128 --ny 128 --max-iter 1000 --output-dir benchmark_results/serial

# Test 2: Small parallel test (4 procs)
echo ""
echo "Test 2: Small parallel 4 procs (256x256)"
mpirun -np 4 "$HYPO" --nx 256 --ny 256 --max-iter 1000 --output-dir benchmark_results/parallel4

# Test 3: Medium scale
echo ""
echo "Test 3: Medium scale 8 procs (512x512)"
mpirun -np 8 "$HYPO" --nx 512 --ny 512 --max-iter 5000 --enable-profiling --output-dir benchmark_results/medium8

# Test 4: With OpenMP
echo ""
echo "Test 4: Medium with OpenMP (4 procs x 2 threads)"
OMP_NUM_THREADS=2 mpirun -np 4 "$HYPO" --nx 512 --ny 512 --omp-threads 2 --max-iter 5000 --output-dir benchmark_results/medium_omp

echo ""
echo "========================================"
echo "Benchmark complete. Results in benchmark_results/"
echo "========================================"
