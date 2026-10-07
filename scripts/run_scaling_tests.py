#!/usr/bin/env python3
"""
Automated scaling test runner for HyPoS.
Runs strong and weak scaling tests with varying MPI process counts.
"""
import subprocess
import json
import os
import sys
from datetime import datetime
import argparse


def run_test(nx, ny, nz, mpi_procs, omp_threads, output_dir, solver="jacobi", max_iter=10000, tol="1e-6"):
    """Run a single HyPoS instance and collect performance metrics."""
    cmd = [
        "mpirun", "-np", str(mpi_procs),
        "./hypos",
        "--nx", str(nx), "--ny", str(ny), "--nz", str(nz),
        "--omp-threads", str(omp_threads),
        "--solver", solver,
        "--max-iter", str(max_iter),
        "--tol", str(tol),
        "--output-format", "json",
        "--output-dir", output_dir
    ]
    
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(omp_threads)
    
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        return None
    
    # Read performance report
    report_path = os.path.join(output_dir, "performance_report.json")
    if os.path.exists(report_path):
        with open(report_path, "r") as f:
            return json.load(f)
    return None


def strong_scaling_test(output_dir, nx=1024, ny=1024, nz=1, max_procs=16, omp_threads_override=0,
                        max_iter=10000, tol="1e-6"):
    """Strong scaling: fixed problem size, increasing processes."""
    print(f"\n{'='*60}")
    print("STRONG SCALING TEST")
    print(f"Fixed grid: {nx}x{ny}x{nz}, varying MPI processes")
    print(f"{'='*60}")
    
    procs_list = [1, 2, 4, 8, 16]
    procs_list = [p for p in procs_list if p <= max_procs]
    
    results = []
    base_time = None
    
    for procs in procs_list:
        if omp_threads_override > 0:
            omp_threads = omp_threads_override
        else:
            omp_threads = max(1, 4 // procs)  # Adjust threads per process
        
        test_dir = os.path.join(output_dir, f"strong_p{procs}")
        os.makedirs(test_dir, exist_ok=True)
        
        print(f"\nRunning with {procs} MPI processes, {omp_threads} OMP threads each...")
        report = run_test(nx, ny, nz, procs, omp_threads, test_dir, max_iter=max_iter, tol=tol)
        
        if report:
            total_time = report["performance"]["total_time_sec"]
            iter_time = report["performance"]["iter_time_ms"]
            
            if base_time is None:
                base_time = total_time
            
            speedup = base_time / total_time if total_time > 0 else 0
            efficiency = speedup / procs
            
            result = {
                "procs": procs,
                "omp_threads": omp_threads,
                "total_time_sec": total_time,
                "iter_time_ms": iter_time,
                "speedup": speedup,
                "efficiency": efficiency,
                "iterations": report["performance"]["iterations"]
            }
            results.append(result)
            
            print(f"  Total time: {total_time:.4f}s, Speedup: {speedup:.2f}, Efficiency: {efficiency:.2%}")
    
    # Save results
    with open(os.path.join(output_dir, "strong_scaling.json"), "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\nStrong scaling results saved to {os.path.join(output_dir, 'strong_scaling.json')}")
    return results


def weak_scaling_test(output_dir, per_proc_n=512, max_procs=16, omp_threads_override=0,
                      max_iter=10000, tol="1e-6"):
    """Weak scaling: fixed per-process workload, increasing processes."""
    print(f"\n{'='*60}")
    print("WEAK SCALING TEST")
    print(f"Fixed per-process grid: {per_proc_n}x{per_proc_n}x1, varying MPI processes")
    print(f"{'='*60}")
    
    procs_list = [1, 2, 4, 8, 16]
    procs_list = [p for p in procs_list if p <= max_procs]
    
    results = []
    base_time = None
    
    for procs in procs_list:
        nx = per_proc_n * procs  # Increase total grid with procs
        ny = per_proc_n
        nz = 1
        omp_threads = omp_threads_override if omp_threads_override > 0 else 1
        
        test_dir = os.path.join(output_dir, f"weak_p{procs}")
        os.makedirs(test_dir, exist_ok=True)
        
        print(f"\nRunning with {procs} MPI processes, total grid {nx}x{ny}x{nz}...")
        report = run_test(nx, ny, nz, procs, omp_threads, test_dir, max_iter=max_iter, tol=tol)
        
        if report:
            total_time = report["performance"]["total_time_sec"]
            iter_time = report["performance"]["iter_time_ms"]
            
            if base_time is None:
                base_time = total_time
            
            efficiency = base_time / total_time if total_time > 0 else 0
            
            result = {
                "procs": procs,
                "total_nx": nx,
                "total_ny": ny,
                "total_time_sec": total_time,
                "iter_time_ms": iter_time,
                "efficiency": efficiency,
                "iterations": report["performance"]["iterations"]
            }
            results.append(result)
            
            print(f"  Total time: {total_time:.4f}s, Efficiency: {efficiency:.2%}")
    
    # Save results
    with open(os.path.join(output_dir, "weak_scaling.json"), "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\nWeak scaling results saved to {os.path.join(output_dir, 'weak_scaling.json')}")
    return results


def print_summary(strong_results, weak_results):
    """Print formatted summary table."""
    print(f"\n{'='*60}")
    print("SCALING TEST SUMMARY")
    print(f"{'='*60}")
    
    if strong_results:
        print("\n--- Strong Scaling ---")
        print(f"{'Procs':<8} {'Threads':<8} {'Time(s)':<12} {'Speedup':<10} {'Efficiency':<12}")
        print("-" * 50)
        for r in strong_results:
            print(f"{r['procs']:<8} {r['omp_threads']:<8} {r['total_time_sec']:<12.4f} "
                  f"{r['speedup']:<10.2f} {r['efficiency']:<12.2%}")
    
    if weak_results:
        print("\n--- Weak Scaling ---")
        print(f"{'Procs':<8} {'Grid':<16} {'Time(s)':<12} {'Efficiency':<12}")
        print("-" * 50)
        for r in weak_results:
            print(f"{r['procs']:<8} {r['total_nx']}x{r['total_ny']:<8} "
                  f"{r['total_time_sec']:<12.4f} {r['efficiency']:<12.2%}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="HyPoS Scaling Test Runner")
    parser.add_argument("--output-dir", default="./scaling_results", help="Output directory for results")
    parser.add_argument("--max-procs", type=int, default=16, help="Maximum MPI processes to test")
    parser.add_argument("--strong-nx", type=int, default=1024, help="Grid size for strong scaling")
    parser.add_argument("--weak-n", type=int, default=512, help="Per-process grid size for weak scaling")
    parser.add_argument("--omp-threads", type=int, default=0,
                        help="Fixed OpenMP threads per process (0=auto: strong=4//procs, weak=1)")
    parser.add_argument("--max-iter", type=int, default=10000, help="Maximum solver iterations per run")
    parser.add_argument("--tol", default="1e-6", help="Solver convergence tolerance")
    parser.add_argument("--strong-only", action="store_true", help="Run only strong scaling")
    parser.add_argument("--weak-only", action="store_true", help="Run only weak scaling")
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    strong_results = None
    weak_results = None
    
    if not args.weak_only:
        strong_results = strong_scaling_test(args.output_dir, 
                                              nx=args.strong_nx, ny=args.strong_nx, 
                                              max_procs=args.max_procs,
                                              omp_threads_override=args.omp_threads,
                                              max_iter=args.max_iter, tol=args.tol)
    
    if not args.strong_only:
        weak_results = weak_scaling_test(args.output_dir, 
                                          per_proc_n=args.weak_n, 
                                          max_procs=args.max_procs,
                                          omp_threads_override=args.omp_threads,
                                          max_iter=args.max_iter, tol=args.tol)
    
    print_summary(strong_results, weak_results)
