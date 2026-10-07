#!/usr/bin/env python3
"""
Plot scaling results from HyPoS scaling tests.
Generates publication-quality scaling curves.
"""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import os


def plot_strong_scaling(results, output_path):
    """Plot strong scaling curve."""
    procs = [r["procs"] for r in results]
    times = [r["total_time_sec"] for r in results]
    speedups = [r["speedup"] for r in results]
    efficiencies = [r["efficiency"] for r in results]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Speedup plot
    ax1.plot(procs, speedups, 'o-', linewidth=2, markersize=8, label='Actual Speedup')
    ax1.plot(procs, procs, '--', linewidth=1, color='gray', label='Ideal Speedup')
    ax1.set_xlabel('MPI Processes', fontsize=12)
    ax1.set_ylabel('Speedup', fontsize=12)
    ax1.set_title('Strong Scaling — Speedup', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xscale('log', base=2)
    ax1.set_yscale('log', base=2)
    
    # Efficiency plot
    ax2.plot(procs, efficiencies, 's-', linewidth=2, markersize=8, color='green')
    ax2.axhline(y=1.0, linestyle='--', color='gray', linewidth=1)
    ax2.set_xlabel('MPI Processes', fontsize=12)
    ax2.set_ylabel('Parallel Efficiency', fontsize=12)
    ax2.set_title('Strong Scaling — Efficiency', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, 1.1])
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Strong scaling plot saved to {output_path}")
    plt.close()


def plot_weak_scaling(results, output_path):
    """Plot weak scaling curve."""
    procs = [r["procs"] for r in results]
    efficiencies = [r["efficiency"] for r in results]
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.plot(procs, efficiencies, 'o-', linewidth=2, markersize=10, color='blue')
    ax.axhline(y=1.0, linestyle='--', color='gray', linewidth=1, label='Ideal Efficiency')
    ax.set_xlabel('MPI Processes', fontsize=12)
    ax.set_ylabel('Parallel Efficiency', fontsize=12)
    ax.set_title('Weak Scaling — Efficiency', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.1])
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Weak scaling plot saved to {output_path}")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot HyPoS scaling results")
    parser.add_argument("--input-dir", default="./scaling_results", help="Directory with JSON results")
    parser.add_argument("--output-dir", default="./scaling_plots", help="Directory for plots")
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Strong scaling
    strong_file = os.path.join(args.input_dir, "strong_scaling.json")
    if os.path.exists(strong_file):
        with open(strong_file, "r") as f:
            strong_results = json.load(f)
        plot_strong_scaling(strong_results, 
                            os.path.join(args.output_dir, "strong_scaling.png"))
    
    # Weak scaling
    weak_file = os.path.join(args.input_dir, "weak_scaling.json")
    if os.path.exists(weak_file):
        with open(weak_file, "r") as f:
            weak_results = json.load(f)
        plot_weak_scaling(weak_results,
                          os.path.join(args.output_dir, "weak_scaling.png"))
