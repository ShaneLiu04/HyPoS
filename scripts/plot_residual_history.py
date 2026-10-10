#!/usr/bin/env python3
"""
Plot residual histories from HyPoS `--residual-history` CSV files.
Overlays one log-scale convergence curve per CSV on a single figure
(AR009, design FP6/D4 — matplotlib convention follows plot_scaling.py).

CSV contract (design §4.3): first line `iteration,residual`, then one
`int,double` row per recorded iteration (row condition iteration % k == 0
for the --residual-check-interval k; the CG family records every
iteration).

Usage:
    python3 plot_residual_history.py --input a.csv b.csv ... \
        [--labels cg pcg ...] [--output-dir ./residual_plots]
"""
import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_residual_csv(path):
    """Read one residual-history CSV into (iterations, residuals) lists."""
    iterations = []
    residuals = []
    with open(path, "r") as f:
        header = f.readline().strip()
        if header != "iteration,residual":
            raise ValueError(
                "%s: bad header '%s' (expected iteration,residual)"
                % (path, header))
        for line in f:
            line = line.strip()
            if not line:
                continue
            it_s, res_s = line.split(",")
            iterations.append(int(it_s))
            residuals.append(float(res_s))
    return iterations, residuals


def plot_histories(files, labels, output_path):
    """Overlay the convergence curves of all CSVs, log-scale y."""
    fig, ax = plt.subplots(figsize=(8, 6))
    for path, label in zip(files, labels):
        iterations, residuals = read_residual_csv(path)
        ax.plot(iterations, residuals, '-', linewidth=2, label=label)

    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('L2 residual', fontsize=12)
    ax.set_title('Residual history', fontsize=14)
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print("Residual history plot saved to %s" % output_path)
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot HyPoS residual histories (--residual-history CSVs)")
    parser.add_argument("--input", nargs='+', required=True,
                        help="One or more residual CSV files")
    parser.add_argument("--labels", nargs='*', default=None,
                        help="Curve labels (defaults to the file names)")
    parser.add_argument("--output-dir", default="./residual_plots",
                        help="Directory for the plot")
    args = parser.parse_args()

    labels = args.labels if args.labels else \
        [os.path.splitext(os.path.basename(p))[0] for p in args.input]
    if len(labels) != len(args.input):
        raise SystemExit("error: --labels count must match --input count")

    os.makedirs(args.output_dir, exist_ok=True)
    plot_histories(args.input, labels,
                   os.path.join(args.output_dir, "residual_history.png"))
