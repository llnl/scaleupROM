#!/usr/bin/env python3
"""Plot the most recent complete scaling campaigns from results.csv."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("--input", default="benchmarks/results.csv")
parser.add_argument("--output-dir", default="benchmarks/plots")
parser.add_argument("--nx", type=int, default=128)
parser.add_argument("--ny", type=int, default=128)
parser.add_argument("--repeats", type=int, default=3,
                    help="Use the last N records for each rank count")
parser.add_argument("--weak-blocks-per-rank", type=int, default=80,
                    help="Weak-scaling campaign to plot (default: 80)")
args = parser.parse_args()

groups = defaultdict(list)
with open(args.input, newline="") as stream:
    for row in csv.DictReader(stream):
        if row["mode"] == "strong" and int(row["nx"]) == args.nx and int(row["ny"]) == args.ny:
            groups[int(row["ranks"])].append(row)

if not groups:
    raise SystemExit(f"No strong-scaling records for {args.nx}x{args.ny}")

ranks = np.array(sorted(groups))
records = {rank: groups[rank][-args.repeats:] for rank in ranks}


def values(field):
    return [np.array([float(row[field]) for row in records[rank]]) for rank in ranks]


def medians(field):
    return np.array([np.median(item) for item in values(field)])


wall = medians("wall_s")
kernel = medians("kernel_s")
assembly = medians("assemble_s")
solve = medians("solve_s")
initialization = wall - kernel
wall_speedup = wall[0] / wall
kernel_speedup = kernel[0] / kernel
wall_efficiency = 100.0 * wall_speedup / ranks
kernel_efficiency = 100.0 * kernel_speedup / ranks
wall_rate = medians("wall_blocks_per_s")
kernel_rate = medians("kernel_blocks_per_s")

output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})

fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), constrained_layout=True)
ax = axes[0, 0]
ax.plot(ranks, wall, "o-", label="End-to-end")
ax.plot(ranks, kernel, "s-", label="ROM kernel")
ax.set(xlabel="MPI ranks", ylabel="Median time (s)", title="Runtime")
ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xticks(ranks, labels=ranks)
ax.legend()

ax = axes[0, 1]
ax.plot(ranks, ranks, "k--", label="Ideal")
ax.plot(ranks, wall_speedup, "o-", label="End-to-end")
ax.plot(ranks, kernel_speedup, "s-", label="ROM kernel")
ax.set(xlabel="MPI ranks", ylabel="Speedup", title="Strong-scaling speedup")
ax.set_xscale("log", base=2)
ax.set_xticks(ranks, labels=ranks)
ax.legend()

ax = axes[1, 0]
ax.plot(ranks, wall_efficiency, "o-", label="End-to-end")
ax.plot(ranks, kernel_efficiency, "s-", label="ROM kernel")
ax.axhline(100, color="k", linestyle="--", linewidth=1)
ax.set(xlabel="MPI ranks", ylabel="Parallel efficiency (%)", title="Parallel efficiency", ylim=(0, 110))
ax.set_xscale("log", base=2)
ax.set_xticks(ranks, labels=ranks)
ax.legend()

ax = axes[1, 1]
ax.plot(ranks, wall_rate, "o-", label="End-to-end")
ax.plot(ranks, kernel_rate, "s-", label="ROM kernel")
ax.set(xlabel="MPI ranks", ylabel="Blocks/s", title="Throughput")
ax.set_xscale("log", base=2)
ax.set_xticks(ranks, labels=ranks)
ax.legend()

fig.suptitle(f"Poisson ROM strong scaling: {args.nx}x{args.ny} blocks")
fig.savefig(output_dir / f"strong_scaling_{args.nx}x{args.ny}.png", dpi=180)
plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
ax.bar(ranks, initialization, label="Initialization and other untimed work")
ax.bar(ranks, assembly, bottom=initialization, label="ROM assembly")
ax.bar(ranks, solve, bottom=initialization + assembly, label="ROM solve")
ax.set(xlabel="MPI ranks", ylabel="Median wall time (s)",
       title=f"Poisson ROM runtime breakdown: {args.nx}x{args.ny} blocks")
ax.set_xticks(ranks)
ax.legend()
fig.savefig(output_dir / f"strong_scaling_breakdown_{args.nx}x{args.ny}.png", dpi=180)
plt.close(fig)

print(f"Wrote plots to {output_dir}")

# Single-rank throughput sweep. Use the latest samples for every square size.
throughput_groups = defaultdict(list)
with open(args.input, newline="") as stream:
    for row in csv.DictReader(stream):
        if row["mode"] == "throughput" and int(row["ranks"]) == 1:
            throughput_groups[(int(row["nx"]), int(row["ny"]))].append(row)

if throughput_groups:
    shapes = sorted(throughput_groups, key=lambda shape: shape[0] * shape[1])
    trecords = {shape: throughput_groups[shape][-args.repeats:] for shape in shapes}

    def tmedians(field):
        return np.array([
            np.median([float(row[field]) for row in trecords[shape]])
            for shape in shapes
        ])

    blocks = np.array([nx * ny for nx, ny in shapes])
    dofs = tmedians("rom_dofs")
    twall = tmedians("wall_s")
    tkernel = tmedians("kernel_s")
    tassembly = tmedians("assemble_s")
    tsolve = tmedians("solve_s")
    tinit = twall - tkernel
    twall_rate = tmedians("wall_blocks_per_s")
    tkernel_rate = tmedians("kernel_blocks_per_s")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(dofs, twall, "o-", label="End-to-end")
    ax.plot(dofs, tkernel, "s-", label="ROM kernel")
    ax.set(xlabel="ROM unknowns", ylabel="Median time (s)", title="Runtime")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()

    ax = axes[0, 1]
    ax.plot(dofs, twall_rate, "o-", label="End-to-end")
    ax.plot(dofs, tkernel_rate, "s-", label="ROM kernel")
    ax.set(xlabel="ROM unknowns", ylabel="Blocks/s", title="Single-rank throughput")
    ax.set_xscale("log")
    ax.legend()

    ax = axes[1, 0]
    ax.plot(dofs, 1000 * twall / blocks, "o-", label="End-to-end")
    ax.plot(dofs, 1000 * tkernel / blocks, "s-", label="ROM kernel")
    ax.set(xlabel="ROM unknowns", ylabel="Time per block (ms)", title="Cost per block")
    ax.set_xscale("log")
    ax.legend()

    ax = axes[1, 1]
    ax.plot(dofs, 100 * tinit / twall, "o-", label="Initialization/other")
    ax.plot(dofs, 100 * tassembly / twall, "s-", label="ROM assembly")
    ax.plot(dofs, 100 * tsolve / twall, "^-", label="ROM solve")
    ax.set(xlabel="ROM unknowns", ylabel="Fraction of wall time (%)",
           title="Runtime composition", ylim=(0, 105))
    ax.set_xscale("log")
    ax.legend()

    fig.suptitle("Poisson ROM single-rank throughput scaling")
    fig.savefig(output_dir / "throughput_scaling.png", dpi=180)
    plt.close(fig)

    labels = [f"{nx}x{ny}" for nx, ny in shapes]
    x = np.arange(len(shapes))
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.bar(x, tinit, label="Initialization and other untimed work")
    ax.bar(x, tassembly, bottom=tinit, label="ROM assembly")
    ax.bar(x, tsolve, bottom=tinit + tassembly, label="ROM solve")
    ax.set(xlabel="Block grid", ylabel="Median wall time (s)",
           title="Poisson ROM single-rank runtime breakdown")
    ax.set_xticks(x, labels)
    ax.legend()
    fig.savefig(output_dir / "throughput_breakdown.png", dpi=180)
    plt.close(fig)

    # A focused view of the practical, complete-run throughput metric.
    wall_rate_samples = [
        np.array([float(row["wall_blocks_per_s"]) for row in trecords[shape]])
        for shape in shapes
    ]
    wall_rate_low = np.array([np.median(sample) - np.min(sample) for sample in wall_rate_samples])
    wall_rate_high = np.array([np.max(sample) - np.median(sample) for sample in wall_rate_samples])
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.errorbar(dofs, twall_rate, yerr=[wall_rate_low, wall_rate_high],
                marker="o", linewidth=2, capsize=4, color="#1f77b4")
    for xvalue, yvalue, label in zip(dofs, twall_rate, labels):
        ax.annotate(label, (xvalue, yvalue), xytext=(0, 9),
                    textcoords="offset points", ha="center", fontsize=9)
    ax.axvspan(50000, 100000, color="#ff7f0e", alpha=0.12,
               label="Expected saturation range")
    ax.set(xlabel="ROM unknowns (10 per component block)",
           ylabel="End-to-end throughput (blocks/s)",
           title="Poisson ROM single-rank end-to-end throughput")
    ax.set_xscale("log")
    ax.set_ylim(bottom=0)
    ax.legend()
    fig.savefig(output_dir / "end_to_end_throughput.png", dpi=180)
    plt.close(fig)

    print(f"Wrote throughput plots for {len(shapes)} problem sizes")

# Weak-scaling sweep. Each rank receives the same number of component blocks.
weak_groups = defaultdict(list)
with open(args.input, newline="") as stream:
    for row in csv.DictReader(stream):
        blocks_per_rank = int(row["blocks"]) // int(row["ranks"])
        if row["mode"] == "weak" and blocks_per_rank == args.weak_blocks_per_rank:
            weak_groups[int(row["ranks"])].append(row)

if weak_groups:
    wranks = np.array(sorted(weak_groups))
    wrecords = {rank: weak_groups[rank][-args.repeats:] for rank in wranks}

    def wsamples(field):
        return [
            np.array([float(row[field]) for row in wrecords[rank]])
            for rank in wranks
        ]

    def wmedians(field):
        return np.array([np.median(sample) for sample in wsamples(field)])

    def wrange(field):
        samples = wsamples(field)
        centers = np.array([np.median(sample) for sample in samples])
        return np.array([
            centers - np.array([np.min(sample) for sample in samples]),
            np.array([np.max(sample) for sample in samples]) - centers,
        ])

    wblocks = wmedians("blocks")
    wwall = wmedians("wall_s")
    wkernel = wmedians("kernel_s")
    wassembly = wmedians("assemble_s")
    wsolve = wmedians("solve_s")
    wother = wwall - wkernel
    wwall_efficiency = 100.0 * wwall[0] / wwall
    wkernel_efficiency = 100.0 * wkernel[0] / wkernel
    wwall_rate_per_rank = wmedians("wall_blocks_per_s") / wranks
    wkernel_rate_per_rank = wmedians("kernel_blocks_per_s") / wranks
    wlabels = [
        f'{int(float(wrecords[rank][-1]["nx"]))}x'
        f'{int(float(wrecords[rank][-1]["ny"]))}'
        for rank in wranks
    ]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), constrained_layout=True)
    ax = axes[0, 0]
    ax.errorbar(wranks, wwall, yerr=wrange("wall_s"), marker="o",
                linewidth=2, capsize=4, label="End-to-end")
    ax.errorbar(wranks, wkernel, yerr=wrange("kernel_s"), marker="s",
                linewidth=2, capsize=4, label="ROM kernel")
    ax.axhline(wwall[0], color="#1f77b4", linestyle="--", alpha=0.5)
    ax.axhline(wkernel[0], color="#ff7f0e", linestyle="--", alpha=0.5)
    ax.set(xlabel="MPI ranks", ylabel="Median time (s)", title="Runtime (ideal is flat)")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(wranks, labels=wranks)
    ax.legend()

    ax = axes[0, 1]
    ax.plot(wranks, wwall_efficiency, "o-", label="End-to-end")
    ax.plot(wranks, wkernel_efficiency, "s-", label="ROM kernel")
    ax.axhline(100, color="k", linestyle="--", linewidth=1, label="Ideal")
    ax.set(xlabel="MPI ranks", ylabel="Weak-scaling efficiency (%)",
           title="Efficiency ($T_1/T_p$)", ylim=(0, 135))
    ax.set_xscale("log", base=2)
    ax.set_xticks(wranks, labels=wranks)
    ax.legend()

    ax = axes[1, 0]
    ax.plot(wranks, wwall_rate_per_rank, "o-", label="End-to-end")
    ax.plot(wranks, wkernel_rate_per_rank, "s-", label="ROM kernel")
    ax.axhline(wwall_rate_per_rank[0], color="#1f77b4", linestyle="--", alpha=0.5)
    ax.axhline(wkernel_rate_per_rank[0], color="#ff7f0e", linestyle="--", alpha=0.5)
    ax.set(xlabel="MPI ranks", ylabel="Blocks/s/rank",
           title="Per-rank throughput (ideal is flat)")
    ax.set_xscale("log", base=2)
    ax.set_xticks(wranks, labels=wranks)
    ax.legend()

    ax = axes[1, 1]
    x = np.arange(len(wranks))
    ax.bar(x, wother, label="Initialization and other untimed work")
    ax.bar(x, wassembly, bottom=wother, label="ROM assembly")
    ax.bar(x, wsolve, bottom=wother + wassembly, label="ROM solve")
    ax.set(xlabel="MPI ranks (block grid)", ylabel="Median wall time (s)",
           title="Runtime breakdown")
    ax.set_xticks(x, [f"{rank}\n{label}" for rank, label in zip(wranks, wlabels)])
    ax.legend(fontsize=8)

    fig.suptitle(
        "Poisson ROM weak scaling: "
        f"{args.weak_blocks_per_rank} blocks and "
        f"{args.weak_blocks_per_rank * 10} ROM unknowns per rank"
    )
    weak_filename = ("weak_scaling.png" if args.weak_blocks_per_rank == 80
                     else f"weak_scaling_{args.weak_blocks_per_rank}.png")
    fig.savefig(output_dir / weak_filename, dpi=180)
    plt.close(fig)

    print(f"Wrote weak-scaling plot for {len(wranks)} rank counts")
