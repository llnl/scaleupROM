#!/usr/bin/env python3
"""Aggregate PARALLEL_PROFILE records from benchmark logs and plot hot phases."""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("--log-dir", default="benchmarks/logs")
parser.add_argument("--results", default="benchmarks/results.csv",
                    help="Benchmark CSV providing end-to-end wall times")
parser.add_argument("--csv", default="benchmarks/profile_summary.csv")
parser.add_argument("--plot-dir", default="benchmarks/plots")
parser.add_argument("--top", type=int, default=8,
                    help="Number of hottest phases to show per campaign")
args = parser.parse_args()

log_name = re.compile(
    r"^(?P<mode>strong|weak|throughput)\.(?P<nx>\d+)x(?P<ny>\d+)"
    r"\.n(?P<ranks>\d+)\.r(?P<repeat>\d+)\.log$"
)
record = re.compile(
    r"PARALLEL_PROFILE,(?P<rank>\d+),(?P<nranks>\d+),"
    r"(?P<phase>[^,]+),(?P<seconds>[^,]+),(?P<work>-?\d+)$"
)

groups = defaultdict(list)
for path in sorted(Path(args.log_dir).glob("*.log")):
    case = log_name.match(path.name)
    if not case:
        continue
    with path.open(errors="replace") as stream:
        for line in stream:
            # Some application progress messages intentionally omit a newline,
            # so a profiling record may not begin at column zero.
            item = record.search(line)
            if item:
                key = (
                    case["mode"], int(case["nx"]), int(case["ny"]),
                    int(case["ranks"]), int(case["repeat"]), item["phase"],
                )
                groups[key].append((
                    int(item["rank"]), float(item["seconds"]), int(item["work"])
                ))

if not groups:
    raise SystemExit(
        "No profiling records found. Run with "
        "'-f time_profile/parallel_scaling=true'."
    )

# Logs are overwritten when a case is rerun, while results are appended. Keep
# the last result for each case/repetition so it corresponds to the current log.
wall_times = {}
with Path(args.results).open(newline="") as stream:
    for item in csv.DictReader(stream):
        key = (
            item["mode"], int(item["nx"]), int(item["ny"]),
            int(item["ranks"]), int(item["repeat"]),
        )
        wall_times[key] = float(item["wall_s"])

rows = []
for key, samples in sorted(groups.items()):
    mode, nx, ny, ranks, repeat, phase = key
    seconds = np.array([sample[1] for sample in samples])
    wall = wall_times.get((mode, nx, ny, ranks, repeat), np.nan)
    max_seconds = seconds.max()
    rows.append({
        "mode": mode, "nx": nx, "ny": ny, "ranks": ranks,
        "repeat": repeat, "phase": phase,
        "rank_samples": len({sample[0] for sample in samples}),
        "min_s": seconds.min(), "mean_s": seconds.mean(),
        "max_s": max_seconds,
        "wall_s": wall,
        "percent_wall": 100.0 * max_seconds / wall,
        "min_work": min(sample[2] for sample in samples),
        "max_work": max(sample[2] for sample in samples),
    })

csv_path = Path(args.csv)
csv_path.parent.mkdir(parents=True, exist_ok=True)
with csv_path.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
print(f"Wrote {csv_path} ({len(rows)} phase/run records)")

plot_dir = Path(args.plot_dir)
plot_dir.mkdir(parents=True, exist_ok=True)
for mode in sorted({row["mode"] for row in rows}):
    mode_rows = [row for row in rows if row["mode"] == mode]
    phase_cost = defaultdict(list)
    for row in mode_rows:
        phase_cost[row["phase"]].append(row["max_s"])
    phases = sorted(phase_cost, key=lambda p: max(phase_cost[p]), reverse=True)[:args.top]
    ranks = sorted({row["ranks"] for row in mode_rows})

    for field, ylabel, suffix in (
        ("max_s", "Median of maximum-rank time (s)", "time"),
        ("percent_wall", "Median percentage of end-to-end wall time", "percent_wall"),
    ):
        fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
        for phase in phases:
            medians = []
            for nranks in ranks:
                values = [row[field] for row in mode_rows
                          if row["phase"] == phase and row["ranks"] == nranks]
                medians.append(np.median(values) if values else np.nan)
            ax.plot(ranks, medians, "o-", label=phase)
        ax.set(xlabel="MPI ranks", ylabel=ylabel,
               title=f"Poisson ROM {mode} profiling hot spots: {suffix.replace('_', ' ')}")
        ax.set_xscale("log", base=2)
        if field == "max_s":
            ax.set_yscale("log")
        else:
            ax.set_ylim(bottom=0)
        ax.set_xticks(ranks, labels=ranks)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        output = plot_dir / f"profile_{mode}_{suffix}.png"
        fig.savefig(output, dpi=180)
        plt.close(fig)
        print(f"Wrote {output}")

# Build an exclusive decomposition for strong and weak scaling. For each rank
# count, use the repetition with the median external wall time and the rank with
# the largest measured process time. Child timers are collapsed into disjoint
# categories; residuals expose work not covered by the detailed timers.
run_phases = defaultdict(lambda: defaultdict(dict))
for key, samples in groups.items():
    mode, nx, ny, ranks, repeat, phase = key
    run = (mode, nx, ny, ranks, repeat)
    for rank, seconds, _ in samples:
        run_phases[run][rank][phase] = seconds


def nonnegative(value):
    return max(0.0, value)


def exclusive_regions(run):
    ranks_data = run_phases[run]
    critical_rank = max(
        ranks_data,
        key=lambda rank: ranks_data[rank].get("app.process_total", -1.0),
    )
    phase = ranks_data[critical_rank]
    get = lambda name: phase.get(name, 0.0)
    wall = wall_times[run]
    process = get("app.process_total")
    workflow = get("app.workflow")

    topology = get("init.solver_and_topology")
    variable_init = get("init.variables_rom_visualization")
    rhs_basis_io = (
        get("init.parameterized_problem")
        + get("setup.rhs_operators_bc_assembly")
        + get("rom.load_reduced_basis")
        + get("rom.load_component_operators")
    )
    local_port = (
        get("rom.assemble_local_domains")
        + get("rom.assemble_global_port_loop")
    )
    block_bookkeeping = nonnegative(
        get("rom.block_assembly") - local_port
    ) + get("rom.block_finalize")
    monolithic = get("rom.set_matrix_monolithic")
    hypre = nonnegative(
        get("rom.create_hypre_total") - get("rom.mumps_setup")
    )
    mumps_setup = get("rom.mumps_setup")
    rhs_projection = get("rom.project_rhs")
    solve = sum(get(name) for name in (
        "rom.solve_extract_local_rhs", "rom.solve_mumps",
        "rom.solve_gather_global_solution", "rom.solve_lift_solution",
    ))
    postprocess = (
        get("postprocess.output_and_comparison")
        + get("cleanup.solver_and_problem")
    )

    detailed_workflow = sum((
        topology, variable_init, rhs_basis_io, local_port, block_bookkeeping,
        monolithic, hypre, mumps_setup, rhs_projection, solve, postprocess,
    ))
    return {
        "Launcher outside process": nonnegative(wall - process),
        "MPI/config/process overhead": nonnegative(process - workflow),
        "Topology and solver init": topology,
        "Variable and ROM init": variable_init,
        "RHS, basis, and operator I/O": rhs_basis_io,
        "Local and port assembly": local_port,
        "Block bookkeeping/finalize": block_bookkeeping,
        "Global monolithic conversion": monolithic,
        "Hypre construction": hypre,
        "MUMPS setup": mumps_setup,
        "RHS projection": rhs_projection,
        "ROM solve and reconstruction": solve,
        "Postprocess and cleanup": postprocess,
        "Other workflow work": nonnegative(workflow - detailed_workflow),
    }


exclusive_rows = []
for mode in ("strong", "weak"):
    mode_runs = [run for run in run_phases
                 if run[0] == mode
                 and any("app.process_total" in phases
                         for phases in run_phases[run].values())]
    if not mode_runs:
        continue
    selected = []
    for nranks in sorted({run[3] for run in mode_runs}):
        candidates = [run for run in mode_runs if run[3] == nranks]
        candidates.sort(key=lambda run: wall_times[run])
        selected.append(candidates[len(candidates) // 2])

    decompositions = [exclusive_regions(run) for run in selected]
    labels = [str(run[3]) for run in selected]
    categories = list(decompositions[0])
    for run, decomposition in zip(selected, decompositions):
        wall = wall_times[run]
        for category, seconds in decomposition.items():
            exclusive_rows.append({
                "mode": mode, "nx": run[1], "ny": run[2],
                "ranks": run[3], "repeat": run[4], "wall_s": wall,
                "region": category, "seconds": seconds,
                "percent_wall": 100.0 * seconds / wall,
            })
    x = np.arange(len(selected))
    bottoms = np.zeros(len(selected))
    fig, ax = plt.subplots(figsize=(12, 7), constrained_layout=True)
    for category in categories:
        values = np.array([item[category] for item in decompositions])
        if np.all(values == 0):
            continue
        ax.bar(x, values, bottom=bottoms, label=category)
        bottoms += values
    for xpos, total in zip(x, bottoms):
        ax.text(xpos, total, f"{total:.2f} s", ha="center", va="bottom")
    ax.set(xlabel="MPI ranks", ylabel="End-to-end wall time (s)",
           title=f"Poisson ROM {mode}: exclusive wall-time decomposition")
    ax.set_xticks(x, labels)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8, ncol=2)
    output = plot_dir / f"profile_{mode}_exclusive_wall.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(f"Wrote {output}")

    bottoms = np.zeros(len(selected))
    fig, ax = plt.subplots(figsize=(12, 7), constrained_layout=True)
    for category in categories:
        values = np.array([
            100.0 * item[category] / wall_times[run]
            for item, run in zip(decompositions, selected)
        ])
        if np.all(values == 0):
            continue
        ax.bar(x, values, bottom=bottoms, label=category)
        bottoms += values
    ax.set(xlabel="MPI ranks", ylabel="Percentage of end-to-end wall time",
           title=f"Poisson ROM {mode}: exclusive wall-time percentages",
           ylim=(0, 100))
    ax.set_xticks(x, labels)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8, ncol=2)
    output = plot_dir / f"profile_{mode}_exclusive_percent_wall.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(f"Wrote {output}")

if exclusive_rows:
    exclusive_csv = csv_path.with_name("profile_exclusive_summary.csv")
    with exclusive_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=exclusive_rows[0].keys())
        writer.writeheader()
        writer.writerows(exclusive_rows)
    print(f"Wrote {exclusive_csv}")
