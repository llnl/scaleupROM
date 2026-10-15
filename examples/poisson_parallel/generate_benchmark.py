#!/usr/bin/env python3
"""Generate rectangular Poisson component configurations for benchmarking."""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "../../../utils/python")
from config import Empty, ObjectInSpace  # noqa: E402
from box_channel_config import BoxChannelConfig  # noqa: E402


def dimensions(spec):
    if "x" in spec.lower():
        nx, ny = map(int, spec.lower().split("x", 1))
    else:
        blocks = int(spec)
        nx = int(blocks**0.5)
        while blocks % nx:
            nx -= 1
        ny = blocks // nx
    if nx < 1 or ny < 1:
        raise argparse.ArgumentTypeError("dimensions must be positive")
    return nx, ny


parser = argparse.ArgumentParser(
    description="Generate random NxM block configurations (for example 64x64 or 4096)."
)
parser.add_argument("sizes", nargs="+", type=dimensions)
parser.add_argument("--seed", type=int, default=1)
parser.add_argument("--output-dir", default="configs/benchmarks")
args = parser.parse_args()

components = [
    Empty(),
    ObjectInSpace("square-circle"),
    ObjectInSpace("square-square"),
    ObjectInSpace("square-triangle"),
    ObjectInSpace("square-star"),
]
output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

for nx, ny in args.sizes:
    # Make each requested shape reproducible, independent of argument order.
    np.random.seed(args.seed + 1000003 * nx + 9176 * ny)
    config = BoxChannelConfig(nx, ny)
    for component in components:
        config.addComponent(component)
    output = output_dir / f"box-channel.{nx}x{ny}.h5"
    config.CreateRandomConfig(str(output))
    print(f"{nx * ny} blocks: {output}")
