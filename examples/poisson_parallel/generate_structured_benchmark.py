#!/usr/bin/env python3
"""Generate structured rectangular Poisson benchmark configurations in O(N)."""

import argparse
from pathlib import Path

import h5py
import numpy as np


COMPONENT_NAMES = [
    "empty",
    "square-circle",
    "square-square",
    "square-triangle",
    "square-star",
]


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


def generate(nx, ny, seed, output):
    np.random.seed(seed + 1000003 * nx + 9176 * ny)
    count = nx * ny

    # Match the random sequence used by BoxChannelConfig.CreateRandomConfig().
    original_types = np.array(
        [np.random.randint(len(COMPONENT_NAMES)) for _ in range(count)],
        dtype=np.int64,
    )
    used_types = sorted(set(original_types))
    type_map = {old: new for new, old in enumerate(used_types)}
    mesh_types = np.array([type_map[value] for value in original_types], dtype=np.int64)
    component_names = [COMPONENT_NAMES[value] for value in used_types]

    configuration = np.zeros((count, 6))
    for i in range(nx):
        start = i * ny
        configuration[start:start + ny, 0] = i
        configuration[start:start + ny, 1] = np.arange(ny)

    ports = {}
    interfaces = []

    def add_interface(mesh1, mesh2, attr1, attr2):
        name1 = COMPONENT_NAMES[original_types[mesh1]]
        name2 = COMPONENT_NAMES[original_types[mesh2]]
        key = (name1, name2, attr1, attr2)
        if key not in ports:
            ports[key] = len(ports)
        interfaces.append([mesh1, mesh2, attr1, attr2, ports[key]])

    # The original insertion order checks the left neighbor before the lower
    # neighbor. Preserve that ordering so seeded small cases compare exactly.
    for i in range(nx):
        for j in range(ny):
            mesh = i * ny + j
            if i > 0:
                add_interface((i - 1) * ny + j, mesh, 2, 4)
            if j > 0:
                add_interface(mesh, mesh - 1, 1, 3)

    boundaries = []
    for i in range(nx):
        for j in range(ny):
            mesh = i * ny + j
            if j == 0:
                boundaries.append([1, mesh, 1])
            if i == nx - 1:
                boundaries.append([2, mesh, 2])
            if j == ny - 1:
                boundaries.append([3, mesh, 3])
            if i == 0:
                boundaries.append([4, mesh, 4])
            if original_types[mesh] != 0:
                boundaries.append([5, mesh, 5])

    interface_data = np.asarray(interfaces, dtype=np.int64).reshape((-1, 5))
    boundary_data = np.asarray(boundaries, dtype=np.int64).reshape((-1, 3))

    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as stream:
        components = stream.create_group("components")
        components.attrs["number_of_components"] = len(component_names)
        for index, name in enumerate(component_names):
            components.attrs[str(index)] = name
        components.create_dataset("meshes", data=mesh_types)
        components.create_dataset("configuration", data=configuration)

        port_group = stream.create_group("ports")
        port_group.attrs["number_of_references"] = len(ports)
        for (name1, name2, attr1, attr2), index in ports.items():
            port_name = f"port{index}"
            port_group.attrs[str(index)] = port_name
            port = port_group.create_group(port_name)
            port.attrs["comp1"] = name1
            port.attrs["comp2"] = name2
            port.attrs["attr1"] = attr1
            port.attrs["attr2"] = attr2
            comp2_configuration = np.zeros(6)
            comp2_configuration[:2] = {
                1: (0, -1), 2: (1, 0), 3: (0, 1), 4: (-1, 0)
            }[attr1]
            port.create_dataset("comp2_configuration", data=comp2_configuration)
        port_group.create_dataset("interface", data=interface_data)
        stream.create_dataset("boundary", data=boundary_data)

    print(f"{count} blocks: {output}")


parser = argparse.ArgumentParser(
    description="Generate structured NxM benchmark configurations in linear time."
)
parser.add_argument("sizes", nargs="+", type=dimensions)
parser.add_argument("--seed", type=int, default=1)
parser.add_argument("--output-dir", default="configs/benchmarks")
args = parser.parse_args()

for shape in args.sizes:
    nx, ny = shape
    generate(
        nx,
        ny,
        args.seed,
        Path(args.output_dir) / f"box-channel.{nx}x{ny}.h5",
    )
