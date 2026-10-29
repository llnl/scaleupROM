#!/usr/bin/env python3
"""Generate the component catalog, training meshes, and paper 16x16 mesh."""

import argparse
import os
import sys
from pathlib import Path


def find_utilities(start: Path) -> Path:
    configured_root = os.environ.get("SCALEUPROM_ROOT")
    if configured_root:
        candidate = Path(configured_root) / "utils" / "python"
        if (candidate / "box_channel_config.py").is_file():
            return candidate
    for parent in (start, *start.parents):
        candidate = parent / "utils" / "python"
        if (candidate / "box_channel_config.py").is_file():
            return candidate
    raise FileNotFoundError("cannot locate scaleupROM/utils/python")


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(find_utilities(HERE)))

from box_channel_config import BoxChannelConfig  # noqa: E402
from config import Empty, ObjectInSpace  # noqa: E402


def components():
    return [
        Empty(),
        ObjectInSpace("square-circle"),
        ObjectInSpace("square-square"),
        ObjectInSpace("square-triangle"),
        ObjectInSpace("square-star"),
    ]


def add_components(config):
    for component in components():
        config.addComponent(component)


def write_component_catalog(path: Path):
    catalog = BoxChannelConfig(1, 5)
    add_components(catalog)

    # Register every oriented component/face pairing. BuildROM uses this
    # compact configuration to construct all reusable component operators.
    for component1 in catalog.comps:
        for displacement, face1 in component1.face_map.items():
            opposite = tuple(-value for value in displacement)
            for component2 in catalog.comps:
                if opposite in component2.face_map:
                    catalog.appendRefPort(
                        component1.name,
                        component2.name,
                        face1,
                        component2.face_map[opposite],
                    )

    for component_index in range(len(catalog.comps)):
        catalog.addMesh(component_index, 0)
    catalog.close()
    catalog.save(str(path))


# Component ordering extracted from the 16x16 configuration used by the paper.
# 0=empty, 1=circle, 2=square, 3=triangle, 4=star.
PAPER_16X16_COMPONENTS = [
    1, 1, 3, 0, 2, 1, 4, 1, 2, 4, 0, 2, 4, 1, 4, 2,
    0, 4, 1, 2, 4, 0, 3, 3, 0, 1, 3, 0, 0, 1, 3, 1,
    0, 3, 4, 0, 1, 1, 1, 2, 1, 4, 0, 4, 1, 4, 4, 3,
    1, 4, 0, 0, 3, 3, 0, 3, 4, 0, 0, 4, 2, 1, 0, 2,
    0, 1, 4, 1, 0, 3, 3, 1, 3, 4, 4, 3, 4, 3, 1, 0,
    3, 3, 0, 2, 0, 2, 0, 0, 2, 1, 0, 2, 3, 3, 1, 1,
    3, 4, 2, 2, 4, 2, 0, 2, 3, 2, 0, 2, 0, 0, 3, 0,
    2, 4, 3, 1, 4, 1, 2, 3, 3, 1, 1, 4, 1, 1, 4, 3,
    3, 2, 3, 4, 4, 0, 0, 1, 2, 1, 2, 3, 1, 2, 1, 0,
    1, 1, 1, 0, 0, 2, 0, 4, 4, 4, 3, 2, 4, 3, 0, 2,
    0, 0, 1, 1, 0, 4, 0, 3, 0, 3, 3, 3, 4, 0, 2, 4,
    3, 4, 4, 3, 4, 4, 3, 1, 4, 2, 1, 0, 4, 3, 4, 0,
    1, 0, 3, 4, 2, 2, 2, 4, 4, 2, 3, 3, 1, 0, 2, 0,
    3, 1, 2, 0, 0, 4, 1, 0, 1, 0, 4, 1, 1, 1, 0, 4,
    1, 4, 0, 0, 1, 3, 3, 4, 2, 2, 3, 0, 0, 2, 1, 3,
    4, 3, 1, 4, 4, 1, 2, 0, 2, 4, 4, 1, 0, 3, 4, 3,
]


def write_paper_target(path: Path):
    target = BoxChannelConfig(16, 16)
    add_components(target)
    for component_index in PAPER_16X16_COMPONENTS:
        target.addMesh(component_index, 0)
    target.close()
    target.save(str(path))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE / "configs")
    args = parser.parse_args()

    output = args.output.resolve()
    training = output / "training"
    training.mkdir(parents=True, exist_ok=True)

    write_component_catalog(output / "component_catalog.h5")

    samples = BoxChannelConfig(2, 2)
    add_components(samples)
    samples.prefix = str(training / "box-channel.config")
    count = samples.GenerateAllConfigs(0)
    if count != 625:
        raise RuntimeError(f"expected 625 training configurations, generated {count}")

    write_paper_target(output / "paper.16x16.h5")
    print(f"Generated {count} training configurations and {output / 'paper.16x16.h5'}")


if __name__ == "__main__":
    main()
