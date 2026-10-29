#!/usr/bin/env python3
"""Print the timing and relative-error datasets from the 16x16 runs."""

from pathlib import Path

import h5py
import numpy as np


def scalar(handle, key):
    return float(np.asarray(handle[key][()]).reshape(-1)[0])


results = Path(__file__).resolve().parent / "results"
with h5py.File(results / "fom_timings.h5", "r") as fom:
    fom_assemble = scalar(fom, "fom_assemble")
    fom_solve = scalar(fom, "fom_solve")

with h5py.File(results / "rom_vs_fom.h5", "r") as rom:
    rom_assemble = scalar(rom, "rom_assemble")
    rom_solve = scalar(rom, "rom_solve")
    error = np.asarray(rom["rel_error"][()]).reshape(-1)

print(f"FOM assembly: {fom_assemble:.6g} s")
print(f"FOM solve:    {fom_solve:.6g} s")
print(f"ROM assembly: {rom_assemble:.6g} s")
print(f"ROM solve:    {rom_solve:.6g} s")
print(f"Velocity relative error: {error[0]:.6%}")
print(f"Pressure relative error: {error[1]:.6%}")

