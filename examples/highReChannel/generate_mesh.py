#!/usr/bin/env python3
# Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: MIT
"""Structured quadrilateral mesh for the OpenLB ``bump2d`` channel.

The fluid domain is the channel ``[0, 25] x [0, 1.5]`` with the rectangular step
``[5, 10] x [0, 0.75]`` removed, matching ``openlb_examples/bump2d_*.cpp``::

    y=1.5 +-------------------------------------------------+
          |  B        |  C        |  E                      |
    y=.75 +-----------+-----------+-------------------------+
          |  A        |###########|  D                      |
    y=0   +-----------+###########+-------------------------+
         x=0         x=5        x=10                      x=25

Five rectangular blocks share whole curves, so transfinite meshing is conforming
by construction and there is no divisibility constraint on the element size.

Boundary physical tags follow the ``backward_facing_step`` parameterized problem
(see ``src/parameterized_problem.cpp``) and ``utils/gmsh/bfs.geo``:

    1 -- dirichlet inflow    (x = 0)
    2 -- no-slip walls       (channel floor/ceiling and all three step faces)
    3 -- neumann outflow     (x = 25)

All five blocks share a single surface tag so that the mesh loads as one
subdomain (``pmesh->attributes.Max() == 1``).
"""

import argparse
import os
import sys

# Channel geometry, in the same physical units as openlb_examples/bump2d_*.cpp.
LX0 = 25.0   # channel length
LY0 = 1.5    # channel height
X0_STEP = 5.0
LX1 = 5.0    # step length
LY1 = 0.75   # step height

XS = (0.0, X0_STEP, X0_STEP + LX1, LX0)   # 0, 5, 10, 25
YS = (0.0, LY1, LY0)                      # 0, 0.75, 1.5

# Inflow peak velocity and the gap-averaged over-speed used for the suggested
# timestep: the channel narrows 1.5 -> 0.75 over the step, so the mean doubles.
U_PEAK = 1.5
U_MAX = 2.0 * U_PEAK
CFL_TARGET = 0.2
VEL_ORDER = 2   # discretization/order: 1 -> uorder = 2 (Q2 velocity)


def fmt(value):
    """Compact fixed-point format, so 0.0625 -> '0.0625' rather than '0.062500'."""
    return ("%.10f" % value).rstrip("0").rstrip(".")


def divisions(length, h):
    return max(1, int(round(length / h)))


def build(h, out):
    import gmsh

    nx = [divisions(XS[i + 1] - XS[i], h) for i in range(3)]
    ny = [divisions(YS[i + 1] - YS[i], h) for i in range(2)]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("bump2d")
    geo = gmsh.model.geo

    # Points, indexed by (column, row) of the 4 x 3 lattice of corners. The two
    # corners interior to the step, (5, 0) and (10, 0), are still needed: they
    # are the bottom of the step's left and right faces.
    p = {}
    for i, x in enumerate(XS):
        for j, y in enumerate(YS):
            p[(i, j)] = geo.addPoint(x, y, 0.0, h)

    # Horizontal curves, left to right.
    h_bot_l = geo.addLine(p[(0, 0)], p[(1, 0)])   # floor, x in [0, 5]
    h_bot_r = geo.addLine(p[(2, 0)], p[(3, 0)])   # floor, x in [10, 25]
    h_mid_l = geo.addLine(p[(0, 1)], p[(1, 1)])   # interior A|B
    h_step = geo.addLine(p[(1, 1)], p[(2, 1)])    # step top (wall)
    h_mid_r = geo.addLine(p[(2, 1)], p[(3, 1)])   # interior D|E
    h_top_l = geo.addLine(p[(0, 2)], p[(1, 2)])   # ceiling, x in [0, 5]
    h_top_c = geo.addLine(p[(1, 2)], p[(2, 2)])   # ceiling, x in [5, 10]
    h_top_r = geo.addLine(p[(2, 2)], p[(3, 2)])   # ceiling, x in [10, 25]

    # Vertical curves, bottom to top.
    v_in_b = geo.addLine(p[(0, 0)], p[(0, 1)])    # inflow, lower half
    v_in_t = geo.addLine(p[(0, 1)], p[(0, 2)])    # inflow, upper half
    v_step_l = geo.addLine(p[(1, 0)], p[(1, 1)])  # step left face (wall)
    v_mid_l = geo.addLine(p[(1, 1)], p[(1, 2)])   # interior B|C
    v_step_r = geo.addLine(p[(2, 0)], p[(2, 1)])  # step right face (wall)
    v_mid_r = geo.addLine(p[(2, 1)], p[(2, 2)])   # interior C|E
    v_out_b = geo.addLine(p[(3, 0)], p[(3, 1)])   # outflow, lower half
    v_out_t = geo.addLine(p[(3, 1)], p[(3, 2)])   # outflow, upper half

    # Blocks, counterclockwise. Each shares whole curves with its neighbours.
    blocks = {
        "A": ([h_bot_l, v_step_l, -h_mid_l, -v_in_b],
              [p[(0, 0)], p[(1, 0)], p[(1, 1)], p[(0, 1)]], nx[0] * ny[0]),
        "B": ([h_mid_l, v_mid_l, -h_top_l, -v_in_t],
              [p[(0, 1)], p[(1, 1)], p[(1, 2)], p[(0, 2)]], nx[0] * ny[1]),
        "C": ([h_step, v_mid_r, -h_top_c, -v_mid_l],
              [p[(1, 1)], p[(2, 1)], p[(2, 2)], p[(1, 2)]], nx[1] * ny[1]),
        "D": ([h_bot_r, v_out_b, -h_mid_r, -v_step_r],
              [p[(2, 0)], p[(3, 0)], p[(3, 1)], p[(2, 1)]], nx[2] * ny[0]),
        "E": ([h_mid_r, v_out_t, -h_top_r, -v_mid_r],
              [p[(2, 1)], p[(3, 1)], p[(3, 2)], p[(2, 2)]], nx[2] * ny[1]),
    }

    surfaces = {}
    for name, (loop, corners, _) in blocks.items():
        surfaces[name] = geo.addPlaneSurface([geo.addCurveLoop(loop)])
        geo.mesh.setTransfiniteSurface(surfaces[name], cornerTags=corners)
        geo.mesh.setRecombine(2, surfaces[name])

    # Node counts = elements + 1. Opposite sides of every block agree because
    # each block spans exactly one (nx, ny) pair.
    for curve, n in ((h_bot_l, nx[0]), (h_mid_l, nx[0]), (h_top_l, nx[0]),
                     (h_step, nx[1]), (h_top_c, nx[1]),
                     (h_bot_r, nx[2]), (h_mid_r, nx[2]), (h_top_r, nx[2]),
                     (v_in_b, ny[0]), (v_step_l, ny[0]),
                     (v_step_r, ny[0]), (v_out_b, ny[0]),
                     (v_in_t, ny[1]), (v_mid_l, ny[1]),
                     (v_mid_r, ny[1]), (v_out_t, ny[1])):
        geo.mesh.setTransfiniteCurve(curve, n + 1)

    geo.synchronize()

    # Interior curves (h_mid_l, h_mid_r, v_mid_l, v_mid_r) are deliberately left
    # out of every physical group, so gmsh does not emit them as boundary
    # elements and MFEM sees exactly three boundary attributes.
    gmsh.model.addPhysicalGroup(1, [v_in_b, v_in_t], 1)
    gmsh.model.addPhysicalGroup(1, [h_bot_l, h_bot_r, h_step, v_step_l, v_step_r,
                                    h_top_l, h_top_c, h_top_r], 2)
    gmsh.model.addPhysicalGroup(1, [v_out_b, v_out_t], 3)
    gmsh.model.addPhysicalGroup(2, sorted(surfaces.values()), 1)

    # MFEM's gmsh reader and the rest of this repo expect format 2.2.
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.model.mesh.generate(2)

    parent = os.path.dirname(os.path.abspath(out))
    if parent:
        os.makedirs(parent, exist_ok=True)
    gmsh.write(out)
    gmsh.finalize()

    counts = {name: blocks[name][2] for name in blocks}
    return nx, ny, counts


def report(h, out, nx, ny, counts):
    dx = [(XS[i + 1] - XS[i]) / nx[i] for i in range(3)]
    dy = [(YS[j + 1] - YS[j]) / ny[j] for j in range(2)]
    h_eff = h / VEL_ORDER
    dt = CFL_TARGET * h_eff / U_MAX

    print("wrote %s" % out)
    print("  requested h    : %s" % fmt(h))
    print("  realized dx    : %s (x in [0,5]), %s ([5,10]), %s ([10,25])"
          % tuple(fmt(v) for v in dx))
    print("  realized dy    : %s (y in [0,.75]), %s ([.75,1.5])"
          % tuple(fmt(v) for v in dy))
    print("  elements       : %s (A=%d B=%d C=%d D=%d E=%d)"
          % (sum(counts.values()), counts["A"], counts["B"], counts["C"],
             counts["D"], counts["E"]))
    print("  velocity-node spacing h_eff = h/%d : %s" % (VEL_ORDER, fmt(h_eff)))
    print("  suggested timestep_size for CFL=%s at u_max=%s : %.3e"
          % (fmt(CFL_TARGET), fmt(U_MAX), dt))


def read_msh22(path):
    """Minimal MSH 2.2 reader: returns (nodes, elements).

    nodes:    {id: (x, y)}
    elements: [(type, physical_tag, [node ids])]
    """
    nodes, elements = {}, []
    with open(path) as f:
        lines = iter(f.read().splitlines())
    for line in lines:
        if line == "$Nodes":
            for _ in range(int(next(lines))):
                parts = next(lines).split()
                nodes[int(parts[0])] = (float(parts[1]), float(parts[2]))
        elif line == "$Elements":
            for _ in range(int(next(lines))):
                parts = [int(v) for v in next(lines).split()]
                etype, ntags = parts[1], parts[2]
                elements.append((etype, parts[3], parts[3 + ntags:]))
    return nodes, elements


def check(path, expected_elements):
    """Structural assertions on the written mesh. No MFEM needed."""
    nodes, elements = read_msh22(path)
    failures = []

    def require(cond, msg):
        if not cond:
            failures.append(msg)

    quads = [e for e in elements if e[0] == 3]
    lines = [e for e in elements if e[0] == 1]
    others = [e for e in elements if e[0] not in (1, 3)]

    require(not others,
            "found %d non-quad/non-line elements (types %s); SteadyNSSolver "
            "requires a single element type"
            % (len(others), sorted({e[0] for e in others})))
    require(len(quads) == expected_elements,
            "expected %d quads, found %d" % (expected_elements, len(quads)))

    # Every node is used: MFEM's gmsh reader does not tolerate orphans.
    used = {n for e in elements for n in e[2]}
    require(used == set(nodes),
            "%d nodes are not referenced by any element" % len(set(nodes) - used))

    # Conformity: every edge is shared by exactly 1 (boundary) or 2 (interior) quads.
    edge_count = {}
    for _, _, ns in quads:
        for k in range(4):
            e = tuple(sorted((ns[k], ns[(k + 1) % 4])))
            edge_count[e] = edge_count.get(e, 0) + 1
    require(all(c in (1, 2) for c in edge_count.values()),
            "non-conforming mesh: %d edges are shared by more than 2 quads"
            % sum(1 for c in edge_count.values() if c > 2))

    # Exactly one surface tag -> one subdomain.
    surf_tags = {e[1] for e in quads}
    require(surf_tags == {1}, "expected surface tag {1}, found %s" % surf_tags)

    # Exactly three boundary tags, each on the right geometry.
    bdr_tags = {e[1] for e in lines}
    require(bdr_tags == {1, 2, 3},
            "expected boundary tags {1,2,3}, found %s -- an interior curve "
            "probably leaked into a physical group" % bdr_tags)

    def on(tag, predicate):
        bad = [e for e in lines if e[1] == tag
               if not all(predicate(*nodes[n]) for n in e[2])]
        return len(bad)

    eps = 1e-9
    require(on(1, lambda x, y: abs(x - 0.0) < eps) == 0,
            "tag-1 (inflow) edges are not all on x = 0")
    require(on(3, lambda x, y: abs(x - LX0) < eps) == 0,
            "tag-3 (outflow) edges are not all on x = 25")

    def is_wall(x, y):
        floor = abs(y) < eps and (x <= XS[1] + eps or x >= XS[2] - eps)
        ceiling = abs(y - LY0) < eps
        step_top = abs(y - LY1) < eps and XS[1] - eps <= x <= XS[2] + eps
        step_side = ((abs(x - XS[1]) < eps or abs(x - XS[2]) < eps)
                     and -eps <= y <= LY1 + eps)
        return floor or ceiling or step_top or step_side

    require(on(2, is_wall) == 0,
            "tag-2 (wall) edges include segments that are not on the channel "
            "floor/ceiling or a step face")

    # The step is actually a hole: no element centroid inside it.
    inside = 0
    for _, _, ns in quads:
        cx = sum(nodes[n][0] for n in ns) / 4.0
        cy = sum(nodes[n][1] for n in ns) / 4.0
        if XS[1] < cx < XS[2] and YS[0] < cy < YS[1]:
            inside += 1
    require(inside == 0, "%d elements lie inside the step [5,10]x[0,0.75]" % inside)

    if failures:
        for msg in failures:
            print("  FAIL: %s" % msg)
        return False

    print("  check          : OK (%d quads, %d boundary edges, tags %s)"
          % (len(quads), len(lines), sorted(bdr_tags)))
    return True


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h", type=float, default=0.0625,
                        help="target element size (default: 0.0625)")
    parser.add_argument("--out", default=None,
                        help="output .msh path (default: meshes/bump2d.h<h>.msh)")
    parser.add_argument("--check", action="store_true",
                        help="re-read the written mesh and assert its invariants")
    args = parser.parse_args()

    if args.h <= 0.0:
        parser.error("--h must be positive")

    out = args.out or os.path.join("meshes", "bump2d.h%s.msh" % fmt(args.h))

    try:
        nx, ny, counts = build(args.h, out)
    except ImportError:
        sys.exit("the gmsh python module is required: pip3 install gmsh")

    report(args.h, out, nx, ny, counts)

    if args.check and not check(out, sum(counts.values())):
        sys.exit(1)


if __name__ == "__main__":
    main()
