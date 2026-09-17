#!/usr/bin/env python3
# Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: MIT
"""Wall-graded quadrilateral meshes for the OpenLB channel cases.

Element size is driven by a gmsh distance-to-wall field: elements are ``h_min``
at every wall and grow to ``h_max = ratio * h_min`` once ``layer`` away from one.
Because the field measures true Euclidean distance to the wall *curves*, the
step's vertical faces and the convex corners are refined exactly like the floor
and ceiling, and nothing that is not a wall is refined at all. The inflow and
outflow curves are simply absent from the field, so they attract no refinement.

    bump2d      channel [0,25] x [0,1.5] less the step [5,10] x [0,0.75]

        y=1.5 +-------------------------------------------------+
              |                                                 |
        y=.75 +-----------+###########+-------------------------+
              |           |###########|                         |
        y=0   +-----------+###########+-------------------------+
             x=0         x=5        x=10                      x=25

    corner2d    channel [0,11] x [0,5] less [0,1]x[0,4] and [6,11]x[1,5]

        y=5 +---+-----------+###############+
            |   |           |###############|
        y=4 +---+-----------+###############+
            |###|           |###############|
        y=1 |###+-----------+---------------+
            |###|                           |
        y=0 +---+---------------------------+
           x=0 x=1         x=6            x=11

``--h`` sets the *degree-of-freedom budget*, not the element size: the mesh is
tuned so its element count matches what a uniform mesh of size ``h`` would have
had, so grading redistributes resolution without changing the cost.

Boundary physical tags depend on which parameterized problem will consume the
mesh, because the two disagree about what tag 2 means (see ``--tags``):

    bfs   backward_facing_step:  1 = inflow, 2 = wall, 3 = outflow
    chan  channel_flow:          4 = inflow, 1 = wall, 2 = outflow
"""

import argparse
import math
import os
import sys

# Blocks are (x0, x1, y0, y1). They no longer drive the mesher -- the mesh is
# unstructured -- but they still define the fluid region, from which the outline,
# the area, the solid holes and the element-count target are all derived.
CASES = {
    "bump2d": {
        "blocks": [(0.0, 5.0, 0.0, 0.75), (0.0, 5.0, 0.75, 1.5),
                   (5.0, 10.0, 0.75, 1.5),
                   (10.0, 25.0, 0.0, 0.75), (10.0, 25.0, 0.75, 1.5)],
        "x_in": 0.0,
        "x_out": 25.0,
        # Distance over which the size ramps from h_min to h_max. Half the gap
        # height, so the size reaches h_max exactly at mid-gap. Tying this to h
        # instead would fail at these resolutions: 6*h_max exceeds the 0.75 gap
        # and the field would never reach h_max anywhere.
        "layer": 0.375,
        # The channel narrows 1.5 -> 0.75 over the step, so by continuity the
        # mean velocity doubles: peak 1.5 -> ~3.0 in the gap.
        "u_max": 3.0,
        # Points that must NOT be refined: they lie on lines a structured mesher
        # would have been forced to refine along their whole length.
        "coarse_probes": [(2.5, 0.75), (5.0, 1.15)],
    },
    "corner2d": {
        "blocks": [(0.0, 1.0, 4.0, 5.0), (1.0, 6.0, 4.0, 5.0),
                   (1.0, 6.0, 1.0, 4.0),
                   (1.0, 6.0, 0.0, 1.0), (6.0, 11.0, 0.0, 1.0)],
        "x_in": 0.0,
        "x_out": 11.0,
        "layer": 0.5,
        # No area contraction: inlet and outlet are both height 1, so the bulk
        # peak stays ~1.5. The margin covers local speed-up at the two convex
        # corners (1,4) and (6,1).
        "u_max": 2.5,
        "coarse_probes": [(3.5, 4.0), (3.5, 1.0)],
    },
}

TAG_SCHEMES = {
    "bfs": {"inflow": 1, "wall": 2, "outflow": 3},    # BackwardFacingStep
    "chan": {"inflow": 4, "wall": 1, "outflow": 2},   # ChannelFlow
}

SURFACE_TAG = 1
CFL_TARGET = 0.2
VEL_ORDER = 2        # discretization/order: 1 -> uorder = 2 (Q2 velocity)
COUNT_TOL = 0.03     # element count must land within 3% of the budget
MAX_TUNE = 6
EPS = 1e-9


def fmt(value):
    """Compact fixed-point format, so 0.0625 -> '0.0625' rather than '0.062500'."""
    return ("%.10f" % value).rstrip("0").rstrip(".")


def divisions(length, h):
    return max(1, int(round(length / h)))


def block_edges(block):
    """The four edges of a block as ((xa, ya), (xb, yb)) corner pairs, CCW."""
    x0, x1, y0, y1 = block
    return [((x0, y0), (x1, y0)), ((x1, y0), (x1, y1)),
            ((x1, y1), (x0, y1)), ((x0, y1), (x0, y0))]


def edge_key(a, b):
    return tuple(sorted((a, b)))


def target_count(case_name, h):
    """Element budget: what a uniform mesh of size h would have had."""
    return sum(divisions(x1 - x0, h) * divisions(y1 - y0, h)
               for x0, x1, y0, y1 in CASES[case_name]["blocks"])


def classify(case, a, b):
    """inflow / outflow / wall for a boundary segment, from its x location."""
    if abs(a[0] - case["x_in"]) < EPS and abs(b[0] - case["x_in"]) < EPS:
        return "inflow"
    if abs(a[0] - case["x_out"]) < EPS and abs(b[0] - case["x_out"]) < EPS:
        return "outflow"
    return "wall"


def outline(case_name):
    """Ordered [(point, role)] cycle around the fluid region.

    Boundary edges are those used by exactly one block. They are chained into a
    cycle, then consecutive collinear segments *of the same role* are merged, so
    a long wall becomes one curve rather than several. Merging never spans a
    role change, which keeps the physical groups exact.
    """
    case = CASES[case_name]
    use_count = {}
    for b in case["blocks"]:
        for a, c in block_edges(b):
            use_count[edge_key(a, c)] = use_count.get(edge_key(a, c), 0) + 1
    if max(use_count.values()) > 2:
        sys.exit("case %s: an edge is shared by more than two blocks" % case_name)

    # Adjacency over boundary edges; every boundary point has exactly two.
    adj = {}
    for (a, b), n in use_count.items():
        if n != 1:
            continue
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    if any(len(v) != 2 for v in adj.values()):
        sys.exit("case %s: fluid boundary is not a simple closed loop" % case_name)

    start = min(adj)
    cycle, prev, cur = [start], None, start
    while True:
        nxt = next(p for p in adj[cur] if p != prev)
        if nxt == start:
            break
        cycle.append(nxt)
        prev, cur = cur, nxt

    # Merge collinear runs that share a role.
    def role_of(i):
        return classify(case, cycle[i], cycle[(i + 1) % len(cycle)])

    merged = True
    while merged and len(cycle) > 3:
        merged = False
        for i in range(len(cycle)):
            p, q, r = cycle[i - 1], cycle[i], cycle[(i + 1) % len(cycle)]
            cross = (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
            if abs(cross) > EPS:
                continue
            if classify(case, p, q) != classify(case, q, r):
                continue
            del cycle[i]
            merged = True
            break

    return [(cycle[i], classify(case, cycle[i], cycle[(i + 1) % len(cycle)]))
            for i in range(len(cycle))]


def perimeter(case_name):
    """Analytic boundary length per role."""
    lengths = {"inflow": 0.0, "wall": 0.0, "outflow": 0.0}
    ring = outline(case_name)
    for i, (p, role) in enumerate(ring):
        q = ring[(i + 1) % len(ring)][0]
        lengths[role] += math.hypot(q[0] - p[0], q[1] - p[1])
    return lengths


def solid_blocks(case_name):
    """Axis-aligned holes: the bounding box less the fluid blocks."""
    case = CASES[case_name]
    xs = sorted({v for b in case["blocks"] for v in b[:2]})
    ys = sorted({v for b in case["blocks"] for v in b[2:]})
    holes = []
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            cx, cy = (xs[i] + xs[i + 1]) / 2, (ys[j] + ys[j + 1]) / 2
            if not any(x0 < cx < x1 and y0 < cy < y1
                       for x0, x1, y0, y1 in case["blocks"]):
                holes.append((xs[i], xs[i + 1], ys[j], ys[j + 1]))
    return holes


def mesh_once(case_name, tags_name, h_min, ratio, out=None):
    """Generate once at the given h_min. Returns the element count."""
    import gmsh

    case = CASES[case_name]
    tags = TAG_SCHEMES[tags_name]
    ring = outline(case_name)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(case_name)
    geo = gmsh.model.geo

    pts = [geo.addPoint(p[0], p[1], 0.0) for p, _ in ring]
    lines, groups = [], {"inflow": [], "wall": [], "outflow": []}
    for i, (_, role) in enumerate(ring):
        tag = geo.addLine(pts[i], pts[(i + 1) % len(pts)])
        lines.append(tag)
        groups[role].append(tag)

    surf = geo.addPlaneSurface([geo.addCurveLoop(lines)])
    geo.synchronize()

    for role, tag_list in groups.items():
        if not tag_list:
            sys.exit("case %s: no %s curves found" % (case_name, role))
        gmsh.model.addPhysicalGroup(1, sorted(tag_list), tags[role])
    gmsh.model.addPhysicalGroup(2, [surf], SURFACE_TAG)

    # Size = h_min at any wall, ramping to h_max once `layer` away. Only the
    # wall curves are listed, so inflow and outflow attract no refinement.
    longest = max(math.hypot(ring[(i + 1) % len(ring)][0][0] - p[0],
                             ring[(i + 1) % len(ring)][0][1] - p[1])
                  for i, (p, _) in enumerate(ring))
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", sorted(groups["wall"]))
    gmsh.model.mesh.field.setNumber(1, "Sampling",
                                    min(5000, max(100, int(2 * longest / h_min))))
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", h_min)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", h_min * ratio)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(2, "DistMax", case["layer"])
    gmsh.model.mesh.field.setAsBackgroundMesh(2)

    # Without these three, gmsh blends point- and curvature-derived sizes into
    # the field and the grading silently degrades.
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    # Frontal-Delaunay for quads + blossom recombination. Verified to give 100%
    # quads for both geometries; --check asserts it, because SteadyNSSolver and
    # the integration rules require a single element type.
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.model.mesh.generate(2)

    types, tags_, _ = gmsh.model.mesh.getElements(2, surf)
    count = sum(len(t) for t in tags_)

    if out:
        parent = os.path.dirname(os.path.abspath(out))
        if parent:
            os.makedirs(parent, exist_ok=True)
        gmsh.write(out)
    gmsh.finalize()
    return count


def build(case_name, tags_name, h, ratio, out, verbose=True):
    """Tune h_min so the element count matches the budget set by h, then write."""
    target = target_count(case_name, h)

    # Start from the uniform size scaled by the mean-spacing rule of thumb,
    # h_min ~ h * 2/(1+ratio), then correct: count ~ 1/h_min^2.
    h_min = h * 2.0 / (1.0 + ratio)
    history = []
    for it in range(MAX_TUNE):
        count = mesh_once(case_name, tags_name, h_min, ratio)
        history.append((h_min, count))
        err = count / target - 1.0
        if verbose:
            print("    tune %d: h_min=%.5f -> %d elements (%+.1f%%)"
                  % (it + 1, h_min, count, 100 * err))
        if abs(err) <= COUNT_TOL:
            break
        h_min *= math.sqrt(count / target)

    count = mesh_once(case_name, tags_name, h_min, ratio, out)
    return h_min, count, target, history


def report(case_name, tags_name, h, ratio, h_min, count, target, out):
    case = CASES[case_name]
    h_max = h_min * ratio
    lengths = perimeter(case_name)
    # The layer spans roughly this many elements, so the per-element growth rate
    # is ratio^(1/n) -- the number boundary-layer meshing practice specifies.
    n_layer = case["layer"] / (0.5 * (h_min + h_max))
    growth = ratio ** (1.0 / n_layer) if n_layer > 1 else float("inf")

    print("wrote %s" % out)
    print("  case / tags    : %s / %s (inflow=%d wall=%d outflow=%d)"
          % (case_name, tags_name, TAG_SCHEMES[tags_name]["inflow"],
             TAG_SCHEMES[tags_name]["wall"], TAG_SCHEMES[tags_name]["outflow"]))
    print("  DOF budget     : --h %s -> %d elements; got %d (%+.1f%%)"
          % (fmt(h), target, count, 100 * (count / target - 1.0)))
    print("  size field     : h_min=%.5f  h_max=%.5f  ratio=%s  layer=%s"
          % (h_min, h_max, fmt(ratio), fmt(case["layer"])))
    print("                   h_min/h_uniform=%.2f, ~%.1f elements across the "
          "layer, growth~%.3f/element" % (h_min / h, n_layer, growth))
    print("  perimeter      : %.4f (inflow %.4f, wall %.4f, outflow %.4f)"
          % (sum(lengths.values()), lengths["inflow"], lengths["wall"],
             lengths["outflow"]))
    print("  velocity-node spacing h_eff = h/%d : %.5f (wall) .. %.5f (core)"
          % (VEL_ORDER, h_min / VEL_ORDER, h_max / VEL_ORDER))
    print("  timestep_size for CFL=%s at u_max=%s : %.3e (was %.3e uniform)"
          % (fmt(CFL_TARGET), fmt(case["u_max"]),
             CFL_TARGET * (h_min / VEL_ORDER) / case["u_max"],
             CFL_TARGET * (h / VEL_ORDER) / case["u_max"]))


def read_msh22(path):
    """Minimal MSH 2.2 reader -> ({id: (x, y)}, [(type, physical_tag, [ids])])."""
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


def check(case_name, tags_name, h, ratio, h_min, path):
    """Structural and grading assertions on the written mesh. No MFEM needed."""
    case = CASES[case_name]
    tags = TAG_SCHEMES[tags_name]
    h_max = h_min * ratio
    target = target_count(case_name, h)
    nodes, elements = read_msh22(path)
    failures = []

    def require(cond, msg):
        if not cond:
            failures.append(msg)

    quads = [e for e in elements if e[0] == 3]
    lines = [e for e in elements if e[0] == 1]
    others = [e for e in elements if e[0] not in (1, 3)]

    # Load-bearing: SteadyNSSolver requires one element type, and ir_nl/ir_face
    # are built from a single GeomType.
    require(not others,
            "found %d elements that are neither quads nor boundary lines "
            "(types %s)" % (len(others), sorted({e[0] for e in others})))
    require(abs(len(quads) / target - 1.0) <= COUNT_TOL,
            "element count %d is more than %.0f%% off the budget %d"
            % (len(quads), 100 * COUNT_TOL, target))

    used = {n for e in elements for n in e[2]}
    require(used == set(nodes),
            "%d nodes are not referenced by any element" % len(set(nodes) - used))

    edge_count = {}
    for _, _, ns in quads:
        for k in range(4):
            e = tuple(sorted((ns[k], ns[(k + 1) % 4])))
            edge_count[e] = edge_count.get(e, 0) + 1
    require(all(c in (1, 2) for c in edge_count.values()),
            "non-conforming mesh: %d edges are shared by more than 2 quads"
            % sum(1 for c in edge_count.values() if c > 2))

    surf_tags = {e[1] for e in quads}
    require(surf_tags == {SURFACE_TAG},
            "expected surface tag {%d}, found %s" % (SURFACE_TAG, surf_tags))

    want_tags = {tags[r] for r in ("inflow", "wall", "outflow")}
    found_tags = {e[1] for e in lines}
    require(found_tags == want_tags,
            "expected boundary tags %s, found %s" % (sorted(want_tags), sorted(found_tags)))

    # Boundary *length* per tag, which stays exact once spacing varies.
    want_len = perimeter(case_name)
    for role in ("inflow", "wall", "outflow"):
        got = sum(math.hypot(nodes[e[2][1]][0] - nodes[e[2][0]][0],
                             nodes[e[2][1]][1] - nodes[e[2][0]][1])
                  for e in lines if e[1] == tags[role])
        require(abs(got - want_len[role]) < 1e-6,
                "tag %d (%s): boundary length %.6f != %.6f"
                % (tags[role], role, got, want_len[role]))

    def all_at_x(role, x):
        return all(abs(nodes[n][0] - x) < EPS
                   for e in lines if e[1] == tags[role] for n in e[2])

    require(all_at_x("inflow", case["x_in"]),
            "tag %d (inflow) edges are not all on x = %s"
            % (tags["inflow"], fmt(case["x_in"])))
    require(all_at_x("outflow", case["x_out"]),
            "tag %d (outflow) edges are not all on x = %s"
            % (tags["outflow"], fmt(case["x_out"])))

    # Geometry: quad centroid and a size proxy the field can be compared against.
    def corners(ns):
        return [nodes[n] for n in ns]

    def area(ns):
        p = corners(ns)
        return abs(sum(p[i][0] * p[(i + 1) % 4][1] - p[(i + 1) % 4][0] * p[i][1]
                       for i in range(4))) / 2.0

    sizes, centroids = [], []
    for _, _, ns in quads:
        sizes.append(math.sqrt(area(ns)))
        p = corners(ns)
        centroids.append((sum(q[0] for q in p) / 4.0, sum(q[1] for q in p) / 4.0))

    holes = solid_blocks(case_name)
    inside = sum(1 for cx, cy in centroids
                 if any(x0 < cx < x1 and y0 < cy < y1 for x0, x1, y0, y1 in holes))
    require(inside == 0, "%d elements lie inside a solid region" % inside)

    # --- grading ---------------------------------------------------------
    wall_nodes = {n for e in lines if e[1] == tags["wall"] for n in e[2]}
    wall_sizes = [s for (_, _, ns), s in zip(quads, sizes)
                  if wall_nodes.intersection(ns)]
    require(bool(wall_sizes), "no element touches a wall")
    wall_sizes.sort()
    wall_med = wall_sizes[len(wall_sizes) // 2]

    if ratio > 1.0 + EPS:
        # Direction: walls must be the FINE end. Transposing SizeMin/SizeMax is
        # otherwise silent.
        require(wall_med < 1.6 * h_min,
                "median wall-adjacent element is %.4f, expected ~h_min=%.4f -- "
                "grading may be inverted" % (wall_med, h_min))

        # Points on lines a structured mesher would have had to refine along
        # their whole length. They must be coarse.
        for px, py in case["coarse_probes"]:
            near = min(zip(centroids, sizes),
                       key=lambda cs: (cs[0][0] - px) ** 2 + (cs[0][1] - py) ** 2)
            require(near[1] > 0.55 * h_max,
                    "interior probe (%s, %s) has element size %.4f, expected "
                    "~h_max=%.4f -- something is refining a non-wall line"
                    % (fmt(px), fmt(py), near[1], h_max))
            require(near[1] > 1.5 * wall_med,
                    "interior probe (%s, %s) is no coarser than the wall "
                    "elements (%.4f vs %.4f)" % (fmt(px), fmt(py), near[1], wall_med))

    if failures:
        for msg in failures:
            print("  FAIL: %s" % msg)
        return False

    print("  check          : OK (%d quads, %d boundary edges, tags %s)"
          % (len(quads), len(lines), sorted(found_tags)))
    print("                   wall element ~%.4f (h_min %.4f), core ~%.4f "
          "(h_max %.4f), measured span %.4f..%.4f"
          % (wall_med, h_min, max(sizes), h_max, min(sizes), max(sizes)))
    return True


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", default="bump2d", choices=sorted(CASES),
                        help="geometry to mesh (default: bump2d)")
    parser.add_argument("--tags", default="bfs", choices=sorted(TAG_SCHEMES),
                        help="boundary tag scheme, matching the parameterized "
                             "problem that will read the mesh (default: bfs)")
    parser.add_argument("--h", type=float, default=0.0625,
                        help="DOF budget: the uniform element size whose element "
                             "count this mesh should match (default: 0.0625)")
    parser.add_argument("--ratio", type=float, default=3.0,
                        help="h_max/h_min, core to wall (default: 3.0; "
                             "1.0 gives a uniform mesh)")
    parser.add_argument("--out", default=None,
                        help="output .msh path "
                             "(default: meshes/<case>.h<h>.<tags>.msh)")
    parser.add_argument("--check", action="store_true",
                        help="re-read the written mesh and assert its invariants")
    args = parser.parse_args()

    if args.h <= 0.0:
        parser.error("--h must be positive")
    if args.ratio < 1.0:
        parser.error("--ratio must be >= 1")

    out = args.out or os.path.join(
        "meshes", "%s.h%s.%s.msh" % (args.case, fmt(args.h), args.tags))

    try:
        h_min, count, target, _ = build(args.case, args.tags, args.h,
                                        args.ratio, out)
    except ImportError:
        sys.exit("the gmsh python module is required: pip3 install gmsh")

    report(args.case, args.tags, args.h, args.ratio, h_min, count, target, out)

    if args.check and not check(args.case, args.tags, args.h, args.ratio,
                                h_min, out):
        sys.exit(1)


if __name__ == "__main__":
    main()
