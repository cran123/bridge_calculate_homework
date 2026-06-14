#!/usr/bin/env python3

import argparse
import math
from collections import Counter, defaultdict, deque


DOF_NAMES = ("UX", "UY", "UZ", "RX", "RY", "RZ")


def read_nonempty(lines, index):
    while index < len(lines) and not lines[index].strip():
        index += 1
    return index


def parse_dat(path):
    lines = open(path, "r", encoding="utf-8", errors="ignore").read().splitlines()
    title = lines[0] if lines else ""
    nnode, ngroup, ncase, mode = map(int, lines[1].split()[:4])
    nodes = {}
    bcodes = {}
    index = 2
    for _ in range(nnode):
        parts = lines[index].split()
        nid = int(parts[0])
        bcodes[nid] = [int(v) for v in parts[1:7]]
        nodes[nid] = tuple(float(v) for v in parts[7:10])
        index += 1

    load_cases = []
    for _ in range(ncase):
        index = read_nonempty(lines, index)
        lcase = int(lines[index].split()[0])
        index += 1
        header = lines[index].split()
        index += 1
        nload = int(header[0])
        load_cases.append((lcase, header))
        index += nload

    mpcs = []
    index = read_nonempty(lines, index)
    maybe = lines[index].split()
    if len(maybe) == 1:
        nmpc = int(maybe[0])
        index += 1
        for _ in range(nmpc):
            nterm = int(lines[index].split()[0])
            index += 1
            terms = []
            for _ in range(nterm):
                node, dof, coef = lines[index].split()[:3]
                terms.append((int(node), int(dof), float(coef)))
                index += 1
            mpcs.append(terms)

    groups = []
    for _ in range(ngroup):
        index = read_nonempty(lines, index)
        etype, nelem, nmat = map(int, lines[index].split()[:3])
        index += 1
        mats = lines[index:index + nmat]
        index += nmat
        elems = []
        for _ in range(nelem):
            parts = lines[index].split()
            eid = int(parts[0])
            vals = [int(float(v)) for v in parts[1:]]
            elems.append((eid, vals))
            index += 1
        groups.append((etype, mats, elems))

    return title, nodes, bcodes, load_cases, mpcs, groups


def active_equations(bcodes):
    eq_to_node_dof = {}
    eq = 0
    for nid, flags in bcodes.items():
        for dof, flag in enumerate(flags, start=1):
            if flag == 0:
                eq += 1
                eq_to_node_dof[eq] = (nid, dof)
    return eq_to_node_dof


def element_used_dofs(groups):
    used = defaultdict(int)
    element_nodes = defaultdict(int)
    by_type = Counter()
    for etype, _mats, elems in groups:
        by_type[etype] += len(elems)
        for _eid, vals in elems:
            if etype in (4, 8):
                conn = vals[:8]
                dofs = (1, 2, 3)
            elif etype == 6:
                conn = vals[:4]
                dofs = (1, 2, 3, 4, 5, 6)
            elif etype == 5:
                conn = vals[:2]
                dofs = (1, 2, 3, 4, 5, 6)
            elif etype == 1:
                conn = vals[:2]
                dofs = (1, 2, 3)
            else:
                conn = []
                dofs = ()
            for nid in conn:
                element_nodes[nid] += 1
                for dof in dofs:
                    used[(nid, dof)] += 1
    return used, element_nodes, by_type


def analyze_mpcs(mpcs, bcodes, node_dof_to_eq):
    slave = {}
    graph = defaultdict(list)
    constrained = Counter()
    repeated = []
    bad_terms = []
    for idx, terms in enumerate(mpcs, start=1):
        positives = []
        active_terms = []
        for nid, dof, coef in terms:
            eq = node_dof_to_eq.get((nid, dof), 0)
            if eq == 0:
                bad_terms.append((idx, nid, dof, coef, "inactive"))
                continue
            active_terms.append((nid, dof, coef, eq))
            constrained[(nid, dof)] += 1
            if abs(coef - 1.0) < 1.0e-12:
                positives.append((nid, dof, eq))
        if positives:
            s = positives[0]
            if (s[0], s[1]) in slave:
                repeated.append((s[0], s[1], slave[(s[0], s[1])], idx))
            slave[(s[0], s[1])] = idx
            for nid, dof, coef, _eq in active_terms:
                if (nid, dof) != (s[0], s[1]):
                    graph[(s[0], s[1])].append((nid, dof))
    return slave, graph, constrained, repeated, bad_terms


def find_chains_and_cycles(graph, slave):
    slave_set = set(slave)
    max_depth = 0
    chain_roots = []
    cycles = []
    for start in slave_set:
        stack = [(start, [start])]
        while stack:
            item, path = stack.pop()
            max_depth = max(max_depth, len(path) - 1)
            next_slaves = [n for n in graph.get(item, []) if n in slave_set]
            if next_slaves:
                chain_roots.append(start)
            for nxt in next_slaves:
                if nxt in path:
                    cycles.append(path + [nxt])
                else:
                    stack.append((nxt, path + [nxt]))
    return max_depth, len(set(chain_roots)), cycles[:10]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dat")
    args = parser.parse_args()

    title, nodes, bcodes, load_cases, mpcs, groups = parse_dat(args.dat)
    eq_to_node_dof = active_equations(bcodes)
    node_dof_to_eq = {v: k for k, v in eq_to_node_dof.items()}
    used, element_nodes, by_type = element_used_dofs(groups)
    slave, graph, constrained, repeated, bad_terms = analyze_mpcs(mpcs, bcodes, node_dof_to_eq)
    max_depth, chain_roots, cycles = find_chains_and_cycles(graph, slave)

    active_dofs = set((nid, dof) for _eq, (nid, dof) in eq_to_node_dof.items())
    element_dofs = set(k for k, v in used.items() if v > 0)
    slave_dofs = set(slave)

    no_element = sorted(active_dofs - element_dofs)
    no_element_not_slave = sorted((active_dofs - element_dofs) - slave_dofs)
    slave_with_no_element = sorted(slave_dofs - element_dofs)
    rotational_active = sorted(k for k in active_dofs if k[1] >= 4)
    rotational_no_element = sorted(k for k in rotational_active if k not in element_dofs)

    print(f"File: {args.dat}")
    print(f"Title: {title}")
    print(f"Nodes: {len(nodes)}")
    print(f"Active equations: {len(eq_to_node_dof)}")
    print(f"Element groups: {dict(by_type)}")
    print(f"MPC constraints: {len(mpcs)}")
    print(f"MPC slave DOFs: {len(slave_dofs)}")
    print(f"MPC terms on inactive DOFs: {len(bad_terms)}")
    print(f"Repeated MPC slave DOFs: {len(repeated)}")
    print(f"MPC chain roots: {chain_roots}")
    print(f"Max MPC dependency depth: {max_depth}")
    print(f"MPC cycles: {len(cycles)}")
    print()

    print(f"Active DOFs not used by any element: {len(no_element)}")
    print(f"Active DOFs not used by any element and not eliminated by MPC: {len(no_element_not_slave)}")
    print(f"MPC slave DOFs not used by any element: {len(slave_with_no_element)}")
    print(f"Active rotational DOFs: {len(rotational_active)}")
    print(f"Active rotational DOFs not used by any element: {len(rotational_no_element)}")
    print()

    def show(label, items, limit=20):
        print(label)
        for nid, dof in items[:limit]:
            xyz = nodes.get(nid, (math.nan, math.nan, math.nan))
            print(f"  node {nid:6d} {DOF_NAMES[dof-1]} xyz=({xyz[0]:.6g}, {xyz[1]:.6g}, {xyz[2]:.6g})")
        if len(items) > limit:
            print(f"  ... {len(items) - limit} more")

    show("Examples: active DOFs not used by element and not MPC-slave", no_element_not_slave)
    show("Examples: MPC slave DOFs not used by element", slave_with_no_element)
    show("Examples: rotational active DOFs not used by element", rotational_no_element)

    if bad_terms[:10]:
        print("Bad MPC term examples:")
        for item in bad_terms[:10]:
            print(" ", item)
    if repeated[:10]:
        print("Repeated slave examples:")
        for item in repeated[:10]:
            print(" ", item)
    if cycles:
        print("Cycle examples:")
        for cycle in cycles:
            print(" ", " -> ".join(f"{n}:{DOF_NAMES[d-1]}" for n, d in cycle))


if __name__ == "__main__":
    main()
