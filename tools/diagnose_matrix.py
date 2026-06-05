#!/usr/bin/env python3

import argparse
import math
import re
from collections import Counter, defaultdict

import numpy as np


DOF_NAMES = ("UX", "UY", "UZ", "RX", "RY", "RZ")


def read_dat(path):
    lines = open(path, "r", encoding="utf-8", errors="ignore").read().splitlines()
    nnode, ngroup, ncase, mode = map(int, lines[1].split()[:4])
    nodes = {}
    flags = {}
    idx = 2
    eq = 0
    node_dof_to_eq = {}
    eq_to_node_dof = {}
    for _ in range(nnode):
        parts = lines[idx].split()
        nid = int(parts[0])
        flags[nid] = [int(v) for v in parts[1:7]]
        nodes[nid] = tuple(float(v) for v in parts[7:10])
        for dof, flag in enumerate(flags[nid], start=1):
            if flag == 0:
                eq += 1
                node_dof_to_eq[(nid, dof)] = eq
                eq_to_node_dof[eq] = (nid, dof)
        idx += 1

    for _ in range(ncase):
        idx += 1
        nload = int(lines[idx].split()[0])
        idx += 1 + nload

    mpcs = []
    maybe = lines[idx].split()
    if len(maybe) == 1:
        nmpc = int(maybe[0])
        idx += 1
        for _ in range(nmpc):
            nterm = int(lines[idx])
            idx += 1
            terms = []
            for _ in range(nterm):
                nid, dof, coef = lines[idx].split()[:3]
                terms.append((int(nid), int(dof), float(coef)))
                idx += 1
            mpcs.append(terms)

    groups = []
    for _ in range(ngroup):
        etype, nelem, nmat = map(int, lines[idx].split()[:3])
        idx += 1
        mats = []
        for _ in range(nmat):
            mats.append([float(v) for v in lines[idx].split()[1:]])
            idx += 1
        elems = []
        for _ in range(nelem):
            vals = lines[idx].split()
            elems.append((int(vals[0]), vals[1:]))
            idx += 1
        groups.append((etype, mats, elems))

    return nodes, flags, node_dof_to_eq, eq_to_node_dof, mpcs, groups


def local_indices(conn, dofs, node_dof_to_eq):
    eqs = []
    for nid in conn:
        for dof in dofs:
            eqs.append(node_dof_to_eq.get((nid, dof), 0))
    return eqs


def add_complete_graph(adjacency, eqs):
    active = [eq for eq in eqs if eq > 0]
    for i, a in enumerate(active):
        s = adjacency[a]
        for b in active[:i]:
            s.add(b)
            adjacency[b].add(a)


def build_dof_adjacency(node_dof_to_eq, groups):
    adjacency = {eq: set() for eq in node_dof_to_eq.values()}
    for etype, _mats, elems in groups:
        for _eid, raw in elems:
            vals = [int(float(v)) for v in raw]
            if etype in (4, 8):
                conn, dofs = vals[:8], (1, 2, 3)
            elif etype == 6:
                conn, dofs = vals[:4], (1, 2, 3, 4, 5, 6)
            elif etype == 5:
                conn, dofs = vals[:2], (1, 2, 3, 4, 5, 6)
            elif etype == 1:
                conn, dofs = vals[:2], (1, 2, 3)
            else:
                continue
            add_complete_graph(adjacency, local_indices(conn, dofs, node_dof_to_eq))
    return adjacency


def build_mpc_relations(mpcs, node_dof_to_eq):
    slave_to_terms = {}
    slave_eqs = set()
    for terms in mpcs:
        active = []
        slave = None
        slave_coef = None
        for nid, dof, coef in terms:
            eq = node_dof_to_eq.get((nid, dof), 0)
            if eq:
                active.append((eq, coef, nid, dof))
                if abs(coef - 1.0) < 1e-12 and slave is None:
                    slave = eq
                    slave_coef = coef
        if slave is None:
            continue
        slave_eqs.add(slave)
        rel = []
        for eq, coef, _nid, _dof in active:
            if eq != slave:
                rel.append((eq, -coef / slave_coef))
        slave_to_terms[slave] = rel
    return slave_to_terms, slave_eqs


def expand_slave_relations(slave_to_terms, slave_eqs):
    cache = {}
    visiting = set()

    def expand(eq):
        if eq not in slave_eqs:
            return {eq: 1.0}
        if eq in cache:
            return cache[eq]
        if eq in visiting:
            raise RuntimeError(f"MPC cycle at eq {eq}")
        visiting.add(eq)
        merged = defaultdict(float)
        for master, coef in slave_to_terms[eq]:
            for base, w in expand(master).items():
                merged[base] += coef * w
        visiting.remove(eq)
        cache[eq] = {k: v for k, v in merged.items() if abs(v) > 1e-14}
        return cache[eq]

    for eq in list(slave_eqs):
        expand(eq)
    return cache


def build_reduced_adjacency(adjacency, slave_to_terms, slave_eqs):
    expanded = expand_slave_relations(slave_to_terms, slave_eqs)
    independent = sorted(eq for eq in adjacency if eq not in slave_eqs)
    reduced_index = {eq: i + 1 for i, eq in enumerate(independent)}
    reduced_adj = {reduced_index[eq]: set() for eq in independent}
    source_eqs = defaultdict(set)

    def map_eq(eq):
        if eq in slave_eqs:
            return [base for base in expanded[eq] if base in reduced_index]
        return [eq]

    for eq, neighs in adjacency.items():
        mapped_i = map_eq(eq)
        for mi in mapped_i:
            ri = reduced_index[mi]
            source_eqs[ri].add(eq)
        for other in neighs:
            mapped_j = map_eq(other)
            for mi in mapped_i:
                ri = reduced_index[mi]
                for mj in mapped_j:
                    rj = reduced_index[mj]
                    if ri != rj:
                        reduced_adj[ri].add(rj)
                        reduced_adj[rj].add(ri)
                    source_eqs[ri].add(eq)
                    source_eqs[rj].add(other)
    return reduced_adj, reduced_index, source_eqs


def greedy_pivots(adjacency, fixed_order=None, max_steps=None):
    remaining = set(adjacency)
    degree_hist = []
    fill_edges = 0
    order = fixed_order or sorted(adjacency)
    if fixed_order is None:
        order = []
        buckets = defaultdict(set)
        for n in remaining:
            buckets[len(adjacency[n] & remaining)].add(n)
        while remaining:
            pivot = min(remaining, key=lambda n: (len(adjacency[n] & remaining), n))
            order.append(pivot)
            nbrs = list(adjacency[pivot] & remaining)
            degree_hist.append((pivot, len(nbrs)))
            for i, a in enumerate(nbrs):
                for b in nbrs[:i]:
                    if b not in adjacency[a]:
                        adjacency[a].add(b)
                        adjacency[b].add(a)
                        fill_edges += 1
            remaining.remove(pivot)
            if max_steps and len(order) >= max_steps:
                break
    else:
        for pivot in order:
            if pivot not in remaining:
                continue
            nbrs = list(adjacency[pivot] & remaining)
            degree_hist.append((pivot, len(nbrs)))
            for i, a in enumerate(nbrs):
                for b in nbrs[:i]:
                    if b not in adjacency[a]:
                        adjacency[a].add(b)
                        adjacency[b].add(a)
                        fill_edges += 1
            remaining.remove(pivot)
            if max_steps and len(degree_hist) >= max_steps:
                break
    return degree_hist, fill_edges


def parse_out_warning(path):
    if not path:
        return []
    result = []
    for line in open(path, "r", encoding="utf-8", errors="ignore"):
        if "Warning" in line or "Error" in line or "SparseLU" in line or "LDLT" in line:
            result.append(line.strip())
    return result[-20:]


def dof_summary(eqs, eq_to_node_dof):
    c = Counter()
    nodes = Counter()
    examples = []
    for eq in eqs:
        nd = eq_to_node_dof.get(eq)
        if not nd:
            continue
        nid, dof = nd
        c[DOF_NAMES[dof - 1]] += 1
        nodes[nid] += 1
        if len(examples) < 20:
            examples.append((eq, nid, dof))
    return c, nodes, examples


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dat")
    parser.add_argument("--out")
    parser.add_argument("--greedy-steps", type=int, default=5000)
    args = parser.parse_args()

    nodes, flags, node_dof_to_eq, eq_to_node_dof, mpcs, groups = read_dat(args.dat)
    adjacency = build_dof_adjacency(node_dof_to_eq, groups)
    slave_to_terms, slave_eqs = build_mpc_relations(mpcs, node_dof_to_eq)
    reduced_adj, reduced_index, source_eqs = build_reduced_adjacency(adjacency, slave_to_terms, slave_eqs)

    degrees = Counter(len(v) for v in adjacency.values())
    reduced_degrees = Counter(len(v) for v in reduced_adj.values())
    isolated = [eq for eq, s in adjacency.items() if not s]
    reduced_isolated = [eq for eq, s in reduced_adj.items() if not s]
    top_reduced = sorted(reduced_adj, key=lambda eq: len(reduced_adj[eq]), reverse=True)[:20]
    source_counts = sorted(source_eqs, key=lambda eq: len(source_eqs[eq]), reverse=True)[:20]

    print(f"File: {args.dat}")
    print(f"Nodes: {len(nodes)}")
    print(f"Full active equations: {len(adjacency)}")
    print(f"MPC constraints: {len(mpcs)}")
    print(f"MPC slave equations eliminated: {len(slave_eqs)}")
    print(f"Reduced equations: {len(reduced_adj)}")
    print()
    print(f"Full isolated equations by graph: {len(isolated)}")
    print(f"Reduced isolated equations by graph: {len(reduced_isolated)}")
    print(f"Full degree min/max: {min(degrees) if degrees else 0}/{max(degrees) if degrees else 0}")
    print(f"Reduced degree min/max: {min(reduced_degrees) if reduced_degrees else 0}/{max(reduced_degrees) if reduced_degrees else 0}")
    print()

    if isolated:
        c, _nodes, examples = dof_summary(isolated, eq_to_node_dof)
        print(f"Isolated full DOF summary: {dict(c)}")
        for eq, nid, dof in examples:
            xyz = nodes[nid]
            print(f"  eq {eq} node {nid} {DOF_NAMES[dof-1]} xyz=({xyz[0]:.6g},{xyz[1]:.6g},{xyz[2]:.6g})")
        print()

    print("Top reduced equations by graph degree:")
    for req in top_reduced:
        source = sorted(source_eqs[req])
        c, node_counts, examples = dof_summary(source, eq_to_node_dof)
        print(f"  red_eq {req}: degree={len(reduced_adj[req])} source_eqs={len(source)} dofs={dict(c)}")
        for eq, nid, dof in examples[:3]:
            xyz = nodes[nid]
            print(f"    source eq {eq} node {nid} {DOF_NAMES[dof-1]} xyz=({xyz[0]:.6g},{xyz[1]:.6g},{xyz[2]:.6g})")
    print()

    print("Reduced equations fed by most original equations:")
    for req in source_counts:
        source = sorted(source_eqs[req])
        c, _node_counts, examples = dof_summary(source, eq_to_node_dof)
        print(f"  red_eq {req}: source_eqs={len(source)} degree={len(reduced_adj[req])} dofs={dict(c)}")
        for eq, nid, dof in examples[:3]:
            xyz = nodes[nid]
            print(f"    source eq {eq} node {nid} {DOF_NAMES[dof-1]} xyz=({xyz[0]:.6g},{xyz[1]:.6g},{xyz[2]:.6g})")
    print()

    if args.greedy_steps:
        graph_copy = {k: set(v) for k, v in reduced_adj.items()}
        hist, fill = greedy_pivots(graph_copy, max_steps=args.greedy_steps)
        if hist:
            max_seen = max(hist, key=lambda x: x[1])
            avg = sum(d for _p, d in hist) / len(hist)
            print(f"Greedy symbolic elimination sampled steps: {len(hist)}")
            print(f"  max pivot degree seen: red_eq {max_seen[0]} degree {max_seen[1]}")
            print(f"  avg pivot degree seen: {avg:.2f}")
            print(f"  fill edges introduced in sample: {fill}")
        print()

    warnings = parse_out_warning(args.out)
    if warnings:
        print("Warnings/errors from output:")
        for line in warnings:
            print(f"  {line}")


if __name__ == "__main__":
    main()
