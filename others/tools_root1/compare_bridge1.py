#!/usr/bin/env python3

import argparse
import csv
import math
import runpy
import statistics
from collections import defaultdict
from pathlib import Path


def read_stap_displacements(path):
    displacements = {}
    in_table = False
    with open(path) as f:
        for line in f:
            if "D I S P L A C E M E N T S" in line:
                in_table = True
                continue
            if not in_table:
                continue
            fields = line.split()
            if len(fields) >= 4 and fields[0].isdigit():
                displacements[int(fields[0])] = (
                    float(fields[1]),
                    float(fields[2]),
                    float(fields[3]),
                )
            elif displacements and (not fields or not fields[0].isdigit()):
                break
    return displacements


def read_abaqus_displacements(path):
    displacements = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            key = (row["instance"].upper(), int(row["node"]))
            displacements[key] = (float(row["U1"]), float(row["U2"]), float(row["U3"]))
    return displacements


def rms(values):
    return math.sqrt(sum(value * value for value in values) / len(values)) if values else 0.0


def main():
    parser = argparse.ArgumentParser(description="Compare Bridge-1 STAPpp .out with Abaqus node CSV.")
    parser.add_argument("--inp", default="data/data-3/Bridge-1.inp")
    parser.add_argument("--abaqus", default="data/data_ref/bridge1_abaqus_nodes.csv")
    parser.add_argument("stap_out")
    parser.add_argument("--top", type=int, default=12)
    args = parser.parse_args()

    converter = Path(__file__).resolve().parents[1] / "others" / "tools" / "abaqus_to_stappp.py"
    module = runpy.run_path(str(converter))
    parsed = module["parse_inp"](args.inp)
    nodes, node_map, _flat_elements = module["flatten_model"](parsed[0], parsed[1])
    inverse_node_map = {gid: (inst.upper(), local_id) for (inst, local_id), gid in node_map.items()}

    abaqus = read_abaqus_displacements(args.abaqus)
    stap = read_stap_displacements(args.stap_out)

    rows = []
    by_instance = defaultdict(list)
    component_errors = [[], [], []]
    for gid, key in inverse_node_map.items():
        if key not in abaqus or gid not in stap:
            continue
        delta = tuple(stap[gid][i] - abaqus[key][i] for i in range(3))
        vector_error = math.sqrt(sum(value * value for value in delta))
        row = (vector_error, delta, gid, key, abaqus[key], stap[gid], nodes[gid])
        rows.append(row)
        by_instance[key[0]].append(row)
        for i in range(3):
            component_errors[i].append(abs(delta[i]))

    if not rows:
        raise SystemExit("No matched nodes. Check input paths and node numbering.")

    worst = max(rows, key=lambda item: item[0])
    print(f"Matched nodes: {len(rows)} / Abaqus {len(abaqus)} / STAP {len(stap)}")
    print(
        "Vector error: "
        f"max {worst[0]:.12g}, mean {statistics.mean(row[0] for row in rows):.12g}, "
        f"rms {rms([row[0] for row in rows]):.12g}"
    )
    print(
        "Worst node: "
        f"{worst[3][0]}:{worst[3][1]} -> STAP {worst[2]}, "
        f"xyz {worst[6]}, Abaqus {worst[4]}, STAP {worst[5]}"
    )
    for name, values in zip(("U1", "U2", "U3"), component_errors):
        print(f"{name} abs error: max {max(values):.12g}, mean {statistics.mean(values):.12g}, rms {rms(values):.12g}")

    print("\nTop instances by mean vector error:")
    for inst, inst_rows in sorted(
        by_instance.items(),
        key=lambda item: statistics.mean(row[0] for row in item[1]),
        reverse=True,
    )[: args.top]:
        inst_worst = max(inst_rows, key=lambda item: item[0])
        values = [row[0] for row in inst_rows]
        print(
            f"{inst:22s} n {len(inst_rows):4d} "
            f"mean {statistics.mean(values):.12g} rms {rms(values):.12g} "
            f"max {inst_worst[0]:.12g} at {inst_worst[3][0]}:{inst_worst[3][1]} "
            f"STAP {inst_worst[2]}"
        )


if __name__ == "__main__":
    main()
