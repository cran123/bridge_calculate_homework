#!/usr/bin/env python3

import argparse
import csv
import math
import statistics


def read_displacements(path):
    rows = {}
    in_table = False
    with open(path) as f:
        for line in f:
            if "D I S P L A C E M E N T S" in line:
                in_table = True
                continue
            if not in_table:
                continue
            fields = line.split()
            if len(fields) >= 7 and fields[0].isdigit():
                rows[int(fields[0])] = tuple(float(value) for value in fields[1:7])
            elif rows and (not fields or not fields[0].isdigit()):
                break
    return rows


def rms(values):
    return math.sqrt(sum(value * value for value in values) / len(values)) if values else 0.0


def main():
    parser = argparse.ArgumentParser(description="Compare two STAP++ .out displacement tables.")
    parser.add_argument("baseline_out")
    parser.add_argument("candidate_out")
    parser.add_argument("--csv", required=True)
    args = parser.parse_args()

    baseline = read_displacements(args.baseline_out)
    candidate = read_displacements(args.candidate_out)
    node_ids = sorted(set(baseline) & set(candidate))
    if not node_ids:
        raise SystemExit("No matching displacement nodes found.")

    labels = ("UX", "UY", "UZ", "RX", "RY", "RZ")
    component_abs = [[] for _ in labels]
    vector_abs = []
    worst = None

    with open(args.csv, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["node"]
        header += [f"baseline_{label}" for label in labels]
        header += [f"candidate_{label}" for label in labels]
        header += [f"d{label}" for label in labels]
        header += ["disp_vector_error", "rot_vector_error"]
        writer.writerow(header)

        for node in node_ids:
            base = baseline[node]
            cand = candidate[node]
            delta = tuple(cand[i] - base[i] for i in range(6))
            disp_error = math.sqrt(sum(delta[i] * delta[i] for i in range(3)))
            rot_error = math.sqrt(sum(delta[i] * delta[i] for i in range(3, 6)))
            vector_abs.append(disp_error)
            for i, value in enumerate(delta):
                component_abs[i].append(abs(value))
            if worst is None or disp_error > worst[0]:
                worst = (disp_error, rot_error, node, base, cand, delta)
            writer.writerow([node] + list(base) + list(cand) + list(delta) + [disp_error, rot_error])

    print(f"Matched nodes: {len(node_ids)} / baseline {len(baseline)} / candidate {len(candidate)}")
    print(f"CSV written: {args.csv}")
    print(f"Displacement vector diff: max {max(vector_abs):.12g}, mean {statistics.mean(vector_abs):.12g}, rms {rms(vector_abs):.12g}")
    for label, values in zip(labels, component_abs):
        print(f"{label} abs diff: max {max(values):.12g}, mean {statistics.mean(values):.12g}, rms {rms(values):.12g}")
    print(
        "Worst displacement node: "
        f"{worst[2]}, baseline {worst[3]}, candidate {worst[4]}, delta {worst[5]}"
    )


if __name__ == "__main__":
    main()
