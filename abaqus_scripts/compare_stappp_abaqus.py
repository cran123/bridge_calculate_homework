#!/usr/bin/env python
"""Compare compact Abaqus CSV files with STAP++ nodal displacement output."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


DISP_ROW = re.compile(
    r"^\s*(\d+)\s+([+-]?(?:\d+\.\d*|\d*\.\d+|\d+)[Ee][+-]?\d+)\s+"
    r"([+-]?(?:\d+\.\d*|\d*\.\d+|\d+)[Ee][+-]?\d+)\s+"
    r"([+-]?(?:\d+\.\d*|\d*\.\d+|\d+)[Ee][+-]?\d+)"
)


def read_stappp_displacements(path: Path) -> dict[int, tuple[float, float, float]]:
    values = {}
    in_displacements = False
    for line in path.read_text(errors="ignore").splitlines():
        if "D I S P L A C E M E N T S" in line:
            in_displacements = True
            continue
        if in_displacements and "S T R E S S" in line:
            break
        if not in_displacements:
            continue
        match = DISP_ROW.match(line)
        if match:
            values[int(match.group(1))] = tuple(float(match.group(i)) for i in range(2, 5))
    return values


def read_abaqus_nodes(path: Path) -> dict[int, tuple[float, float, float]]:
    values = {}
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            values[int(row["node"])] = (float(row["u1"]), float(row["u2"]), float(row["u3"]))
    return values


def compare(stap: dict[int, tuple[float, float, float]], aba: dict[int, tuple[float, float, float]]):
    rows = []
    for node in sorted(set(stap) & set(aba)):
        diffs = [abs(stap[node][i] - aba[node][i]) for i in range(3)]
        denom = max(max(abs(x) for x in aba[node]), 1.0e-14)
        rows.append((node, max(diffs), max(diffs) / denom))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", default=".")
    args = parser.parse_args()
    root = Path(args.repo_root).resolve()
    out_rows = []
    detail_rows = []
    for nodes_csv in root.glob("cases/*/abaqus/*_nodes.csv"):
        case = nodes_csv.parents[1].name
        stap_outs = list((root / "cases" / case / "stappp").glob("*.out"))
        if not stap_outs:
            continue
        rows = compare(read_stappp_displacements(stap_outs[0]), read_abaqus_nodes(nodes_csv))
        if not rows:
            continue
        max_abs = max(row[1] for row in rows)
        max_rel = max(row[2] for row in rows)
        max_node = max(rows, key=lambda row: row[1])[0]
        out_rows.append({"case": case, "nodes": len(rows), "max_node": max_node, "max_abs": max_abs, "max_rel": max_rel})
        stap = read_stappp_displacements(stap_outs[0])
        aba = read_abaqus_nodes(nodes_csv)
        for node in sorted(set(stap) & set(aba)):
            diffs = [stap[node][i] - aba[node][i] for i in range(3)]
            detail_rows.append(
                {
                    "case": case,
                    "node": node,
                    "stap_u1": stap[node][0],
                    "stap_u2": stap[node][1],
                    "stap_u3": stap[node][2],
                    "abaqus_u1": aba[node][0],
                    "abaqus_u2": aba[node][1],
                    "abaqus_u3": aba[node][2],
                    "du1": diffs[0],
                    "du2": diffs[1],
                    "du3": diffs[2],
                    "max_abs": max(abs(x) for x in diffs),
                }
            )

    target = root / "cases" / "abaqus_comparison.csv"
    with target.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["case", "nodes", "max_node", "max_abs", "max_rel"])
        writer.writeheader()
        for row in out_rows:
            writer.writerow(row)
    detail_target = root / "cases" / "abaqus_comparison_details.csv"
    with detail_target.open("w", newline="", encoding="utf-8") as stream:
        fieldnames = [
            "case",
            "node",
            "stap_u1",
            "stap_u2",
            "stap_u3",
            "abaqus_u1",
            "abaqus_u2",
            "abaqus_u3",
            "du1",
            "du2",
            "du3",
            "max_abs",
        ]
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in detail_rows:
            writer.writerow(row)
    print(target)


if __name__ == "__main__":
    main()
