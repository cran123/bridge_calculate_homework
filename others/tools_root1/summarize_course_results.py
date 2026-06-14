#!/usr/bin/env python
"""Summarize STAP++ output files for the course-design report."""

from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


TIME_PATTERNS = {
    "input_s": re.compile(r"TIME FOR INPUT PHASE\s*=\s*([0-9.Ee+-]+)"),
    "stiffness_s": re.compile(r"TIME FOR CALCULATION OF STIFFNESS MATRIX\s*=\s*([0-9.Ee+-]+)"),
    "factorization_s": re.compile(r"TIME FOR FACTORIZATION\s*=\s*([0-9.Ee+-]+)"),
    "load_solution_s": re.compile(r"TIME FOR LOAD CASE SOLUTIONS\s*=\s*([0-9.Ee+-]+)"),
    "output_s": re.compile(r"TIME FOR RESULT OUTPUT\s*=\s*([0-9.Ee+-]+)"),
    "total_s": re.compile(r"T O T A L\s+S O L U T I O N\s+T I M E\s*=\s*([0-9.Ee+-]+)"),
}


def parse_dat(path: Path) -> dict[str, int | str]:
    lines = [line.strip() for line in path.read_text(errors="ignore").splitlines() if line.strip()]
    if len(lines) < 2:
        return {"title": path.stem, "nodes": 0, "element_groups": 0, "load_cases": 0, "elements": 0}
    title = lines[0]
    parts = lines[1].split()
    nodes, groups, load_cases = map(int, parts[:3])
    idx = 2 + nodes
    for _ in range(load_cases):
        idx += 1
        load_header = lines[idx].split()
        nloads = int(load_header[0])
        idx += 1 + nloads
        if len(load_header) >= 6:
            idx += int(load_header[5])
    if idx < len(lines):
        maybe_mpc = lines[idx].split()
        if len(maybe_mpc) == 1:
            nmpc = int(float(maybe_mpc[0]))
            idx += 1
            for _ in range(nmpc):
                nterms = int(lines[idx].split()[0])
                idx += 1 + nterms
    elements = 0
    for _ in range(groups):
        if idx >= len(lines):
            break
        etype, nume, nummat = map(int, lines[idx].split()[:3])
        elements += nume
        idx += 1 + nummat + nume
    return {
        "title": title,
        "nodes": nodes,
        "element_groups": groups,
        "load_cases": load_cases,
        "elements": elements,
    }


def parse_out(path: Path) -> dict[str, float]:
    text = path.read_text(errors="ignore")
    result: dict[str, float] = {}
    for name, pattern in TIME_PATTERNS.items():
        match = pattern.search(text)
        if match:
            result[name] = float(match.group(1))
    return result


def main() -> None:
    rows = []
    for dat in sorted((ROOT / "cases").glob("*/stappp/*.dat")):
        out = dat.with_suffix(".out")
        info = parse_dat(dat)
        timing = parse_out(out) if out.exists() else {}
        rows.append(
            {
                "case": dat.parents[1].name,
                "input": str(dat.relative_to(ROOT)),
                **info,
                **timing,
                "out_exists": out.exists(),
                "vtu_count": len(list(dat.parent.glob(dat.stem + "_lc*.vtu"))),
            }
        )
    target = ROOT / "cases" / "result_summary.csv"
    target.parent.mkdir(exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as stream:
        fieldnames = [
            "case",
            "input",
            "title",
            "nodes",
            "elements",
            "element_groups",
            "load_cases",
            "input_s",
            "stiffness_s",
            "factorization_s",
            "load_solution_s",
            "output_s",
            "total_s",
            "out_exists",
            "vtu_count",
        ]
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(target)


if __name__ == "__main__":
    main()
