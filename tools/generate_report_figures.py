"""Generate report figures from existing STAP++/Abaqus results.

This script deliberately separates figures backed by existing data from missing
items. It writes plots and a coverage matrix under cases/report_figures.
"""

from __future__ import annotations

import math
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "cases" / "report_figures"
CONV_WORK = OUT / "convergence" / "generated"

ELEMENT_CASES = {
    "bar": {
        "label": "Bar",
        "cases": ["bar_patch"],
        "patch_glob": "01_bar_type*_ux.png",
        "vtu": ROOT / "cases/bar_patch/stappp/bar_patch_lc1.vtu",
    },
    "beam": {
        "label": "Beam3D2",
        "cases": ["beam_cantilever"],
        "patch_glob": "03_beam_type*_ux.png",
        "vtu": ROOT / "cases/beam_cantilever/stappp/beam_cantilever_tip_load_lc1.vtu",
    },
    "plate4": {
        "label": "Plate4",
        "cases": ["plate4_membrane"],
        "patch_glob": "02_plate_type*_ux.png",
        "vtu": ROOT / "cases/plate4_membrane/stappp/plate4_membrane_tension_lc1.vtu",
    },
    "hex8": {
        "label": "Hex8",
        "cases": ["hex8_patch"],
        "patch_glob": "04_hex8_type*_ux.png",
        "vtu": ROOT / "cases/hex8_patch/stappp/hex8_uniaxial_patch_lc1.vtu",
    },
    "bbarhex8": {
        "label": "BbarHex8",
        "cases": ["bbarhex8_locking"],
        "patch_glob": "05_bbarhex8_type*_ux.png",
        "vtu": ROOT / "cases/bbarhex8_locking/stappp/bbarhex8_cantilever_lc1.vtu",
    },
}

CASE_LABELS = {
    case: meta["label"]
    for meta in ELEMENT_CASES.values()
    for case in meta["cases"]
}
CASE_LABELS["mixed_beam_plate_hex"] = "Mixed Beam/Plate/Hex"

DISP_RE = re.compile(
    r"^\s*(\d+)\s+([+-]?(?:\d+\.\d*|\d*\.\d+|\d+)[Ee][+-]?\d+)\s+"
    r"([+-]?(?:\d+\.\d*|\d*\.\d+|\d+)[Ee][+-]?\d+)\s+"
    r"([+-]?(?:\d+\.\d*|\d*\.\d+|\d+)[Ee][+-]?\d+)"
)


def ensure_dirs() -> None:
    for subdir in ["abaqus", "patch", "validation", "convergence", "coverage"]:
        (OUT / subdir).mkdir(parents=True, exist_ok=True)
    CONV_WORK.mkdir(parents=True, exist_ok=True)


def find_stappp() -> Path | None:
    candidates = [
        ROOT / "build_pardiso/Release/stap++.exe",
        ROOT / "build_eigen_run/Release/stap++.exe",
        ROOT / "build_pardiso_check/Release/stap++.exe",
        ROOT / "build/Debug/stap++.exe",
        ROOT / "build-merge/Release/stap++.exe",
    ]
    return next((path for path in candidates if path.exists()), None)


def write_lines(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def run_stappp(input_path: Path) -> Path | None:
    exe = find_stappp()
    if not exe:
        return None
    result = subprocess.run(
        [str(exe), str(input_path)],
        cwd=str(input_path.parent),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    if result.returncode != 0:
        return None
    out_path = input_path.with_suffix(".out")
    return out_path if out_path.exists() else None


def h8_node(ix: int, iy: int, iz: int, nx: int, ny: int, nz: int) -> int:
    return 1 + ix + (nx + 1) * (iy + (ny + 1) * iz)


def plate4_node(ix: int, iy: int, n: int) -> int:
    return 1 + ix + (n + 1) * iy


def savefig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()


def read_displacements(path: Path) -> pd.DataFrame:
    rows = []
    in_disp = False
    for line in path.read_text(errors="ignore").splitlines():
        if "D I S P L A C E M E N T S" in line:
            in_disp = True
            continue
        if in_disp and "S T R E S S" in line:
            break
        if not in_disp:
            continue
        match = DISP_RE.match(line)
        if match:
            ux, uy, uz = (float(match.group(i)) for i in range(2, 5))
            rows.append(
                {
                    "node": int(match.group(1)),
                    "ux": ux,
                    "uy": uy,
                    "uz": uz,
                    "umag": math.sqrt(ux * ux + uy * uy + uz * uz),
                }
            )
    return pd.DataFrame(rows)


def plot_abaqus_case(case: str, details: pd.DataFrame) -> Path:
    df = details[details["case"] == case].sort_values("node")
    element = next(
        (key for key, meta in ELEMENT_CASES.items() if case in meta["cases"]),
        case,
    )
    label = CASE_LABELS.get(case, case)
    x = np.arange(len(df))

    fig, axes = plt.subplots(2, 1, figsize=(9.2, 7.0), sharex=True)
    components = [("u1", "U1"), ("u2", "U2"), ("u3", "U3")]
    for suffix, label_comp in components:
        axes[0].plot(
            x,
            df[f"stap_{suffix}"],
            marker="o",
            linewidth=1.6,
            label=f"STAP++ {label_comp}",
        )
        axes[0].plot(
            x,
            df[f"abaqus_{suffix}"],
            marker="x",
            linestyle="--",
            linewidth=1.2,
            label=f"Abaqus {label_comp}",
        )
        axes[1].plot(
            x,
            df[f"d{suffix}"].abs(),
            marker="o",
            linewidth=1.5,
            label=f"|d{label_comp}|",
        )

    axes[0].set_title(f"{label}: STAP++ vs Abaqus nodal displacement")
    axes[0].set_ylabel("Displacement")
    axes[0].grid(True, alpha=0.28)
    axes[0].legend(ncol=3, fontsize=7.5)

    axes[1].set_yscale("symlog", linthresh=1.0e-12)
    axes[1].set_xlabel("Matched node index")
    axes[1].set_ylabel("Absolute component error")
    axes[1].grid(True, alpha=0.28)
    axes[1].legend(ncol=3, fontsize=8)

    path = OUT / "abaqus" / f"{element}_abaqus_displacement_comparison.png"
    savefig(path)
    return path


def plot_abaqus_summary(summary: pd.DataFrame) -> Path:
    labels = []
    max_abs = []
    max_rel = []
    for _, row in summary.iterrows():
        labels.append(CASE_LABELS.get(row["case"], row["case"]))
        max_abs.append(float(row["max_abs"]))
        max_rel.append(float(row["max_rel"]))

    x = np.arange(len(labels))
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.8), sharex=True)
    axes[0].bar(x, max_abs, color="#4477aa")
    axes[0].set_yscale("symlog", linthresh=1.0e-14)
    axes[0].set_ylabel("Max abs displacement error")
    axes[0].grid(True, axis="y", alpha=0.28)
    axes[1].bar(x, max_rel, color="#cc6677")
    axes[1].set_yscale("symlog", linthresh=1.0e-12)
    axes[1].set_ylabel("Max relative error")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=18, ha="right")
    axes[1].grid(True, axis="y", alpha=0.28)
    axes[0].set_title("Abaqus displacement-comparison error summary")
    path = OUT / "abaqus" / "abaqus_error_summary.png"
    savefig(path)
    return path


def copy_patch_images() -> list[Path]:
    rendered = []
    source_dir = ROOT / "cases" / "patch_test_paraview"
    for element, meta in ELEMENT_CASES.items():
        glob = meta["patch_glob"]
        if not glob:
            continue
        target_dir = OUT / "patch" / element
        if target_dir.exists():
            shutil.rmtree(target_dir)
        target_dir.mkdir(parents=True, exist_ok=True)
        for source in sorted(source_dir.glob(glob)):
            target = target_dir / source.name
            shutil.copy2(source, target)
            rendered.append(target)
    return rendered


def plot_validation_case(element: str, case: str) -> Path | None:
    out_path = ROOT / "cases" / case / "stappp" / f"{Path(case).name}.out"
    if case == "beam_cantilever":
        out_path = ROOT / "cases/beam_cantilever/stappp/beam_cantilever_tip_load.out"
    elif case == "plate4_membrane":
        out_path = ROOT / "cases/plate4_membrane/stappp/plate4_membrane_tension.out"
    elif case == "hex8_patch":
        out_path = ROOT / "cases/hex8_patch/stappp/hex8_uniaxial_patch.out"
    elif case == "bbarhex8_locking":
        out_path = ROOT / "cases/bbarhex8_locking/stappp/bbarhex8_cantilever.out"
    elif case == "bar_patch":
        out_path = ROOT / "cases/bar_patch/stappp/bar_patch.out"

    if not out_path.exists():
        return None
    df = read_displacements(out_path)
    if df.empty:
        return None
    x = np.arange(len(df))
    plt.figure(figsize=(8.6, 4.8))
    plt.plot(x, df["ux"], marker="o", label="UX")
    plt.plot(x, df["uy"], marker="s", label="UY")
    plt.plot(x, df["uz"], marker="^", label="UZ")
    plt.plot(x, df["umag"], marker="x", linestyle="--", label="|U|")
    plt.title(f"{ELEMENT_CASES[element]['label']}: validation nodal displacement")
    plt.xlabel("Node index")
    plt.ylabel("Displacement")
    plt.grid(True, alpha=0.28)
    plt.legend(ncol=4)
    path = OUT / "validation" / f"{element}_validation_displacement.png"
    savefig(path)
    return path


def parse_plate4_convergence() -> Path | None:
    n = 8
    generated = CONV_WORK / "plate4_cantilever_gravity_8x8.dat"
    if not generated.exists():
        write_plate4_gravity_dat(generated, n)
    if not generated.with_suffix(".out").exists():
        run_stappp(generated)

    rows = []
    out_paths = list((ROOT / "tests/convergence").glob("plate4_cantilever_gravity_*x*.out"))
    out_paths.extend(CONV_WORK.glob("plate4_cantilever_gravity_*x*.out"))
    for out_path in sorted(out_paths):
        match = re.search(r"_(\d+)x\1\.out$", out_path.name)
        if not match:
            continue
        n = int(match.group(1))
        df = read_displacements(out_path)
        if df.empty:
            continue
        tip = df.loc[df["uz"].abs().idxmax()]
        rows.append({"n": n, "h": 1.0 / n, "uz": float(tip["uz"])})
    if len(rows) < 2:
        return None
    conv = pd.DataFrame(rows).sort_values("n")
    reference = float(conv.iloc[-1]["uz"])
    conv["abs_error_vs_finest"] = (conv["uz"] - reference).abs()
    rate = estimate_convergence_rate(conv, "abs_error_vs_finest")
    conv["estimated_order"] = rate
    conv.to_csv(OUT / "convergence" / "plate4_existing_convergence.csv", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.2), sharex=True)
    axes[0].plot(conv["n"], conv["uz"].abs(), marker="o", label="max |UZ|")
    axes[0].set_ylabel("Displacement magnitude")
    axes[0].set_title("Plate4 existing convergence: cantilever gravity sequence")
    axes[0].grid(True, alpha=0.28)
    axes[0].legend()
    axes[1].plot(conv["n"], conv["abs_error_vs_finest"], marker="s", color="#cc6677")
    axes[1].set_yscale("symlog", linthresh=1.0e-14)
    axes[1].set_xlabel("Mesh divisions per side")
    axes[1].set_ylabel("Abs error vs finest")
    axes[1].grid(True, alpha=0.28)
    if not np.isnan(rate):
        axes[1].text(
            0.04,
            0.9,
            f"fitted order p = {rate:.3g}",
            transform=axes[1].transAxes,
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
        )
    path = OUT / "convergence" / "plate4_existing_convergence.png"
    savefig(path)
    return path


def write_plate4_gravity_dat(path: Path, n: int) -> None:
    lines = [f"Plate4 cantilever gravity convergence {n}x{n}", f"{(n + 1) * (n + 1)} 1 1 1"]
    for iy in range(n + 1):
        for ix in range(n + 1):
            node = plate4_node(ix, iy, n)
            if ix == 0:
                bc = "1 1 1 1 1 1"
            else:
                bc = "1 1 0 0 0 1"
            lines.append(f"{node} {bc} {ix / n:.12g} {iy / n:.12g} 0.0")
    lines.extend(["1", "0 1 0.0 0.0 -9.81", f"6 {n * n} 1", "1 3.0e10 0.2 0.1 2500.0"])
    eid = 1
    for iy in range(n):
        for ix in range(n):
            n1 = plate4_node(ix, iy, n)
            n2 = plate4_node(ix + 1, iy, n)
            n3 = plate4_node(ix + 1, iy + 1, n)
            n4 = plate4_node(ix, iy + 1, n)
            lines.append(f"{eid} {n1} {n2} {n3} {n4} 1")
            eid += 1
    write_lines(path, lines)


def generate_bar_convergence() -> Path | None:
    rows = []
    exact = 1.0 / (1.0e3 * math.pi)
    for n in [1, 2, 4, 8, 16, 32]:
        lines = [
            f"Bar sine-load convergence N{n}",
            f"{n + 1} 1 1 1",
        ]
        for i in range(n + 1):
            bc = 1 if i == 0 else 0
            lines.append(f"{i + 1} {bc} 1 1 1 1 1 {i / n:.12g} 0.0 0.0")
        h = 1.0 / n
        loads = []
        for i in range(1, n + 1):
            x = i / n
            weight = 0.5 if i == n else 1.0
            loads.append((i + 1, math.sin(math.pi * x) * h * weight))
        lines.extend(["1", str(len(loads))])
        lines.extend(f"{node} 1 {load:.12g}" for node, load in loads)
        lines.extend([f"1 {n} 1", "1 1.0E3 1.0 0.0"])
        for i in range(1, n + 1):
            lines.append(f"{i} {i} {i + 1} 1")
        dat = CONV_WORK / f"bar_convergence_n{n}.dat"
        write_lines(dat, lines)
        out = run_stappp(dat)
        if not out:
            continue
        df = read_displacements(out)
        if df.empty:
            continue
        value = float(df[df["node"] == n + 1]["ux"].iloc[0])
        rows.append(
            {
                "n": n,
                "h": 1.0 / n,
                "tip_ux": value,
                "relative_error": abs(value - exact) / exact,
            }
        )
    return plot_generated_convergence(
        "bar", "Bar sine axial-load convergence against analytical tip UX", rows, "tip_ux"
    )


def generate_beam_convergence() -> Path | None:
    rows = []
    exact = -1000.0 * 1.0**4 / (8.0 * 2.10e11 * 8.333333333e-6)
    for n in [1, 2, 4, 8, 16]:
        lines = [f"Beam distributed-load convergence N{n}", f"{n + 1} 1 1 1"]
        for i in range(n + 1):
            bc = "1 1 1 1 1 1" if i == 0 else "0 0 0 0 0 0"
            lines.append(f"{i + 1} {bc} {i / n:.12g} 0.0 0.0")
        loads = []
        element_load = -1000.0 / n
        for i in range(n + 1):
            share = 0.5 if i in (0, n) else 1.0
            if i == 0:
                continue
            loads.append((i + 1, element_load * share))
        lines.extend(
            [
                "1",
                str(len(loads)),
                *[f"{node} 3 {load:.12g}" for node, load in loads],
                f"5 {n} 1",
                "1 2.10E11 0.30 1.0E-2 8.333333333E-6 8.333333333E-6 1.666666667E-5 7850.0",
            ]
        )
        for i in range(1, n + 1):
            lines.append(f"{i} {i} {i + 1} 1 0.0 0.0 1.0")
        dat = CONV_WORK / f"beam_convergence_n{n}.dat"
        write_lines(dat, lines)
        out = run_stappp(dat)
        if not out:
            continue
        df = read_displacements(out)
        if df.empty:
            continue
        value = float(df[df["node"] == n + 1]["uz"].iloc[0])
        rows.append(
            {
                "n": n,
                "h": 1.0 / n,
                "tip_uz": value,
                "relative_error": abs(value - exact) / abs(exact),
            }
        )
    return plot_generated_convergence(
        "beam", "Beam3D2 distributed-load convergence against analytical tip UZ", rows, "tip_uz"
    )


def write_h8_convergence_dat(path: Path, n: int, etype: int) -> None:
    ny = nz = 2
    length = 4.0
    lines = [f"H8 convergence etype {etype} N{n}", f"{(n + 1) * (ny + 1) * (nz + 1)} 1 1 1"]
    for iz in range(nz + 1):
        for iy in range(ny + 1):
            for ix in range(n + 1):
                node = h8_node(ix, iy, iz, n, ny, nz)
                bc = "1 1 1 1 1 1" if ix == 0 else "0 0 0 1 1 1"
                lines.append(
                    f"{node} {bc} {length * ix / n:.12g} {iy / ny:.12g} {iz / nz:.12g}"
                )
    loaded = [
        h8_node(n, 0, 0, n, ny, nz),
        h8_node(n, ny, 0, n, ny, nz),
        h8_node(n, 0, nz, n, ny, nz),
        h8_node(n, ny, nz, n, ny, nz),
    ]
    lines.extend(["1", str(len(loaded))])
    for node in loaded:
        lines.append(f"{node} 3 -25.0")
    element_count = n * ny * nz
    lines.extend([f"{etype} {element_count} 1", "1 1.0E5 0.499 0.0"])
    eid = 1
    for ix in range(n):
        for iz in range(nz):
            for iy in range(ny):
                n1 = h8_node(ix, iy, iz, n, ny, nz)
                n2 = h8_node(ix + 1, iy, iz, n, ny, nz)
                n3 = h8_node(ix + 1, iy + 1, iz, n, ny, nz)
                n4 = h8_node(ix, iy + 1, iz, n, ny, nz)
                n5 = h8_node(ix, iy, iz + 1, n, ny, nz)
                n6 = h8_node(ix + 1, iy, iz + 1, n, ny, nz)
                n7 = h8_node(ix + 1, iy + 1, iz + 1, n, ny, nz)
                n8 = h8_node(ix, iy + 1, iz + 1, n, ny, nz)
                lines.append(f"{eid} {n1} {n2} {n3} {n4} {n5} {n6} {n7} {n8} 1")
                eid += 1
    write_lines(path, lines)


def generate_h8_convergence(element: str, etype: int, title: str) -> Path | None:
    rows = []
    for n in [1, 2, 4, 8]:
        dat = CONV_WORK / f"{element}_convergence_n{n}.dat"
        write_h8_convergence_dat(dat, n, etype)
        out = run_stappp(dat)
        if not out:
            continue
        df = read_displacements(out)
        if df.empty:
            continue
        value = float(df["uz"].min())
        rows.append({"n": n, "h": 1.0 / n, "max_downward_uz": value, "relative_error": np.nan})
    if len(rows) >= 2:
        reference = rows[-1]["max_downward_uz"]
        for row in rows:
            row["relative_error"] = abs(row["max_downward_uz"] - reference) / max(
                abs(reference), 1.0e-14
            )
    return plot_generated_convergence(element, title, rows, "max_downward_uz")


def plot_generated_convergence(
    element: str, title: str, rows: list[dict[str, float]], value_column: str
) -> Path | None:
    if len(rows) < 2:
        return None
    conv = pd.DataFrame(rows).sort_values("n")
    if "h" not in conv.columns:
        conv["h"] = 1.0 / conv["n"]
    rate = estimate_convergence_rate(conv, "relative_error")
    conv["estimated_order"] = rate
    csv_path = OUT / "convergence" / f"{element}_generated_convergence.csv"
    conv.to_csv(csv_path, index=False)
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.2), sharex=True)
    axes[0].plot(conv["n"], conv[value_column].abs(), marker="o")
    axes[0].set_ylabel(f"|{value_column}|")
    axes[0].grid(True, alpha=0.28)
    axes[0].set_title(title)
    axes[1].plot(conv["n"], conv["relative_error"], marker="s", color="#cc6677")
    axes[1].set_yscale("symlog", linthresh=1.0e-14)
    axes[1].set_xlabel("Mesh divisions / elements along length")
    axes[1].set_ylabel("Relative error")
    axes[1].grid(True, alpha=0.28)
    if not np.isnan(rate):
        axes[1].text(
            0.04,
            0.9,
            f"fitted order p = {rate:.3g}",
            transform=axes[1].transAxes,
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
        )
    path = OUT / "convergence" / f"{element}_generated_convergence.png"
    savefig(path)
    return path


def estimate_convergence_rate(conv: pd.DataFrame, error_column: str) -> float:
    if "h" in conv:
        h = conv["h"].astype(float)
    elif "n" in conv:
        h = 1.0 / conv["n"].astype(float)
    else:
        return float("nan")
    error = conv[error_column].astype(float).abs()
    mask = np.isfinite(h) & np.isfinite(error) & (h > 0.0) & (error > 0.0)
    if int(mask.sum()) < 2:
        return float("nan")
    coeff = np.polyfit(np.log(h[mask]), np.log(error[mask]), 1)
    return float(coeff[0])


def write_convergence_rate_summary() -> Path | None:
    specs = [
        ("Bar", OUT / "convergence" / "bar_generated_convergence.csv", "relative_error"),
        ("Beam3D2", OUT / "convergence" / "beam_generated_convergence.csv", "relative_error"),
        ("Plate4", OUT / "convergence" / "plate4_existing_convergence.csv", "abs_error_vs_finest"),
        ("Hex8", OUT / "convergence" / "hex8_generated_convergence.csv", "relative_error"),
        ("BbarHex8", OUT / "convergence" / "bbarhex8_generated_convergence.csv", "relative_error"),
    ]
    rows = []
    for element, path, error_column in specs:
        if not path.exists():
            continue
        conv = pd.read_csv(path)
        if error_column not in conv:
            continue
        if "h" not in conv and "n" in conv:
            conv["h"] = 1.0 / conv["n"]
        rate = estimate_convergence_rate(conv, error_column)
        usable = conv[
            np.isfinite(conv["h"].astype(float))
            & np.isfinite(conv[error_column].astype(float))
            & (conv["h"].astype(float) > 0.0)
            & (conv[error_column].astype(float).abs() > 0.0)
        ]
        rows.append(
            {
                "element": element,
                "error_column": error_column,
                "fit_points": len(usable),
                "fitted_order": rate,
                "source_csv": str(path.relative_to(ROOT)),
            }
        )
    if not rows:
        return None
    summary = pd.DataFrame(rows)
    path = OUT / "convergence" / "convergence_rate_summary.csv"
    summary.to_csv(path, index=False)
    return path


def read_vtu_arrays(path: Path) -> dict[str, np.ndarray]:
    root = ET.parse(path).getroot()
    piece = root.find(".//Piece")
    if piece is None:
        return {}
    arrays: dict[str, np.ndarray] = {}
    points_text = piece.find("./Points/DataArray").text or ""
    arrays["Points"] = np.array([float(x) for x in points_text.split()]).reshape((-1, 3))
    for section_name in ["PointData", "CellData"]:
        section = piece.find(f"./{section_name}")
        if section is None:
            continue
        for data in section.findall("./DataArray"):
            name = data.attrib.get("Name")
            if not name:
                continue
            text = data.text or ""
            values = np.array([float(x) for x in text.split()], dtype=float)
            components = int(data.attrib.get("NumberOfComponents", "1"))
            if components > 1:
                values = values.reshape((-1, components))
            arrays[name] = values
    return arrays


def generate_patch_numeric_summary() -> Path | None:
    vtu_dir = ROOT / "cases/patch_test_paraview/subdivided_vtu"
    patch_cases = [
        ("Bar", "bar_patch_6_segments", "SXX", 1.0),
        ("Beam3D2", "beam_patch_6_segments", "AxialForce", 1.0),
        ("Plate4", "plate4_patch_3x3", "SXX", 1.0e5),
        ("Hex8", "hex8_patch_2x2x2", "SXX", 1.0),
        ("BbarHex8", "bbarhex8_patch_2x2x2", "SXX", 1.0),
    ]
    trial_cases = [
        (
            "type1_all_node_displacement",
            "all_node_displacement_markers",
            "prescribed displacement on every node",
        ),
        (
            "type2_boundary_displacement",
            "boundary_displacement_markers",
            "prescribed displacement on boundary nodes",
        ),
        (
            "type3_equivalent_force_boundary",
            "equivalent_force_boundary_markers",
            "equivalent force boundary condition",
        ),
    ]
    rows = []
    for element, stem, stress_name, expected_stress in patch_cases:
        base_path = vtu_dir / f"{stem}.vtu"
        if not base_path.exists():
            continue
        arrays = read_vtu_arrays(base_path)
        points = arrays.get("Points")
        ux = arrays.get("UX")
        stress = arrays.get(stress_name)
        if points is None or ux is None:
            continue
        x = points[:, 0]
        coeff = np.polyfit(x, np.asarray(ux, dtype=float), 1)
        predicted = np.polyval(coeff, x)
        max_abs_ux_error = float(np.max(np.abs(ux - predicted)))
        ux_scale = max(float(np.max(np.abs(ux))), 1.0e-14)

        stress_mean = np.nan
        stress_std = np.nan
        stress_max_abs_error = np.nan
        if stress is not None and len(stress):
            stress = np.asarray(stress, dtype=float)
            stress_mean = float(np.mean(stress))
            stress_std = float(np.std(stress))
            stress_max_abs_error = float(np.max(np.abs(stress - expected_stress)))

        for trial_id, marker_suffix, trial_label in trial_cases:
            marker_path = vtu_dir / f"{stem}_{marker_suffix}.vtu"
            marker_arrays = read_vtu_arrays(marker_path) if marker_path.exists() else {}
            marker_points = marker_arrays.get("Points")
            rows.append(
                {
                    "element": element,
                    "trial": trial_id,
                    "trial_label": trial_label,
                    "mesh_nodes": len(points),
                    "marker_nodes": 0 if marker_points is None else len(marker_points),
                    "ux_slope": float(coeff[0]),
                    "ux_intercept": float(coeff[1]),
                    "max_abs_ux_fit_error": max_abs_ux_error,
                    "max_rel_ux_fit_error": max_abs_ux_error / ux_scale,
                    "stress_quantity": stress_name,
                    "expected_stress": expected_stress,
                    "stress_mean": stress_mean,
                    "stress_std": stress_std,
                    "stress_max_abs_error": stress_max_abs_error,
                    "passes_linear_ux": max_abs_ux_error <= 1.0e-12,
                    "passes_constant_stress": bool(
                        np.isnan(stress_std) or stress_std <= max(abs(expected_stress), 1.0) * 1.0e-12
                    ),
                    "source_vtu": str(base_path.relative_to(ROOT)),
                    "marker_vtu": (
                        str(marker_path.relative_to(ROOT)) if marker_path.exists() else ""
                    ),
                }
            )
    if not rows:
        return None
    summary = pd.DataFrame(rows)
    csv_path = OUT / "patch" / "patch_numeric_summary.csv"
    summary.to_csv(csv_path, index=False)

    element_summary = (
        summary.groupby("element", as_index=False)
        .agg(
            mesh_nodes=("mesh_nodes", "first"),
            max_marker_nodes=("marker_nodes", "max"),
            max_abs_ux_fit_error=("max_abs_ux_fit_error", "max"),
            max_rel_ux_fit_error=("max_rel_ux_fit_error", "max"),
            stress_quantity=("stress_quantity", "first"),
            expected_stress=("expected_stress", "first"),
            stress_mean=("stress_mean", "first"),
            stress_std=("stress_std", "max"),
            stress_max_abs_error=("stress_max_abs_error", "max"),
            all_linear_ux_pass=("passes_linear_ux", "all"),
            all_constant_stress_pass=("passes_constant_stress", "all"),
        )
        .sort_values("element")
    )
    element_summary.to_csv(OUT / "patch" / "patch_numeric_element_summary.csv", index=False)

    pivot = summary.pivot(
        index="element", columns="trial", values="max_abs_ux_fit_error"
    ).reset_index()
    pivot.to_csv(OUT / "patch" / "patch_numeric_error_pivot.csv", index=False)

    labels = [f"{row.element}\n{row.trial.split('_')[0]}" for row in summary.itertuples()]
    x = np.arange(len(summary))
    fig, axes = plt.subplots(3, 1, figsize=(10.8, 8.8), sharex=True)
    axes[0].bar(x, summary["max_abs_ux_fit_error"], color="#4477aa")
    axes[0].set_yscale("symlog", linthresh=1.0e-18)
    axes[0].set_ylabel("Max |UX - linear fit|")
    axes[0].set_title("Patch-test numerical summary by boundary form")
    axes[0].grid(True, axis="y", alpha=0.28)
    stress_pass = summary["passes_constant_stress"].astype(float)
    colors = ["#228833" if passed else "#cc6677" for passed in summary["passes_constant_stress"]]
    axes[1].bar(x, stress_pass, color=colors)
    axes[1].set_ylim(0.0, 1.15)
    axes[1].set_yticks([0.0, 1.0])
    axes[1].set_yticklabels(["fail", "pass"])
    axes[1].set_ylabel("Constant stress/force check")
    axes[1].grid(True, axis="y", alpha=0.28)
    max_stress_std = float(summary["stress_std"].fillna(0.0).max())
    axes[1].text(
        0.02,
        0.88,
        f"max std = {max_stress_std:.3g}",
        transform=axes[1].transAxes,
        fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
    )
    axes[2].bar(x, summary["marker_nodes"], color="#cc6677")
    axes[2].set_ylabel("Marked boundary nodes")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, rotation=35, ha="right")
    axes[2].grid(True, axis="y", alpha=0.28)
    path = OUT / "patch" / "patch_numeric_summary.png"
    savefig(path)
    return path


def build_coverage(rendered: dict[str, list[Path]]) -> Path:
    rows = []
    for element, meta in ELEMENT_CASES.items():
        cases = set(meta["cases"])
        rows.append(
            {
                "element": meta["label"],
                "convergence": "available" if rendered.get(f"{element}:convergence") else "missing",
                "patch_test": "available" if rendered.get(f"{element}:patch") else "missing",
                "abaqus_displacement_comparison": (
                    "available" if rendered.get(f"{element}:abaqus") else "missing"
                ),
                "validation_case": "available" if rendered.get(f"{element}:validation") else "missing",
                "source_cases": ", ".join(sorted(cases)),
            }
        )
    coverage = pd.DataFrame(rows)
    csv_path = OUT / "coverage" / "figure_coverage_matrix.csv"
    coverage.to_csv(csv_path, index=False)

    status = coverage[
        ["convergence", "patch_test", "abaqus_displacement_comparison", "validation_case"]
    ].replace({"available": 1, "missing": 0}).astype(float)
    fig, ax = plt.subplots(figsize=(8.8, 4.4))
    ax.imshow(status.to_numpy(), cmap="RdYlGn", vmin=0, vmax=1)
    ax.set_xticks(np.arange(status.shape[1]))
    ax.set_xticklabels(status.columns, rotation=25, ha="right")
    ax.set_yticks(np.arange(len(coverage)))
    ax.set_yticklabels(coverage["element"])
    for i in range(status.shape[0]):
        for j in range(status.shape[1]):
            ax.text(
                j,
                i,
                "OK" if status.iat[i, j] else "MISS",
                ha="center",
                va="center",
                fontsize=8,
                color="black",
            )
    ax.set_title("Figure coverage by element")
    path = OUT / "coverage" / "figure_coverage_matrix.png"
    savefig(path)
    return csv_path


def write_index(rendered: dict[str, list[Path]], coverage_csv: Path) -> None:
    lines = ["# Report figure index", ""]
    for key in sorted(rendered):
        lines.append(f"## {key}")
        for path in rendered[key]:
            lines.append(f"- `{path.relative_to(ROOT)}`")
        lines.append("")
    lines.append(f"Coverage matrix: `{coverage_csv.relative_to(ROOT)}`")
    lines.append("")
    lines.append("Missing entries are not inferred from single-mesh data.")
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ensure_dirs()
    rendered: dict[str, list[Path]] = {}

    detail_path = ROOT / "cases/abaqus_comparison_details.csv"
    summary_path = ROOT / "cases/abaqus_comparison.csv"
    if detail_path.exists() and summary_path.exists():
        details = pd.read_csv(detail_path)
        summary = pd.read_csv(summary_path)
        rendered.setdefault("abaqus:summary", []).append(plot_abaqus_summary(summary))
        for case in sorted(details["case"].unique()):
            element = next(
                (key for key, meta in ELEMENT_CASES.items() if case in meta["cases"]),
                None,
            )
            if element:
                path = plot_abaqus_case(case, details)
                rendered.setdefault(f"{element}:abaqus", []).append(path)

    for path in copy_patch_images():
        element = path.parent.name
        rendered.setdefault(f"{element}:patch", []).append(path)
    patch_summary = generate_patch_numeric_summary()
    if patch_summary:
        rendered.setdefault("patch:numeric_summary", []).append(patch_summary)

    for element, meta in ELEMENT_CASES.items():
        for case in meta["cases"]:
            path = plot_validation_case(element, case)
            if path:
                rendered.setdefault(f"{element}:validation", []).append(path)

    path = parse_plate4_convergence()
    if path:
        rendered.setdefault("plate4:convergence", []).append(path)
    generated_convergence = [
        ("bar", generate_bar_convergence()),
        ("beam", generate_beam_convergence()),
        ("hex8", generate_h8_convergence("hex8", 4, "Hex8 cantilever convergence")),
        (
            "bbarhex8",
            generate_h8_convergence("bbarhex8", 8, "BbarHex8 cantilever convergence"),
        ),
    ]
    for element, path in generated_convergence:
        if path:
            rendered.setdefault(f"{element}:convergence", []).append(path)
    rate_summary = write_convergence_rate_summary()
    if rate_summary:
        rendered.setdefault("convergence:rate_summary", []).append(rate_summary)

    coverage_csv = build_coverage(rendered)
    write_index(rendered, coverage_csv)


if __name__ == "__main__":
    main()
