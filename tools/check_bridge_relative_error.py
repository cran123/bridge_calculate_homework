"""Check if bridge node errors are within 5% of Abaqus displacement."""
import re, csv, math, sys
from pathlib import Path

ROOT = Path(r"e:\学习\大三下\有限元基础\大作业\tests")

# ---- Config ----
BRIDGES = {
    "Bridge-1": {
        "stap_out": ROOT / "data" / "data-3" / "Bridge-1.rcm.generated.out",
        "aba_full": ROOT / "data" / "data-3" / "Bridge-1_nodes_full.csv",
    },
    "Bridge-2": {
        "stap_out": ROOT / "data" / "data-3" / "Bridge-2.eigen_mpc.out",
        "aba_full": ROOT / "data" / "data-3" / "Bridge-2_nodes_full.csv",
    },
    "Bridge-3": {
        "stap_out": ROOT / "cases" / "bridge3_run" / "stappp" / "bridge3_run.out",
        "aba_full": ROOT / "cases" / "bridge3_run" / "abaqus" / "Bridge-3_nodes_full.csv",
    },
}


def read_stappp_coords_disp(path: Path):
    """Read STAP++ node coordinates and displacements from .out."""
    coords = {}
    disp = {}
    section = None
    if not path.exists():
        return coords, disp
    with path.open(errors="ignore") as f:
        for line in f:
            if "N O D A L   P O I N T   D A T A" in line:
                section = "nodes"; continue
            if "D I S P L A C E M E N T S" in line:
                section = "disp"; continue
            if section == "disp" and ("S T R E S S" in line or "R E A C T I O N" in line):
                break
            if section == "nodes" and ("L O A D" in line or "E L E M E N T" in line):
                section = None; continue
            if section == "nodes":
                parts = line.split()
                if len(parts) >= 10:
                    try:
                        node = int(parts[0])
                        coords[node] = (float(parts[7]), float(parts[8]), float(parts[9]))
                    except: continue
            if section == "disp":
                m = re.match(r"^\s*(\d+)\s+([+-]?\d+\.\d+[Ee][+-]\d+)\s+([+-]?\d+\.\d+[Ee][+-]\d+)\s+([+-]?\d+\.\d+[Ee][+-]\d+)", line)
                if m:
                    disp[int(m.group(1))] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
    return coords, disp


def read_abaqus_full(path: Path):
    """Read Abaqus nodes_full.csv → {(x,y,z)_rounded: list of (u1,u2,u3)}."""
    if not path.exists():
        return {}
    by_coord = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            try:
                x, y, z = float(row["x"]), float(row["y"]), float(row["z"])
                u1, u2, u3 = float(row["u1"]), float(row["u2"]), float(row["u3"])
            except: continue
            rc = (round(x, 3), round(y, 3), round(z, 3))
            if rc not in by_coord:
                by_coord[rc] = []
            by_coord[rc].append((u1, u2, u3))
    return by_coord


def analyze(bridge_name, cfg):
    print(f"\n{'='*70}")
    print(f"  {bridge_name}")
    print(f"{'='*70}")

    stap_c, stap_d = read_stappp_coords_disp(cfg["stap_out"])
    aba = read_abaqus_full(cfg["aba_full"])

    print(f"  STAP++ nodes: {len(stap_c)} coords, {len(stap_d)} disp")
    print(f"  Abaqus unique coords: {len(aba)}")

    # Match and compute errors
    matched = 0
    abs_errs = []    # absolute vector errors
    rel_errs = []    # relative errors = |du|/|u_aba|, thresholded
    abaqus_mags = [] # |u_abaqus|
    stap_mags = []   # |u_stappp|
    max_rel = {"node": 0, "rel": 0.0, "stap_mag": 0.0, "aba_mag": 0.0}

    TINY = 1e-12  # treat displacements smaller than this as zero

    for stap_node, coord in stap_c.items():
        rc = (round(coord[0], 3), round(coord[1], 3), round(coord[2], 3))
        entries = aba.get(rc, [])
        if not entries:
            continue
        u_aba = entries[0]  # (u1,u2,u3)
        u_stap = stap_d.get(stap_node)
        if u_stap is None:
            continue

        du1 = u_stap[0] - u_aba[0]
        du2 = u_stap[1] - u_aba[1]
        du3 = u_stap[2] - u_aba[2]
        abs_err = math.sqrt(du1**2 + du2**2 + du3**2)

        mag_aba = math.sqrt(u_aba[0]**2 + u_aba[1]**2 + u_aba[2]**2)
        mag_stap = math.sqrt(u_stap[0]**2 + u_stap[1]**2 + u_stap[2]**2)

        # Relative error: for nodes with very small displacement, use absolute threshold
        if mag_aba > TINY:
            rel_err = abs_err / mag_aba
        else:
            # Node has ~zero displacement → check if STAP++ also gives ~zero
            rel_err = abs_err if abs_err > TINY else 0.0

        matched += 1
        abs_errs.append(abs_err)
        rel_errs.append(rel_err)
        abaqus_mags.append(mag_aba)
        stap_mags.append(mag_stap)

        if rel_err > max_rel["rel"]:
            max_rel = {"node": stap_node, "rel": rel_err, "stap_mag": mag_stap, "aba_mag": mag_aba,
                       "abs_err": abs_err, "u_stap": u_stap, "u_aba": u_aba}

    # Statistics
    n = len(rel_errs)
    if n == 0:
        print("  No matched nodes!")
        return

    within_1pct = sum(1 for r in rel_errs if r < 0.01)
    within_5pct = sum(1 for r in rel_errs if r < 0.05)
    within_10pct = sum(1 for r in rel_errs if r < 0.10)
    within_50pct = sum(1 for r in rel_errs if r < 0.50)

    max_aba_mag = max(abaqus_mags)
    mean_aba_mag = sum(abaqus_mags) / n
    median_aba_mag = sorted(abaqus_mags)[n // 2]

    print(f"\n  匹配节点数: {matched}")
    print(f"  Abaqus 位移量级: max={max_aba_mag:.4f}, mean={mean_aba_mag:.4e}, median={median_aba_mag:.4e}")
    print(f"\n  --- 相对误差分布 ---")
    print(f"  < 1%:    {within_1pct:>6d} / {n}  ({100*within_1pct/n:.1f}%)")
    print(f"  < 5%:    {within_5pct:>6d} / {n}  ({100*within_5pct/n:.1f}%)")
    print(f"  < 10%:   {within_10pct:>6d} / {n}  ({100*within_10pct/n:.1f}%)")
    print(f"  < 50%:   {within_50pct:>6d} / {n}  ({100*within_50pct/n:.1f}%)")
    print(f"  Max rel: {max_rel['rel']*100:.2f}% @ node {max_rel['node']} (|u_aba|={max_rel['aba_mag']:.4e}, |u_stap|={max_rel['stap_mag']:.4e})")
    
    # Mean/median relative error (only for nodes with meaningful displacement)
    meaningful = [(r, m) for r, m in zip(rel_errs, abaqus_mags) if m > 1e-6]
    if meaningful:
        m_rel = sorted(r for r, _ in meaningful)
        n_m = len(m_rel)
        print(f"\n  有意义位移节点(>1e-6): {n_m}/{n}")
        print(f"  平均相对误差: {sum(r for r,_ in meaningful)/n_m*100:.2f}%")
        print(f"  中位相对误差: {m_rel[n_m//2]*100:.2f}%")

    # Show worst nodes
    sorted_by_rel = sorted(zip(rel_errs, abs_errs, abaqus_mags), key=lambda x: -x[0])
    print(f"\n  --- 相对误差最大的 5 个节点 ---")
    for i, (r, a, m) in enumerate(sorted_by_rel[:5]):
        print(f"  #{i+1}: rel={r*100:.1f}%, |du|={a:.4e}, |u_aba|={m:.4e}")


def main():
    for name in ["Bridge-1", "Bridge-2"]:
        analyze(name, BRIDGES[name])
    # Bridge-3 needs full CSV — check if it exists, else skip
    b3_full = BRIDGES["Bridge-3"]["aba_full"]
    if b3_full.exists():
        analyze("Bridge-3", BRIDGES["Bridge-3"])
    else:
        print("\n[Bridge-3] _nodes_full.csv 不存在，跳过。运行 abaqus python extract_odb_coords.py 生成。")

if __name__ == "__main__":
    main()
