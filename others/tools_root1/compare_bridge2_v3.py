"""Bridge-2 comparison using ODB global coordinates + displacements."""
import re, csv, math
from pathlib import Path

ROOT = Path(r"e:\学习\大三下\有限元基础\大作业\tests")

# ===== 1. Read STAP++ coords + disp =====
stap_out = ROOT / "data" / "data-3" / "Bridge-2.eigen_mpc.out"
print("Reading STAP++...")

stap_disp: dict[int, tuple[float, float, float]] = {}
stap_coords: dict[int, tuple[float, float, float]] = {}
section = None

with stap_out.open(errors="ignore") as f:
    for line in f:
        if "N O D A L   P O I N T   D A T A" in line:
            section = "nodes"
            continue
        if "D I S P L A C E M E N T S" in line:
            section = "disp"
            continue
        if section == "disp" and ("S T R E S S" in line or "R E A C T I O N" in line):
            break
        if section == "nodes" and ("L O A D" in line or "E L E M E N T" in line):
            section = None
            continue

        if section == "nodes":
            parts = line.split()
            if len(parts) >= 10:
                try:
                    node = int(parts[0])
                    x, y, z = float(parts[7]), float(parts[8]), float(parts[9])
                    stap_coords[node] = (x, y, z)
                except (ValueError, IndexError):
                    continue

        if section == "disp":
            m = re.match(
                r"^\s*(\d+)\s+([+-]?\d+\.\d+[Ee][+-]\d+)\s+([+-]?\d+\.\d+[Ee][+-]\d+)\s+([+-]?\d+\.\d+[Ee][+-]\d+)",
                line,
            )
            if m:
                stap_disp[int(m.group(1))] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))

print(f"  STAP++ coords: {len(stap_coords)}, disp: {len(stap_disp)}")

# ===== 2. Read Abaqus full CSV (with global coords!) =====
aba_full = ROOT / "data" / "data-3" / "Bridge-2_nodes_full.csv"
print(f"Reading Abaqus full CSV ({aba_full.stat().st_size / 1e6:.1f} MB)...")

# Build: rounded coordinates → list of (instance, node, u1, u2, u3)
# Because multiple instances might share same coordinates (Tie)
abaqus_by_coord: dict[tuple[float, float, float], list[tuple[str, int, float, float, float]]] = {}

with aba_full.open(newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            x = float(row["x"])
            y = float(row["y"])
            z = float(row["z"])
            u1 = float(row["u1"])
            u2 = float(row["u2"])
            u3 = float(row["u3"])
        except (ValueError, KeyError):
            continue
        rc = (round(x, 3), round(y, 3), round(z, 3))
        inst = row["instance"]
        node_id = int(row["node"])
        if rc not in abaqus_by_coord:
            abaqus_by_coord[rc] = []
        abaqus_by_coord[rc].append((inst, node_id, u1, u2, u3))

print(f"  Abaqus unique coords: {len(abaqus_by_coord)}")
print(f"  Abaqus total entries: {sum(len(v) for v in abaqus_by_coord.values())}")

# ===== 3. Match =====
print("Matching...")
matched = 0
errs_u1, errs_u2, errs_u3, errs_vec = [], [], [], []
max_err = {"stap": 0, "aba_inst": "", "aba_node": 0, "u1": 0., "u2": 0., "u3": 0., "vec": 0.}

for stap_node, coord in stap_coords.items():
    rc = (round(coord[0], 3), round(coord[1], 3), round(coord[2], 3))
    aba_entries = abaqus_by_coord.get(rc, [])
    if not aba_entries:
        continue
    # Use first match
    inst, aba_node, u1_aba, u2_aba, u3_aba = aba_entries[0]
    u_stap = stap_disp.get(stap_node)
    if u_stap is None:
        continue

    du1 = abs(u_stap[0] - u1_aba)
    du2 = abs(u_stap[1] - u2_aba)
    du3 = abs(u_stap[2] - u3_aba)
    vec = math.sqrt(du1**2 + du2**2 + du3**2)
    errs_u1.append(du1)
    errs_u2.append(du2)
    errs_u3.append(du3)
    errs_vec.append(vec)
    matched += 1
    if vec > max_err["vec"]:
        max_err = {"stap": stap_node, "aba_inst": inst, "aba_node": aba_node,
                   "u1": du1, "u2": du2, "u3": du3, "vec": vec}

def rms(lst):
    return math.sqrt(sum(e * e for e in lst) / len(lst))

print(f"\n{'='*60}")
print(f"  Bridge-2 STAP++ vs Abaqus")
print(f"{'='*60}")
print(f"  匹配节点数:        {matched}")
print(f"  位移向量 RMSE:     {rms(errs_vec):.6e}")
print(f"  U1 RMSE:           {rms(errs_u1):.6e}")
print(f"  U2 RMSE:           {rms(errs_u2):.6e}")
print(f"  U3 RMSE:           {rms(errs_u3):.6e}")
print(f"  Max |U1|:          {max(errs_u1):.6e}")
print(f"  Max |U2|:          {max(errs_u2):.6e}")
print(f"  Max |U3|:          {max(errs_u3):.6e}")
print(f"  Max |U|:           {max(errs_vec):.6e}")
print(f"     @ STAP++ node   {max_err['stap']}")
print(f"     @ Abaqus        {max_err['aba_inst']}:{max_err['aba_node']}")

# Quick sanity: show a few displacements
print(f"\n  Sample displacements (first 5 matched):")
count = 0
for stap_node, coord in sorted(stap_coords.items()):
    rc = (round(coord[0], 3), round(coord[1], 3), round(coord[2], 3))
    aba_entries = abaqus_by_coord.get(rc, [])
    if aba_entries and stap_node in stap_disp:
        u_s = stap_disp[stap_node]
        u_a = aba_entries[0]
        print(f"  node {stap_node}: STAP=({u_s[0]:.4e},{u_s[1]:.4e},{u_s[2]:.4e})  ABAQ=({u_a[2]:.4e},{u_a[3]:.4e},{u_a[4]:.4e})")
        count += 1
        if count >= 5:
            break

print(f"\nDone!")
