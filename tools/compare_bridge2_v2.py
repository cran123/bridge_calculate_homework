"""Bridge-2 STAP++ vs Abaqus comparison — coordinate-based matching."""
import re, csv, math
from pathlib import Path

ROOT = Path(r"e:\学习\大三下\有限元基础\大作业\tests")

# ===== 1. Read STAP++ displacements & coordinates from .out =====
stap_out = ROOT / "data" / "data-3" / "Bridge-2.eigen_mpc.out"
print(f"Reading STAP++ output ({stap_out.stat().st_size / 1e6:.1f} MB)...")

stap_disp: dict[int, tuple[float, float, float]] = {}
stap_coords: dict[int, tuple[float, float, float]] = {}
section = None  # "nodes" or "disp"

with stap_out.open(errors="ignore") as f:
    for line in f:
        # Detect sections
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
                stap_disp[int(m.group(1))] = (
                    float(m.group(2)),
                    float(m.group(3)),
                    float(m.group(4)),
                )

print(f"  STAP++ nodes with coords: {len(stap_coords)}")
print(f"  STAP++ nodes with disp:   {len(stap_disp)}")

# ===== 2. Read Abaqus INP — build coordinate → instance:node map =====
abaqus_inp = ROOT / "data" / "data-3" / "Bridge-2.inp"
print(f"Reading Abaqus INP ({abaqus_inp.stat().st_size / 1e6:.1f} MB)...")

abaqus_coords: dict[tuple[float, float, float], str] = {}  # rounded coord → "instance:node"
current_part = None
in_nodes = False

with abaqus_inp.open(errors="ignore") as f:
    for line in f:
        s = line.strip()
        u = s.upper()

        if u.startswith("*PART"):
            current_part = None
            for token in s.split(","):
                token = token.strip()
                if "=" in token:
                    k, v = token.split("=", 1)
                    if k.strip().upper() == "NAME":
                        current_part = v.strip()
                        break
            in_nodes = False
            continue

        if u.startswith("*END PART"):
            current_part = None
            in_nodes = False
            continue

        if u.startswith("*NODE"):
            in_nodes = bool(current_part)
            continue

        if in_nodes and s.startswith("*"):
            in_nodes = False
            continue

        if in_nodes and current_part:
            parts = s.split(",")
            if len(parts) >= 4:
                try:
                    lid = int(parts[0].strip())
                    x = float(parts[1].strip())
                    y = float(parts[2].strip())
                    z = float(parts[3].strip())
                    key = f"{current_part}:{lid}"
                    # Round to 3 decimals to permit tiny float differences
                    rc = (round(x, 3), round(y, 3), round(z, 3))
                    abaqus_coords[rc] = key
                except (ValueError, IndexError):
                    continue

print(f"  Abaqus INP unique coords: {len(abaqus_coords)}")

# ===== 3. Read Abaqus displacements =====
abaqus_csv = ROOT / "data" / "data-3" / "Bridge-2_nodes.csv"
print(f"Reading Abaqus CSV...")

abaqus_disp: dict[str, tuple[float, float, float]] = {}
with abaqus_csv.open(newline="") as f:
    for row in csv.DictReader(f):
        key = f"{row['instance']}:{row['node']}"
        abaqus_disp[key] = (float(row["u1"]), float(row["u2"]), float(row["u3"]))

print(f"  Abaqus displacement rows: {len(abaqus_disp)}")

# ===== 4. Match by coordinate =====
print("Matching by coordinate...")
matched = 0
errs_u1, errs_u2, errs_u3, errs_vec = [], [], [], []
max_err = {"stap_node": 0, "aba_key": "", "u1": 0.0, "u2": 0.0, "u3": 0.0, "vec": 0.0}

for stap_node, coord in stap_coords.items():
    rc = (round(coord[0], 3), round(coord[1], 3), round(coord[2], 3))
    aba_key = abaqus_coords.get(rc)
    if aba_key is None:
        continue
    u_aba = abaqus_disp.get(aba_key)
    u_stap = stap_disp.get(stap_node)
    if u_aba is None or u_stap is None:
        continue

    du1 = abs(u_stap[0] - u_aba[0])
    du2 = abs(u_stap[1] - u_aba[1])
    du3 = abs(u_stap[2] - u_aba[2])
    vec = math.sqrt(du1**2 + du2**2 + du3**2)
    errs_u1.append(du1)
    errs_u2.append(du2)
    errs_u3.append(du3)
    errs_vec.append(vec)
    matched += 1
    if vec > max_err["vec"]:
        max_err = {
            "stap_node": stap_node,
            "aba_key": aba_key,
            "u1": du1,
            "u2": du2,
            "u3": du3,
            "vec": vec,
        }


def rms(lst):
    return math.sqrt(sum(e * e for e in lst) / len(lst))


print(f"\n{'='*60}")
print(f"  Bridge-2 对比结果")
print(f"{'='*60}")
print(f"  匹配节点数: {matched}")
print(f"  位移向量 RMSE: {rms(errs_vec):.6e}")
print(f"  分量 RMSE: U1={rms(errs_u1):.6e}  U2={rms(errs_u2):.6e}  U3={rms(errs_u3):.6e}")
print(f"  最大绝对误差:")
print(f"    U1={max(errs_u1):.6e}, U2={max(errs_u2):.6e}, U3={max(errs_u3):.6e}")
print(f"    |U|={max(errs_vec):.6e} @ STAP++ node {max_err['stap_node']} ({max_err['aba_key']})")
