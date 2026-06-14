"""Quick Bridge-2 comparison: STAP++ vs Abaqus."""
import re, csv, math
from pathlib import Path

ROOT = Path(r"e:\学习\大三下\有限元基础\大作业\tests")

# --- Read STAP++ displacements ---
stap_out = ROOT / "data" / "data-3" / "Bridge-2.eigen_mpc.out"
text = stap_out.read_text(errors="ignore")
stap_disp: dict[int, tuple[float, float, float]] = {}
in_disp = False
for line in text.splitlines():
    if "D I S P L A C E M E N T S" in line:
        in_disp = True
        continue
    if in_disp and ("S T R E S S" in line or "R E A C T I O N" in line):
        break
    if not in_disp:
        continue
    m = re.match(
        r"^\s*(\d+)\s+([+-]?\d+\.\d+[Ee][+-]\d+)\s+([+-]?\d+\.\d+[Ee][+-]\d+)\s+([+-]?\d+\.\d+[Ee][+-]\d+)",
        line,
    )
    if m:
        stap_disp[int(m.group(1))] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
print(f"STAP++ nodes: {len(stap_disp)}")

# --- Read Abaqus CSV ---
abaqus_csv = ROOT / "data" / "data-3" / "Bridge-2_nodes.csv"
aba_disp: list[tuple[float, float, float]] = []
with abaqus_csv.open(newline="") as f:
    for row in csv.DictReader(f):
        aba_disp.append((float(row["u1"]), float(row["u2"]), float(row["u3"])))
print(f"Abaqus rows: {len(aba_disp)}")

# --- Direct comparison (same order) ---
n = min(len(stap_disp), len(aba_disp))
stap_list = [stap_disp[i + 1] for i in range(n)]

errs_u1, errs_u2, errs_u3, errs_vec = [], [], [], []
max_err = {"node": 0, "u1": 0.0, "u2": 0.0, "u3": 0.0, "vec": 0.0}
for i in range(n):
    du1 = abs(stap_list[i][0] - aba_disp[i][0])
    du2 = abs(stap_list[i][1] - aba_disp[i][1])
    du3 = abs(stap_list[i][2] - aba_disp[i][2])
    vec = math.sqrt(du1**2 + du2**2 + du3**2)
    errs_u1.append(du1)
    errs_u2.append(du2)
    errs_u3.append(du3)
    errs_vec.append(vec)
    if vec > max_err["vec"]:
        max_err = {"node": i + 1, "u1": du1, "u2": du2, "u3": du3, "vec": vec}

def rms(lst):
    return math.sqrt(sum(e * e for e in lst) / len(lst))

print(f"\n对比节点数: {n}")
print(f"Vec RMSE: {rms(errs_vec):.6e}")
print(f"U1 RMSE:  {rms(errs_u1):.6e}")
print(f"U2 RMSE:  {rms(errs_u2):.6e}")
print(f"U3 RMSE:  {rms(errs_u3):.6e}")
print(f"Max |U1|: {max(errs_u1):.6e}")
print(f"Max |U2|: {max(errs_u2):.6e}")
print(f"Max |U3|: {max(errs_u3):.6e}")
print(f"Max |U|:  {max(errs_vec):.6e} @ node {max_err['node']}")

# Check a few individual values for sanity
print(f"\nSanity check (first 5 nodes):")
for i in range(min(5, n)):
    du = math.sqrt(
        (stap_list[i][0] - aba_disp[i][0]) ** 2
        + (stap_list[i][1] - aba_disp[i][1]) ** 2
        + (stap_list[i][2] - aba_disp[i][2]) ** 2
    )
    print(f"  node {i+1}: STAP=({stap_list[i][0]:.6e},{stap_list[i][1]:.6e},{stap_list[i][2]:.6e})")
    print(f"           ABAQ=({aba_disp[i][0]:.6e},{aba_disp[i][1]:.6e},{aba_disp[i][2]:.6e})")
    print(f"           |du|={du:.6e}")
