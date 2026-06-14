"""Check max-displacement-point error (现场位移最大值误差)."""
import re, csv, math
from pathlib import Path

ROOT = Path(r"e:\学习\大三下\有限元基础\大作业\tests")

CONFIGS = {
    "Bridge-1": {
        "stap": ROOT / "data" / "data-3" / "Bridge-1.rcm.generated.out",
        "aba": ROOT / "data" / "data-3" / "Bridge-1_nodes_full.csv",
    },
    "Bridge-2 (refrot+norot)": {
        "stap": ROOT / "data" / "data-3" / "Bridge-2.refrot_norot.out",
        "aba": ROOT / "data" / "data-3" / "Bridge-2_nodes_full.csv",
    },
    "Bridge-3": {
        "stap": ROOT / "bridge3" / "Bridge-3.rcm.refrot.norot.generated.out",
        "aba": ROOT / "cases" / "bridge3_run" / "abaqus" / "Bridge-3_nodes_full.csv",
    },
}


def read_stap_coords_disp(path: Path):
    coords, disp = {}, {}
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
                p = line.split()
                if len(p) >= 10:
                    try: coords[int(p[0])] = (float(p[7]), float(p[8]), float(p[9]))
                    except: continue
            if section == "disp":
                m = re.match(r"^\s*(\d+)\s+([+-]?\d+\.\d+[Ee][+-]\d+)\s+([+-]?\d+\.\d+[Ee][+-]\d+)\s+([+-]?\d+\.\d+[Ee][+-]\d+)", line)
                if m:
                    disp[int(m.group(1))] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
    return coords, disp


def read_aba_coords_disp(path: Path, has_coords=True):
    """Read Abaqus CSV. has_coords=True for _full.csv, False for _nodes.csv."""
    by_coord = {}
    if not path.exists():
        return by_coord
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                if has_coords:
                    x, y, z = float(row["x"]), float(row["y"]), float(row["z"])
                else:
                    x = y = z = 0.0  # no coords in plain _nodes.csv
                u1, u2, u3 = float(row["u1"]), float(row["u2"]), float(row["u3"])
            except: continue
            rc = (round(x, 3), round(y, 3), round(z, 3))
            if rc not in by_coord:
                by_coord[rc] = []
            by_coord[rc].append((u1, u2, u3))
    return by_coord


def read_aba_direct(path: Path):
    """Read _nodes.csv directly, return list of (u1,u2,u3) in order."""
    if not path.exists():
        return []
    result = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            try:
                result.append((float(row["u1"]), float(row["u2"]), float(row["u3"])))
            except: continue
    return result


for name, cfg in CONFIGS.items():
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")

    stap_c, stap_d = read_stap_coords_disp(cfg["stap"])
    if not stap_d:
        print(f"  STAP++ 输出为空或无法读取")
        continue

    # Bridge-3 uses _nodes.csv (no coords), others use _full.csv
    use_full = "_full" in str(cfg["aba"])

    if use_full:
        aba = read_aba_coords_disp(cfg["aba"], has_coords=True)
        print(f"  STAP: {len(stap_c)} coords, {len(stap_d)} disp | ABA: {len(aba)} unique coords")

        # ---------- 坐标匹配法 ----------
        # Find Abaqus node(s) with max displacement
        max_aba_mag = 0.0
        max_aba_coord = None
        max_aba_u = None
        for rc, entries in aba.items():
            for u in entries:
                mag = math.sqrt(u[0]**2 + u[1]**2 + u[2]**2)
                if mag > max_aba_mag:
                    max_aba_mag = mag
                    max_aba_coord = rc
                    max_aba_u = u

        # Find matching STAP++ node
        max_stap_u = None
        max_stap_node = None
        for sn, sc in stap_c.items():
            rc = (round(sc[0], 3), round(sc[1], 3), round(sc[2], 3))
            if rc == max_aba_coord:
                max_stap_u = stap_d.get(sn)
                max_stap_node = sn
                break

        if max_stap_u is None:
            print(f"  ✗ 无法匹配最大位移节点")
            continue

        abs_err = math.sqrt(
            (max_stap_u[0] - max_aba_u[0])**2 +
            (max_stap_u[1] - max_aba_u[1])**2 +
            (max_stap_u[2] - max_aba_u[2])**2
        )
        rel_err = abs_err / max_aba_mag * 100

        print(f"\n  --- 最大位移节点 ---")
        print(f"  STAP++ 节点: {max_stap_node}")
        print(f"  Abaqus |u|:  {max_aba_mag:.6e}")
        print(f"  STAP++ |u|:  {math.sqrt(sum(x**2 for x in max_stap_u)):.6e}")
        print(f"  |du|:        {abs_err:.6e}")
        print(f"  |du|/max|u|: {rel_err:.2f}%  {'✅ <5%' if rel_err < 5 else '❌ >5%'}")

        # Also compute normalized RMSE
        matched = 0
        errs = []
        for sn, sc in stap_c.items():
            rc = (round(sc[0], 3), round(sc[1], 3), round(sc[2], 3))
            entries = aba.get(rc, [])
            if not entries: continue
            u_s = stap_d.get(sn)
            if u_s is None: continue
            u_a = entries[0]
            du = math.sqrt((u_s[0]-u_a[0])**2 + (u_s[1]-u_a[1])**2 + (u_s[2]-u_a[2])**2)
            errs.append(du)
            matched += 1
        rmse = math.sqrt(sum(e**2 for e in errs)/len(errs))
        print(f"\n  归一化 RMSE (RMSE/max|u|): {rmse/max_aba_mag*100:.2f}%")
        print(f"  匹配节点: {matched}")

    else:
        # Bridge-3: use direct index matching (STAP nodes 1..N, ABA rows 1..N)
        aba_list = read_aba_direct(cfg["aba"])
        print(f"  STAP disp: {len(stap_d)} | ABA rows: {len(aba_list)}")

        n = min(len(stap_d), len(aba_list))
        # Find ABA max displacement
        max_aba_mag = 0.0
        max_aba_idx = 0
        max_aba_u = None
        for i in range(n):
            u = aba_list[i]
            mag = math.sqrt(u[0]**2 + u[1]**2 + u[2]**2)
            if mag > max_aba_mag:
                max_aba_mag = mag
                max_aba_idx = i
                max_aba_u = u

        # Corresponding STAP node
        stap_node = max_aba_idx + 1  # 1-based
        max_stap_u = stap_d.get(stap_node)
        if max_stap_u is None:
            print(f"  ✗ STAP++ node {stap_node} not found")
            continue

        abs_err = math.sqrt(
            (max_stap_u[0] - max_aba_u[0])**2 +
            (max_stap_u[1] - max_aba_u[1])**2 +
            (max_stap_u[2] - max_aba_u[2])**2
        )
        rel_err = abs_err / max_aba_mag * 100

        print(f"\n  --- 最大位移节点 (index-based) ---")
        print(f"  Abaqus row:  {max_aba_idx+1}")
        print(f"  Abaqus |u|:  {max_aba_mag:.6e}")
        print(f"  STAP++ |u|:  {math.sqrt(sum(x**2 for x in max_stap_u)):.6e}")
        print(f"  |du|:        {abs_err:.6e}")
        print(f"  |du|/max|u|: {rel_err:.2f}%  {'✅ <5%' if rel_err < 5 else '❌ >5%'}")

        # Compute normalized RMSE
        errs = []
        for i in range(n):
            u_s = stap_d.get(i+1)
            u_a = aba_list[i]
            if u_s is None: continue
            du = math.sqrt((u_s[0]-u_a[0])**2 + (u_s[1]-u_a[1])**2 + (u_s[2]-u_a[2])**2)
            errs.append(du)
        rmse = math.sqrt(sum(e**2 for e in errs)/len(errs))
        print(f"\n  归一化 RMSE (RMSE/max|u|): {rmse/max_aba_mag*100:.2f}%")
        print(f"  对比节点: {len(errs)}")

    # Also show top-5 displacement components comparison
    print(f"\n  --- Abaqus 位移最大的 5 个节点及其 STAP++ 误差 ---")
    all_errs = []
    if use_full:
        for sn, sc in stap_c.items():
            rc = (round(sc[0], 3), round(sc[1], 3), round(sc[2], 3))
            entries = aba.get(rc, [])
            if not entries: continue
            u_s = stap_d.get(sn)
            if u_s is None: continue
            u_a = entries[0]
            mag_a = math.sqrt(u_a[0]**2 + u_a[1]**2 + u_a[2]**2)
            du = math.sqrt((u_s[0]-u_a[0])**2 + (u_s[1]-u_a[1])**2 + (u_s[2]-u_a[2])**2)
            all_errs.append((mag_a, du, sn, u_s, u_a))
    else:
        for i in range(n):
            u_s = stap_d.get(i+1)
            u_a = aba_list[i]
            if u_s is None: continue
            mag_a = math.sqrt(u_a[0]**2 + u_a[1]**2 + u_a[2]**2)
            du = math.sqrt((u_s[0]-u_a[0])**2 + (u_s[1]-u_a[1])**2 + (u_s[2]-u_a[2])**2)
            all_errs.append((mag_a, du, i+1, u_s, u_a))

    all_errs.sort(key=lambda x: -x[0])
    for j, (mag_a, du, nd, u_s, u_a) in enumerate(all_errs[:5]):
        rel = du / mag_a * 100 if mag_a > 1e-12 else 0
        print(f"  #{j+1}: node {nd}, |u_aba|={mag_a:.4e}, |du|={du:.4e}, |du|/|u|={rel:.2f}%")

print("\nDone!")
