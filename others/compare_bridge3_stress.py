import csv
import importlib.util
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ABAQUS_DIR = Path(r"C:\Users\Lenovo\Desktop\abaqus_bridge3_results")
STAP_DIR = Path(r"C:\Users\Lenovo\Desktop\Bridge3_results")


def load_converter():
    path = ROOT / "others/tools/abaqus_to_stappp.py"
    spec = importlib.util.spec_from_file_location("abaqus_to_stappp", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_stap_stress(out_path):
    group_stress = defaultdict(dict)
    current_group = None
    in_group = False
    for line in out_path.open(errors="ignore"):
        if "S T R E S S" in line and "G R O U P" in line:
            current_group = int(line.split()[-1])
            in_group = True
            continue
        if in_group:
            p = line.split()
            if p and p[0].isdigit():
                try:
                    group_stress[current_group][int(p[0])] = [float(x) for x in p[1:]]
                except ValueError:
                    pass
            elif "S T R E S S" in line and "G R O U P" in line:
                current_group = int(line.split()[-1])
    return group_stress


def stats(values):
    n = len(values)
    if n == 0:
        return None
    abs_vals = [abs(x) for x in values]
    return {
        "n": n,
        "max_abs": max(abs_vals),
        "mean_abs": sum(abs_vals) / n,
        "rmse": math.sqrt(sum(x * x for x in values) / n),
    }


conv = load_converter()
parts, instances, materials, assembly_nsets, assembly_surfaces, ties, boundaries, loads, gravity = conv.parse_inp(
    str(ABAQUS_DIR / "Bridge-3.inp")
)
nodes, node_map, flat_elements = conv.flatten_model(parts, instances)
nodes, node_map, flat_elements = conv.renumber_nodes(nodes, node_map, flat_elements, "rcm")

element_map = {}
for group, keys in ((1, ("C3D8R", "C3D8")),):
    idx = 1
    for etype in keys:
        for eid, conn, part_name, elset in flat_elements[etype]:
            # Need instance name, not part name. Reconstruct by matching element connectivity order from flatten_model
            idx += 1

# Rebuild element map directly from converter parse data and instance order.
group_index = {"solid": 1, "plate": 2, "beam": 3, "truss": 4}
group_counts = defaultdict(int)
for inst in instances:
    part = parts[inst["part"]]
    for etype, kind in (("C3D8R", "solid"), ("C3D8", "solid"), ("S4R", "plate"), ("B31", "beam"), ("T3D2", "truss")):
        for eid, _conn, _elset in part["elements"].get(etype, []):
            group_counts[kind] += 1
            element_map[(inst["name"].upper(), eid)] = (kind, group_index[kind], group_counts[kind])

stap_stress = parse_stap_stress(STAP_DIR / "Bridge-3.rcm.generated.out")

diffs = defaultdict(list)
matched = defaultdict(int)
missing = 0
worst = {}

with (ABAQUS_DIR / "Bridge-3_stress.csv").open(newline="", errors="ignore") as f:
    for row in csv.DictReader(f):
        key = (row["instance"].upper(), int(row["element"]))
        m = element_map.get(key)
        if not m:
            missing += 1
            continue
        kind, group, idx = m
        vals = stap_stress[group].get(idx)
        if vals is None:
            missing += 1
            continue
        matched[kind] += 1
        if kind == "solid":
            # STAP: SXX,SYY,SZZ,SXY,SYZ,SXZ,VonMises
            pairs = [
                ("solid_s11", vals[0], float(row["s11"])),
                ("solid_s22", vals[1], float(row["s22"])),
                ("solid_s33", vals[2], float(row["s33"])),
                ("solid_s12", vals[3], float(row["s12"])),
                ("solid_mises", vals[6], float(row["mises"])),
            ]
        elif kind == "truss":
            # STAP: Force, Stress. Abaqus truss stress is s11.
            pairs = [("truss_stress", vals[1], float(row["s11"]))]
        elif kind == "plate":
            # STAP plate output is engineering membrane/moment resultants, not direct Abaqus shell stresses.
            pairs = [("plate_sx_vs_s11", vals[0], float(row["s11"])), ("plate_sy_vs_s22", vals[1], float(row["s22"]))]
        else:
            continue
        for name, st, abq in pairs:
            d = st - abq
            diffs[name].append(d)
            if name not in worst or abs(d) > abs(worst[name]["diff"]):
                worst[name] = {
                    "diff": d,
                    "stap": st,
                    "abaqus": abq,
                    "instance": row["instance"],
                    "element": row["element"],
                    "stap_group_element": idx,
                }

print(f"matched={dict(matched)} missing={missing}")
for name in sorted(diffs):
    s = stats(diffs[name])
    w = worst[name]
    print(
        f"{name}: n={s['n']} max_abs_diff={s['max_abs']:.6g} mean_abs_diff={s['mean_abs']:.6g} rmse={s['rmse']:.6g}; "
        f"worst instance={w['instance']} elem={w['element']} stap_elem={w['stap_group_element']} "
        f"abaqus={w['abaqus']:.6g} stap={w['stap']:.6g} diff={w['diff']:.6g}"
    )
