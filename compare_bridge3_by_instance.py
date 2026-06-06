import csv
import importlib.util
import math
import os
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ABAQUS_DIR = Path(r"C:\Users\Lenovo\Desktop\abaqus_bridge3_results")
STAP_DIR = Path(os.environ.get("STAP_DIR", r"C:\Users\Lenovo\Desktop\Bridge3_results"))
STAP_OUT = Path(os.environ.get("STAP_OUT", str(STAP_DIR / "Bridge-3.rcm.generated.out")))


def load_converter():
    path = ROOT / "others/tools/abaqus_to_stappp.py"
    spec = importlib.util.spec_from_file_location("abaqus_to_stappp", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_stap_displacements(path):
    disp = {}
    in_disp = False
    for line in path.open(errors="ignore"):
        if "D I S P L A C E M E N T S" in line:
            in_disp = True
            continue
        if in_disp:
            p = line.split()
            if len(p) >= 7 and p[0].isdigit():
                disp[int(p[0])] = tuple(float(x) for x in p[1:4])
            elif "S T R E S S" in line:
                break
    return disp


def prefix(name):
    u = name.upper()
    for p in ("PART-FLOOR", "PART-PIER", "PART-RIVERBANK", "PART-SUPPORTBEAM", "PART-CABLE"):
        if u.startswith(p):
            return p
    return u


def norm(v):
    return math.sqrt(sum(x * x for x in v))


conv = load_converter()
parts, instances, materials, assembly_nsets, assembly_surfaces, ties, boundaries, loads, gravity = conv.parse_inp(
    str(ABAQUS_DIR / "Bridge-3.inp")
)
nodes, node_map, flat_elements = conv.flatten_model(parts, instances)
nodes, node_map, flat_elements = conv.renumber_nodes(nodes, node_map, flat_elements, "rcm")
node_map = {(a.upper(), b): c for (a, b), c in node_map.items()}
stap = parse_stap_displacements(STAP_OUT)

stats = defaultdict(lambda: {
    "n": 0,
    "sum_abs": [0.0, 0.0, 0.0],
    "sum_sq": [0.0, 0.0, 0.0],
    "max_abs": [0.0, 0.0, 0.0],
    "vec_sq": 0.0,
    "max_vec": (0.0, None),
    "abq_umag_max": (0.0, None),
    "stap_umag_max": (0.0, None),
})

with (ABAQUS_DIR / "Bridge-3_nodes.csv").open(newline="", errors="ignore") as f:
    for row in csv.DictReader(f):
        inst = row["instance"].upper()
        local = int(row["node"])
        sid = node_map[(inst, local)]
        abq = (float(row["u1"]), float(row["u2"]), float(row["u3"]))
        st = stap[sid]
        d = tuple(st[i] - abq[i] for i in range(3))
        for key in (inst, prefix(inst), "ALL"):
            s = stats[key]
            s["n"] += 1
            for i in range(3):
                s["sum_abs"][i] += abs(d[i])
                s["sum_sq"][i] += d[i] * d[i]
                s["max_abs"][i] = max(s["max_abs"][i], abs(d[i]))
            ve = norm(d)
            s["vec_sq"] += ve * ve
            if ve > s["max_vec"][0]:
                s["max_vec"] = (ve, (inst, local, sid, nodes[sid], abq, st, d))
            au = norm(abq)
            su = norm(st)
            if au > s["abq_umag_max"][0]:
                s["abq_umag_max"] = (au, (inst, local, sid, nodes[sid], abq, st, d))
            if su > s["stap_umag_max"][0]:
                s["stap_umag_max"] = (su, (inst, local, sid, nodes[sid], abq, st, d))


def row(key, s):
    n = s["n"]
    mean = [x / n for x in s["sum_abs"]]
    rmse = [math.sqrt(x / n) for x in s["sum_sq"]]
    vrmse = math.sqrt(s["vec_sq"] / n)
    print(
        f"{key}, n={n}, max_abs=({s['max_abs'][0]:.6g},{s['max_abs'][1]:.6g},{s['max_abs'][2]:.6g}), "
        f"mean_abs=({mean[0]:.6g},{mean[1]:.6g},{mean[2]:.6g}), "
        f"rmse=({rmse[0]:.6g},{rmse[1]:.6g},{rmse[2]:.6g}), vec_rmse={vrmse:.6g}, max_vec={s['max_vec'][0]:.6g}"
    )


print("PREFIX_STATS")
for key in ("ALL", "PART-FLOOR", "PART-PIER", "PART-RIVERBANK", "PART-SUPPORTBEAM", "PART-CABLE"):
    if key in stats:
        row(key, stats[key])
print()

print("WORST_INSTANCES_BY_VEC_RMSE")
items = []
for key, s in stats.items():
    if key.startswith("PART-") and key.count("-") >= 2 and s["n"] > 0:
        items.append((math.sqrt(s["vec_sq"] / s["n"]), key, s))
for _, key, s in sorted(items, reverse=True)[:12]:
    row(key, s)
