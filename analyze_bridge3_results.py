from pathlib import Path
import math
import os

root = Path(os.environ.get("BRIDGE3_ROOT", str(Path.home() / "bridge3_stappp_sparse")))
dat_path = Path(os.environ.get("BRIDGE3_DAT", str(root / "data/data-3/Bridge-3.rcm.generated.dat")))
out_path = Path(os.environ.get("BRIDGE3_OUT", str(root / "data/data-3/Bridge-3.rcm.generated.out")))


def parse_dat(path):
    nodes = {}
    groups = []
    with path.open(errors="ignore") as f:
        title = f.readline().strip()
        nnode, ngroup, ncase, mode = map(int, f.readline().split()[:4])
        for _ in range(nnode):
            p = f.readline().split()
            nid = int(p[0])
            nodes[nid] = tuple(float(x) for x in p[7:10])

        for _ in range(ncase):
            f.readline()
            header = f.readline().split()
            nloads = int(header[0])
            for _ in range(nloads):
                f.readline()
            if len(header) >= 6:
                for _ in range(int(header[5])):
                    f.readline()

        line = f.readline()
        while line and not line.strip():
            line = f.readline()
        buffered = line
        if len(line.split()) == 1:
            nmpc = int(float(line.split()[0]))
            for _ in range(nmpc):
                nterms = int(float(f.readline().split()[0]))
                for _ in range(nterms):
                    f.readline()
            buffered = None

        for gi in range(1, ngroup + 1):
            if buffered is None:
                line = f.readline()
                while line and not line.strip():
                    line = f.readline()
            else:
                line = buffered
                buffered = None
            etype, nelem, nmat = map(int, line.split()[:3])
            mats = [f.readline().split() for _ in range(nmat)]
            elems = []
            for _ in range(nelem):
                vals = f.readline().split()
                eid = int(vals[0])
                if etype in (4, 8):
                    conn = [int(x) for x in vals[1:9]]
                elif etype == 6:
                    conn = [int(x) for x in vals[1:5]]
                elif etype in (1, 5):
                    conn = [int(x) for x in vals[1:3]]
                else:
                    conn = []
                elems.append((eid, conn))
            groups.append({"group": gi, "etype": etype, "mats": mats, "elems": elems})
    return title, nodes, groups


def centroid(conn, nodes):
    if not conn:
        return (math.nan, math.nan, math.nan)
    xs = [nodes[n][0] for n in conn]
    ys = [nodes[n][1] for n in conn]
    zs = [nodes[n][2] for n in conn]
    n = len(conn)
    return (sum(xs) / n, sum(ys) / n, sum(zs) / n)


def update_abs(best, name, value, payload):
    if name not in best or abs(value) > abs(best[name]["value"]):
        best[name] = {"value": value, **payload}


def parse_out(path, nodes, groups):
    disp = {}
    disp_ext = {
        "ux_min": (math.inf, None),
        "ux_max": (-math.inf, None),
        "uy_min": (math.inf, None),
        "uy_max": (-math.inf, None),
        "uz_min": (math.inf, None),
        "uz_max": (-math.inf, None),
        "umag_max": (-math.inf, None),
        "rx_abs": (0.0, None),
        "ry_abs": (0.0, None),
        "rz_abs": (0.0, None),
    }
    stress_best = {}

    group_by_index = {g["group"]: g for g in groups}
    current_group = None
    current_type = None
    pending_group = None
    in_disp = False

    with path.open(errors="ignore") as f:
        for line in f:
            if "D I S P L A C E M E N T S" in line:
                in_disp = True
                continue
            if in_disp:
                p = line.split()
                if len(p) >= 7 and p[0].isdigit():
                    nid = int(p[0])
                    vals = [float(x) for x in p[1:7]]
                    disp[nid] = vals
                    ux, uy, uz, rx, ry, rz = vals
                    umag = math.sqrt(ux * ux + uy * uy + uz * uz)
                    for key, val in (("ux_min", ux), ("uy_min", uy), ("uz_min", uz)):
                        if val < disp_ext[key][0]:
                            disp_ext[key] = (val, nid)
                    for key, val in (("ux_max", ux), ("uy_max", uy), ("uz_max", uz)):
                        if val > disp_ext[key][0]:
                            disp_ext[key] = (val, nid)
                    if umag > disp_ext["umag_max"][0]:
                        disp_ext["umag_max"] = (umag, nid)
                    for key, val in (("rx_abs", rx), ("ry_abs", ry), ("rz_abs", rz)):
                        if abs(val) > abs(disp_ext[key][0]):
                            disp_ext[key] = (val, nid)
                    continue
                if "S T R E S S" in line:
                    in_disp = False

            if "S T R E S S" in line and "G R O U P" in line:
                parts = line.split()
                pending_group = int(parts[-1])
                current_group = None
                current_type = None
                continue

            if pending_group is not None and current_group is None:
                grp = group_by_index[pending_group]
                current_group = grp
                current_type = grp["etype"]

            p = line.split()
            if current_group and p and p[0].isdigit():
                eid = int(p[0])
                elems = current_group["elems"]
                if not (1 <= eid <= len(elems)):
                    continue
                conn = elems[eid - 1][1]
                c = centroid(conn, nodes)
                payload = {
                    "group": current_group["group"],
                    "etype": current_type,
                    "element": eid,
                    "centroid": c,
                    "conn": conn,
                }
                try:
                    vals = [float(x) for x in p[1:]]
                except ValueError:
                    continue

                if current_type in (4, 8) and len(vals) >= 7:
                    names = ["SXX", "SYY", "SZZ", "SXY", "SYZ", "SXZ", "VonMises"]
                    for name, value in zip(names, vals[:7]):
                        update_abs(stress_best, f"H8_{name}_abs", value, payload)
                    if "H8_VonMises_max" not in stress_best or vals[6] > stress_best["H8_VonMises_max"]["value"]:
                        stress_best["H8_VonMises_max"] = {"value": vals[6], **payload}
                elif current_type == 6 and len(vals) >= 6:
                    names = ["SX", "SY", "SXY", "MX", "MY", "MXY"]
                    for name, value in zip(names, vals[:6]):
                        update_abs(stress_best, f"Plate_{name}_abs", value, payload)
                elif current_type == 5 and len(vals) >= 6:
                    names = ["Axial", "ShearY", "ShearZ", "Torque", "MomentY", "MomentZ"]
                    for name, value in zip(names, vals[:6]):
                        update_abs(stress_best, f"Beam_{name}_abs", value, payload)
                elif current_type == 1 and len(vals) >= 2:
                    update_abs(stress_best, "Truss_Force_abs", vals[0], payload)
                    update_abs(stress_best, "Truss_Stress_abs", vals[1], payload)

    return disp, disp_ext, stress_best


def node_line(label, item, nodes):
    value, nid = item
    xyz = nodes.get(nid, (math.nan, math.nan, math.nan))
    print(f"{label}: value={value:.12g}, node={nid}, xyz=({xyz[0]:.12g}, {xyz[1]:.12g}, {xyz[2]:.12g})")


def best_line(name, item):
    c = item["centroid"]
    print(
        f"{name}: value={item['value']:.12g}, group={item['group']}, etype={item['etype']}, "
        f"element={item['element']}, centroid=({c[0]:.12g}, {c[1]:.12g}, {c[2]:.12g}), conn={item['conn'][:8]}"
    )


title, nodes, groups = parse_dat(dat_path)
disp, disp_ext, stress_best = parse_out(out_path, nodes, groups)

print("model")
print(f"title={title}")
print(f"nodes={len(nodes)}")
print("groups=" + ", ".join(f"g{g['group']}:etype{g['etype']}:{len(g['elems'])}" for g in groups))
print()

print("displacement_extrema")
for key in ("ux_min", "ux_max", "uy_min", "uy_max", "uz_min", "uz_max", "umag_max", "rx_abs", "ry_abs", "rz_abs"):
    node_line(key, disp_ext[key], nodes)
print()

print("element_extrema")
for key in (
    "H8_VonMises_max",
    "H8_SXX_abs",
    "H8_SYY_abs",
    "H8_SZZ_abs",
    "Plate_SX_abs",
    "Plate_SY_abs",
    "Plate_MX_abs",
    "Plate_MY_abs",
    "Beam_Axial_abs",
    "Beam_ShearY_abs",
    "Beam_ShearZ_abs",
    "Beam_Torque_abs",
    "Beam_MomentY_abs",
    "Beam_MomentZ_abs",
    "Truss_Force_abs",
    "Truss_Stress_abs",
):
    if key in stress_best:
        best_line(key, stress_best[key])
