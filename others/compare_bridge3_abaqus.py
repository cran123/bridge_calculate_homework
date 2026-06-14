import csv
import importlib.util
import math
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ABAQUS_DIR = Path(r"C:\Users\Lenovo\Desktop\abaqus_bridge3_results")
STAP_DIR = Path(os.environ.get("STAP_DIR", r"C:\Users\Lenovo\Desktop\Bridge3_results"))

INP = ABAQUS_DIR / "Bridge-3.inp"
ABAQUS_NODES = ABAQUS_DIR / "Bridge-3_nodes.csv"
ABAQUS_STRESS = ABAQUS_DIR / "Bridge-3_stress.csv"
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
                disp[int(p[0])] = tuple(float(x) for x in p[1:7])
            elif "S T R E S S" in line:
                break
    return disp


def vec_norm(v):
    return math.sqrt(sum(x * x for x in v))


def main():
    conv = load_converter()
    parts, instances, materials, assembly_nsets, assembly_surfaces, ties, boundaries, loads, gravity = conv.parse_inp(str(INP))
    nodes, node_map, flat_elements = conv.flatten_model(parts, instances)
    nodes, node_map, flat_elements = conv.renumber_nodes(nodes, node_map, flat_elements, "rcm")
    normalized_node_map = {(inst.upper(), local): sid for (inst, local), sid in node_map.items()}
    stap_disp = parse_stap_displacements(STAP_OUT)

    rows = 0
    mapped = 0
    missing_map = 0
    missing_stap = 0
    comp_abs_max = [0.0, 0.0, 0.0]
    comp_abs_sum = [0.0, 0.0, 0.0]
    comp_sq_sum = [0.0, 0.0, 0.0]
    vector_sq_sum = 0.0
    max_vec_err = (-1.0, None)
    max_rel_vec_err = (-1.0, None)
    max_aba_umag = (-1.0, None)
    max_stap_umag_on_abq_rows = (-1.0, None)
    uz_min_abq = (math.inf, None)
    uz_min_stap_on_abq_rows = (math.inf, None)
    samples = []

    with ABAQUS_NODES.open(newline="", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows += 1
            inst = row["instance"]
            local = int(row["node"])
            abq = (float(row["u1"]), float(row["u2"]), float(row["u3"]))
            sid = normalized_node_map.get((inst.upper(), local))
            if sid is None:
                missing_map += 1
                continue
            st = stap_disp.get(sid)
            if st is None:
                missing_stap += 1
                continue
            st3 = st[:3]
            xyz = nodes[sid]
            mapped += 1
            diff = tuple(st3[i] - abq[i] for i in range(3))
            absdiff = tuple(abs(x) for x in diff)
            for i in range(3):
                comp_abs_max[i] = max(comp_abs_max[i], absdiff[i])
                comp_abs_sum[i] += absdiff[i]
                comp_sq_sum[i] += diff[i] * diff[i]
            vec_err = vec_norm(diff)
            vector_sq_sum += vec_err * vec_err
            abq_umag = vec_norm(abq)
            st_umag = vec_norm(st3)
            rel_vec = vec_err / max(abq_umag, 1e-12)

            payload = {
                "instance": inst,
                "local_node": local,
                "stap_node": sid,
                "xyz": xyz,
                "abaqus": abq,
                "stap": st3,
                "diff": diff,
                "vec_err": vec_err,
                "rel_vec": rel_vec,
                "abq_umag": abq_umag,
                "stap_umag": st_umag,
            }
            if vec_err > max_vec_err[0]:
                max_vec_err = (vec_err, payload)
            if rel_vec > max_rel_vec_err[0] and abq_umag > 1e-8:
                max_rel_vec_err = (rel_vec, payload)
            if abq_umag > max_aba_umag[0]:
                max_aba_umag = (abq_umag, payload)
            if st_umag > max_stap_umag_on_abq_rows[0]:
                max_stap_umag_on_abq_rows = (st_umag, payload)
            if abq[2] < uz_min_abq[0]:
                uz_min_abq = (abq[2], payload)
            if st3[2] < uz_min_stap_on_abq_rows[0]:
                uz_min_stap_on_abq_rows = (st3[2], payload)

    rmse = [math.sqrt(x / mapped) for x in comp_sq_sum]
    mean_abs = [x / mapped for x in comp_abs_sum]
    vec_rmse = math.sqrt(vector_sq_sum / mapped)

    def print_payload(name, item):
        value, p = item
        print(name)
        print(f"  metric={value:.12g}")
        print(f"  instance={p['instance']} local_node={p['local_node']} stap_node={p['stap_node']}")
        x, y, z = p["xyz"]
        print(f"  xyz=({x:.12g}, {y:.12g}, {z:.12g})")
        print(f"  abaqus=({p['abaqus'][0]:.12g}, {p['abaqus'][1]:.12g}, {p['abaqus'][2]:.12g}), |U|={p['abq_umag']:.12g}")
        print(f"  stap=({p['stap'][0]:.12g}, {p['stap'][1]:.12g}, {p['stap'][2]:.12g}), |U|={p['stap_umag']:.12g}")
        print(f"  diff=({p['diff'][0]:.12g}, {p['diff'][1]:.12g}, {p['diff'][2]:.12g}), vec_err={p['vec_err']:.12g}, rel={p['rel_vec']:.12g}")

    print("node_displacement_comparison")
    print(f"abaqus_rows={rows}")
    print(f"mapped_rows={mapped}")
    print(f"missing_map={missing_map}")
    print(f"missing_stap={missing_stap}")
    print(f"component_max_abs=U1:{comp_abs_max[0]:.12g}, U2:{comp_abs_max[1]:.12g}, U3:{comp_abs_max[2]:.12g}")
    print(f"component_mean_abs=U1:{mean_abs[0]:.12g}, U2:{mean_abs[1]:.12g}, U3:{mean_abs[2]:.12g}")
    print(f"component_rmse=U1:{rmse[0]:.12g}, U2:{rmse[1]:.12g}, U3:{rmse[2]:.12g}")
    print(f"vector_rmse={vec_rmse:.12g}")
    print()
    print_payload("max_vector_error", max_vec_err)
    print()
    print_payload("max_relative_vector_error_abq_umag_gt_1e-8", max_rel_vec_err)
    print()
    print_payload("abaqus_max_umag", max_aba_umag)
    print()
    print_payload("stap_max_umag_on_abaqus_rows", max_stap_umag_on_abq_rows)
    print()
    print_payload("abaqus_min_u3", uz_min_abq)
    print()
    print_payload("stap_min_u3_on_abaqus_rows", uz_min_stap_on_abq_rows)
    print()

    # Abaqus stress summary is extracted by instance/element, while STAP++ has renumbered element groups.
    # Summarize global Abaqus extrema for coarse magnitude checks.
    stress_rows = 0
    max_mises = (-1.0, None)
    max_abs = {k: (-1.0, None) for k in ("s11", "s22", "s33", "s12", "s13", "s23")}
    with ABAQUS_STRESS.open(newline="", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            stress_rows += 1
            mises = float(row["mises"])
            if mises > max_mises[0]:
                max_mises = (mises, row)
            for k in max_abs:
                v = abs(float(row[k]))
                if v > max_abs[k][0]:
                    max_abs[k] = (v, row)
    print("abaqus_stress_global_extrema")
    print(f"stress_rows={stress_rows}")
    print(f"max_mises={max_mises[0]:.12g}, instance={max_mises[1]['instance']}, element={max_mises[1]['element']}")
    for k, (v, row) in max_abs.items():
        print(f"max_abs_{k}={v:.12g}, instance={row['instance']}, element={row['element']}")


if __name__ == "__main__":
    main()
