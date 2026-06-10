#!/usr/bin/env python3

import argparse
import sys


def parse_options(line):
    parts = [p.strip() for p in line.split(",")]
    keyword = parts[0].strip().upper()
    opts = {}
    for part in parts[1:]:
        if "=" in part:
            k, v = part.split("=", 1)
            opts[k.strip().lower()] = v.strip()
        elif part:
            opts[part.strip().lower()] = True
    return keyword, opts


def parse_inp(path):
    nodes = {}
    elements = {}
    materials = {}
    solid_sections = {}
    truss_sections = {}
    boundaries = {}
    loads = []
    gravity = None

    current = None
    current_element = None
    current_elset = None
    current_material = None

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("**"):
                continue
            if line.startswith("*"):
                keyword, opts = parse_options(line[1:])
                current = keyword
                current_element = None
                current_elset = None

                if keyword == "ELEMENT":
                    current_element = opts.get("type", "").upper()
                    current_elset = opts.get("elset", None)
                elif keyword == "MATERIAL":
                    current_material = opts.get("name", None)
                    if current_material:
                        materials.setdefault(current_material, {"E": None, "nu": None, "rho": None})
                elif keyword == "SOLID SECTION":
                    elset = opts.get("elset")
                    mat = opts.get("material")
                    if elset and mat:
                        solid_sections[elset] = mat
                elif keyword == "TRUSS SECTION":
                    elset = opts.get("elset")
                    mat = opts.get("material")
                    area = opts.get("area")
                    if elset and mat and area:
                        try:
                            truss_sections[elset] = {"material": mat, "area": float(area)}
                        except ValueError:
                            pass
                continue

            if current == "NODE":
                parts = [p.strip() for p in line.split(",")]
                if len(parts) >= 4:
                    try:
                        nid = int(parts[0])
                        x = float(parts[1])
                        y = float(parts[2])
                        z = float(parts[3])
                        nodes[nid] = (x, y, z)
                    except ValueError:
                        continue

            elif current == "ELEMENT" and current_element:
                parts = [p.strip() for p in line.split(",")]
                if len(parts) >= 3:
                    try:
                        eid = int(parts[0])
                        conn = [int(p) for p in parts[1:] if p]
                        elements.setdefault(current_element, []).append((eid, conn, current_elset))
                    except ValueError:
                        continue

            elif current == "ELASTIC" and current_material:
                parts = [p.strip() for p in line.split(",")]
                if len(parts) >= 2:
                    try:
                        materials[current_material]["E"] = float(parts[0])
                        materials[current_material]["nu"] = float(parts[1])
                    except ValueError:
                        continue

            elif current == "DENSITY" and current_material:
                parts = [p.strip() for p in line.split(",")]
                if parts:
                    try:
                        materials[current_material]["rho"] = float(parts[0])
                    except ValueError:
                        continue

            elif current == "BOUNDARY":
                parts = [p.strip() for p in line.split(",")]
                if len(parts) >= 2:
                    try:
                        nid = int(parts[0])
                    except ValueError:
                        continue
                    dof1 = int(parts[1])
                    dof2 = dof1
                    if len(parts) >= 3 and parts[2]:
                        try:
                            dof2 = int(parts[2])
                        except ValueError:
                            dof2 = dof1
                    bc = boundaries.setdefault(nid, [0, 0, 0])
                    for d in range(dof1, dof2 + 1):
                        if 1 <= d <= 3:
                            bc[d - 1] = 1

            elif current == "CLOAD":
                parts = [p.strip() for p in line.split(",")]
                if len(parts) >= 3:
                    try:
                        nid = int(parts[0])
                        dof = int(parts[1])
                        val = float(parts[2])
                        if 1 <= dof <= 3:
                            loads.append((nid, dof, val))
                    except ValueError:
                        continue

            elif current == "DLOAD":
                parts = [p.strip() for p in line.split(",")]
                if len(parts) >= 2:
                    tag = parts[1].upper()
                    if tag == "GRAV" or parts[0].upper() == "GRAV":
                        nums = [p for p in parts if p and p.upper() != "GRAV"]
                        try:
                            if len(nums) >= 4:
                                g = float(nums[0])
                                gx = float(nums[1])
                                gy = float(nums[2])
                                gz = float(nums[3])
                                gravity = (g * gx, g * gy, g * gz)
                        except ValueError:
                            continue

    return nodes, elements, materials, solid_sections, truss_sections, boundaries, loads, gravity


def build_material_sets(elements, materials, solid_sections, truss_sections):
    h8_materials = []
    h8_map = {}
    truss_materials = []
    truss_map = {}

    for eid, conn, elset in elements.get("C3D8R", []) + elements.get("C3D8", []):
        mat_name = solid_sections.get(elset)
        if not mat_name or mat_name not in materials:
            continue
        mat = materials[mat_name]
        key = (mat["E"], mat["nu"], mat.get("rho") or 0.0)
        if key not in h8_map:
            h8_map[key] = len(h8_materials) + 1
            h8_materials.append(key)

    for eid, conn, elset in elements.get("T3D2", []):
        sec = truss_sections.get(elset)
        if not sec:
            continue
        mat_name = sec["material"]
        if mat_name not in materials:
            continue
        mat = materials[mat_name]
        key = (mat["E"], sec["area"], mat.get("rho") or 0.0)
        if key not in truss_map:
            truss_map[key] = len(truss_materials) + 1
            truss_materials.append(key)

    return h8_materials, h8_map, truss_materials, truss_map


def write_dat(path, nodes, elements, materials, solid_sections, truss_sections, boundaries, loads, gravity, cli_gravity):
    node_ids = sorted(nodes.keys())
    node_map = {nid: i + 1 for i, nid in enumerate(node_ids)}

    h8_elements = elements.get("C3D8R", []) + elements.get("C3D8", [])
    truss_elements = elements.get("T3D2", [])

    h8_materials, h8_map, truss_materials, truss_map = build_material_sets(elements, materials, solid_sections, truss_sections)

    groups = []
    if h8_elements:
        groups.append((4, h8_elements, h8_materials, h8_map))
    if truss_elements:
        groups.append((1, truss_elements, truss_materials, truss_map))

    nlcase = 1
    mode = 1

    if gravity is None and cli_gravity is not None:
        gravity = cli_gravity

    gravity_flag = 1 if gravity is not None else 0
    if gravity is None:
        gravity = (0.0, 0.0, 0.0)

    with open(path, "w", encoding="utf-8") as f:
        f.write("Converted from Abaqus inp\n")
        f.write(f"{len(node_ids)} {len(groups)} {nlcase} {mode}\n")

        for nid in node_ids:
            bc = boundaries.get(nid, [0, 0, 0])
            x, y, z = nodes[nid]
            f.write(f"{node_map[nid]} {bc[0]} {bc[1]} {bc[2]} {x} {y} {z}\n")

        f.write("1\n")
        f.write(f"{len(loads)} {gravity_flag} {gravity[0]} {gravity[1]} {gravity[2]}\n")
        for nid, dof, val in loads:
            if nid in node_map:
                f.write(f"{node_map[nid]} {dof} {val}\n")

        for element_type, ele_list, mats, mat_map in groups:
            f.write(f"{element_type} {len(ele_list)} {len(mats)}\n")

            if element_type == 4:
                for i, (E, nu, rho) in enumerate(mats, start=1):
                    f.write(f"{i} {E} {nu} {rho}\n")
                for idx, (eid, conn, elset) in enumerate(ele_list, start=1):
                    mat_name = solid_sections.get(elset)
                    mat = materials.get(mat_name, {})
                    key = (mat.get("E"), mat.get("nu"), mat.get("rho") or 0.0)
                    mset = mat_map.get(key, 1)
                    conn_new = [node_map[n] for n in conn]
                    f.write(f"{idx} " + " ".join(str(n) for n in conn_new) + f" {mset}\n")

            elif element_type == 1:
                for i, (E, area, rho) in enumerate(mats, start=1):
                    f.write(f"{i} {E} {area} {rho}\n")
                for idx, (eid, conn, elset) in enumerate(ele_list, start=1):
                    sec = truss_sections.get(elset, {})
                    mat_name = sec.get("material")
                    mat = materials.get(mat_name, {})
                    key = (mat.get("E"), sec.get("area"), mat.get("rho") or 0.0)
                    mset = mat_map.get(key, 1)
                    conn_new = [node_map[n] for n in conn]
                    f.write(f"{idx} {conn_new[0]} {conn_new[1]} {mset}\n")

    return {
        "nodes": len(node_ids),
        "h8": len(h8_elements),
        "truss": len(truss_elements),
        "materials": len(materials),
        "bcs": sum(sum(bc) for bc in boundaries.values()),
        "load_cases": nlcase,
    }


def main():
    parser = argparse.ArgumentParser(description="Convert Abaqus .inp to STAPpp .dat")
    parser.add_argument("inp", help="Input .inp file")
    parser.add_argument("-o", "--output", help="Output .dat file", default=None)
    parser.add_argument("--gravity", nargs=3, type=float, metavar=("GX", "GY", "GZ"),
                        help="Gravity vector if inp has no *Dload,GRAV")

    args = parser.parse_args()

    output = args.output
    if output is None:
        if args.inp.lower().endswith(".inp"):
            output = args.inp[:-4] + ".dat"
        else:
            output = args.inp + ".dat"

    nodes, elements, materials, solid_sections, truss_sections, boundaries, loads, gravity = parse_inp(args.inp)

    stats = write_dat(output, nodes, elements, materials, solid_sections, truss_sections,
                      boundaries, loads, gravity, tuple(args.gravity) if args.gravity else None)

    print(f"Nodes: {stats['nodes']}")
    print(f"C3D8R/H8: {stats['h8']}")
    print(f"T3D2/Truss3D2: {stats['truss']}")
    print(f"Materials: {stats['materials']}")
    print(f"Boundary DOFs: {stats['bcs']}")
    print(f"Load cases: {stats['load_cases']}")


if __name__ == "__main__":
    main()
