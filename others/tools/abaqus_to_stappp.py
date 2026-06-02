#!/usr/bin/env python3

import argparse
import math


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


def numbers(line, cast=float):
    vals = []
    for part in line.split(","):
        part = part.strip()
        if not part:
            continue
        vals.append(cast(part))
    return vals


def parse_ints(line):
    return numbers(line, int)


def add_generate_ids(target, values):
    if len(values) >= 3:
        start, end, step = values[:3]
        target.extend(range(start, end + 1, step))


def matmul_vec(m, v):
    return (
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    )


def rotation_matrix(axis, angle_deg):
    ax, ay, az = axis
    n = math.sqrt(ax * ax + ay * ay + az * az)
    if n == 0.0:
        return ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))

    x, y, z = ax / n, ay / n, az / n
    a = math.radians(angle_deg)
    c, s, t = math.cos(a), math.sin(a), 1.0 - math.cos(a)
    return (
        (t * x * x + c, t * x * y - s * z, t * x * z + s * y),
        (t * x * y + s * z, t * y * y + c, t * y * z - s * x),
        (t * x * z - s * y, t * y * z + s * x, t * z * z + c),
    )


def rotate_about(point, axis_p1, matrix):
    rel = (point[0] - axis_p1[0], point[1] - axis_p1[1], point[2] - axis_p1[2])
    rot = matmul_vec(matrix, rel)
    return (rot[0] + axis_p1[0], rot[1] + axis_p1[1], rot[2] + axis_p1[2])


def box_section_properties(values):
    if len(values) < 6:
        return (1.0, 1.0, 1.0, 1.0)

    b, h, t_top, t_bottom, t_left, t_right = values[:6]
    bi = max(b - t_left - t_right, 1.0e-12)
    hi = max(h - t_top - t_bottom, 1.0e-12)
    area = b * h - bi * hi
    iy = (b * h ** 3 - bi * hi ** 3) / 12.0
    iz = (h * b ** 3 - hi * bi ** 3) / 12.0
    return (area, iy, iz, iy + iz)


def section_material(part, elset):
    return (
        part["solid_sections"].get(elset)
        or part["shell_sections"].get(elset)
        or part["beam_sections"].get(elset)
    )


def resolve_element_elset(part, etype, eid, explicit_elset):
    if explicit_elset:
        return explicit_elset

    if etype in ("C3D8", "C3D8R", "T3D2"):
        candidates = part["solid_sections"]
    elif etype == "S4R":
        candidates = part["shell_sections"]
    elif etype == "B31":
        candidates = part["beam_sections"]
    else:
        candidates = {}

    for elset in candidates:
        if eid in part["elsets"].get(elset, []):
            return elset
    return explicit_elset


def parse_inp(path):
    parts = {}
    materials = {}
    instances = []
    assembly_nsets = {}
    assembly_surfaces = {}
    ties = []
    boundaries = []
    loads = []
    gravity = None

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    i = 0
    current_part = None
    current_material = None
    in_assembly = False

    while i < len(lines):
        raw = lines[i]
        line = raw.strip()
        i += 1
        if not line or line.startswith("**"):
            continue

        if not line.startswith("*"):
            continue

        keyword, opts = parse_options(line[1:])

        if keyword == "PART":
            name = opts.get("name")
            current_part = {
                "name": name,
                "nodes": {},
                "elements": {},
                "nsets": {},
                "elsets": {},
                "solid_sections": {},
                "shell_sections": {},
                "beam_sections": {},
            }
            parts[name] = current_part
            continue

        if keyword == "END PART":
            current_part = None
            continue

        if keyword == "ASSEMBLY":
            in_assembly = True
            continue

        if keyword == "END ASSEMBLY":
            in_assembly = False
            continue

        if keyword == "MATERIAL":
            current_material = opts.get("name")
            if current_material:
                materials.setdefault(current_material, {"E": None, "nu": None, "rho": 0.0})
            continue

        if keyword == "DENSITY" and current_material:
            while i < len(lines) and not lines[i].strip().startswith("*"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    vals = numbers(data)
                    if vals:
                        materials[current_material]["rho"] = vals[0]
                    break
            continue

        if keyword == "ELASTIC" and current_material:
            while i < len(lines) and not lines[i].strip().startswith("*"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    vals = numbers(data)
                    if len(vals) >= 2:
                        materials[current_material]["E"] = vals[0]
                        materials[current_material]["nu"] = vals[1]
                    break
            continue

        if current_part and keyword == "NODE":
            while i < len(lines) and not lines[i].strip().startswith("*"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    vals = numbers(data)
                    if len(vals) >= 4:
                        current_part["nodes"][int(vals[0])] = (vals[1], vals[2], vals[3])
            continue

        if current_part and keyword == "ELEMENT":
            etype = opts.get("type", "").upper()
            elset = opts.get("elset")
            elems = current_part["elements"].setdefault(etype, [])
            while i < len(lines) and not lines[i].strip().startswith("*"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    vals = parse_ints(data)
                    if len(vals) >= 3:
                        elems.append((vals[0], vals[1:], elset))
                        if elset:
                            current_part["elsets"].setdefault(elset, []).append(vals[0])
            continue

        if current_part and keyword in ("NSET", "ELSET"):
            name = opts.get("nset") or opts.get("elset")
            target = current_part["nsets" if keyword == "NSET" else "elsets"].setdefault(name, [])
            while i < len(lines) and not lines[i].strip().startswith("*"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    vals = parse_ints(data)
                    if "generate" in opts:
                        add_generate_ids(target, vals)
                    else:
                        target.extend(vals)
            continue

        if current_part and keyword == "SOLID SECTION":
            elset = opts.get("elset")
            mat = opts.get("material")
            area = None
            while i < len(lines) and not lines[i].strip().startswith("*"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    vals = numbers(data)
                    area = vals[0] if vals else None
                    break
            if elset and mat:
                current_part["solid_sections"][elset] = {"material": mat, "area": area}
            continue

        if current_part and keyword == "SHELL SECTION":
            elset = opts.get("elset")
            mat = opts.get("material")
            thickness = None
            while i < len(lines) and not lines[i].strip().startswith("*"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    vals = numbers(data)
                    thickness = vals[0] if vals else None
                    break
            if elset and mat:
                current_part["shell_sections"][elset] = {"material": mat, "thickness": thickness or 1.0}
            continue

        if current_part and keyword == "BEAM SECTION":
            elset = opts.get("elset")
            mat = opts.get("material")
            profile_values = []
            ref = (0.0, 0.0, 1.0)
            data_lines = []
            while i < len(lines) and not lines[i].strip().startswith("*"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    data_lines.append(data)
            if data_lines:
                profile_values = numbers(data_lines[0])
            if len(data_lines) > 1:
                vals = numbers(data_lines[1])
                if len(vals) >= 3:
                    ref = (vals[0], vals[1], vals[2])
            if elset and mat:
                area, iy, iz, j = box_section_properties(profile_values)
                current_part["beam_sections"][elset] = {
                    "material": mat,
                    "area": area,
                    "iy": iy,
                    "iz": iz,
                    "j": j,
                    "ref": ref,
                }
            continue

        if in_assembly and keyword == "INSTANCE":
            name = opts.get("name")
            part_name = opts.get("part")
            transform_lines = []
            while i < len(lines) and not lines[i].strip().upper().startswith("*END INSTANCE"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    transform_lines.append(numbers(data))
            if i < len(lines) and lines[i].strip().upper().startswith("*END INSTANCE"):
                i += 1

            translation = (0.0, 0.0, 0.0)
            rot_p1 = rot_matrix = None
            if transform_lines:
                vals = transform_lines[0]
                if len(vals) >= 3:
                    translation = (vals[0], vals[1], vals[2])
            if len(transform_lines) > 1 and len(transform_lines[1]) >= 7:
                vals = transform_lines[1]
                rot_p1 = (vals[0], vals[1], vals[2])
                rot_p2 = (vals[3], vals[4], vals[5])
                axis = (rot_p2[0] - rot_p1[0], rot_p2[1] - rot_p1[1], rot_p2[2] - rot_p1[2])
                rot_matrix = rotation_matrix(axis, vals[6])

            instances.append({
                "name": name,
                "part": part_name,
                "translation": translation,
                "rot_p1": rot_p1,
                "rot_matrix": rot_matrix,
            })
            continue

        if in_assembly and keyword == "NSET":
            name = opts.get("nset")
            inst = opts.get("instance")
            key = name
            target = assembly_nsets.setdefault(key, [])
            while i < len(lines) and not lines[i].strip().startswith("*"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    vals = parse_ints(data)
                    ids = []
                    if "generate" in opts:
                        add_generate_ids(ids, vals)
                    else:
                        ids.extend(vals)
                    target.extend((inst, nid) for nid in ids)
            continue

        if in_assembly and keyword == "SURFACE":
            name = opts.get("name")
            surface_type = opts.get("type", "").upper()
            refs = []
            while i < len(lines) and not lines[i].strip().startswith("*"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    parts_line = [p.strip() for p in data.split(",")]
                    if parts_line and parts_line[0]:
                        refs.append(parts_line[0])
            if name and surface_type == "NODE":
                assembly_surfaces[name] = refs
            continue

        if in_assembly and keyword == "TIE":
            while i < len(lines) and not lines[i].strip().startswith("*"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    parts_line = [p.strip() for p in data.split(",")]
                    if len(parts_line) >= 2:
                        ties.append((parts_line[0], parts_line[1]))
                    break
            continue

        if keyword == "BOUNDARY":
            while i < len(lines) and not lines[i].strip().startswith("*"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    parts_line = [p.strip() for p in data.split(",")]
                    if len(parts_line) >= 2:
                        ref = parts_line[0]
                        dof1 = int(parts_line[1])
                        dof2 = int(parts_line[2]) if len(parts_line) >= 3 and parts_line[2] else dof1
                        boundaries.append((ref, dof1, dof2))
            continue

        if keyword == "CLOAD":
            while i < len(lines) and not lines[i].strip().startswith("*"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    parts_line = [p.strip() for p in data.split(",")]
                    if len(parts_line) >= 3:
                        loads.append((parts_line[0], int(parts_line[1]), float(parts_line[2])))
            continue

        if keyword == "DLOAD":
            while i < len(lines) and not lines[i].strip().startswith("*"):
                data = lines[i].strip()
                i += 1
                if data and not data.startswith("**"):
                    parts_line = [p.strip() for p in data.split(",") if p.strip()]
                    if "GRAV" in [p.upper() for p in parts_line]:
                        vals = [float(p) for p in parts_line if p.upper() != "GRAV"]
                        if len(vals) >= 4:
                            gravity = (vals[0] * vals[1], vals[0] * vals[2], vals[0] * vals[3])
            continue

    return parts, instances, materials, assembly_nsets, assembly_surfaces, ties, boundaries, loads, gravity


def flatten_model(parts, instances):
    nodes = {}
    node_map = {}
    flat_elements = {"C3D8R": [], "C3D8": [], "S4R": [], "B31": [], "T3D2": []}
    next_node = 1

    for inst in instances:
        part = parts[inst["part"]]
        trans = inst["translation"]
        for local_id in sorted(part["nodes"]):
            x, y, z = part["nodes"][local_id]
            p = (x + trans[0], y + trans[1], z + trans[2])
            if inst["rot_matrix"] is not None:
                p = rotate_about(p, inst["rot_p1"], inst["rot_matrix"])
            node_map[(inst["name"], local_id)] = next_node
            nodes[next_node] = p
            next_node += 1

        for etype, elems in part["elements"].items():
            if etype not in flat_elements:
                continue
            for eid, conn, elset in elems:
                flat_conn = [node_map[(inst["name"], nid)] for nid in conn]
                resolved_elset = resolve_element_elset(part, etype, eid, elset)
                flat_elements[etype].append((eid, flat_conn, inst["part"], resolved_elset))

    return nodes, node_map, flat_elements


def nearest_node(source, candidates, nodes):
    sx, sy, sz = nodes[source]
    best = None
    best_d2 = None
    for candidate in candidates:
        cx, cy, cz = nodes[candidate]
        d2 = (sx - cx) ** 2 + (sy - cy) ** 2 + (sz - cz) ** 2
        if best_d2 is None or d2 < best_d2:
            best = candidate
            best_d2 = d2
    return best


def merge_tied_nodes(nodes, node_map, flat_elements, assembly_nsets, assembly_surfaces, ties):
    parent = {nid: nid for nid in nodes}

    def find(nid):
        while parent[nid] != nid:
            parent[nid] = parent[parent[nid]]
            nid = parent[nid]
        return nid

    def union(slave, master):
        parent[find(slave)] = find(master)

    def surface_nodes(surface_name):
        result = []
        for set_name in assembly_surfaces.get(surface_name, []):
            for inst, local_id in assembly_nsets.get(set_name, []):
                gid = node_map.get((inst, local_id))
                if gid is not None:
                    result.append(gid)
        return result

    merged_pairs = 0
    for slave_surface, master_surface in ties:
        slave_nodes = surface_nodes(slave_surface)
        master_nodes = surface_nodes(master_surface)
        if not slave_nodes or not master_nodes:
            continue
        for slave in slave_nodes:
            master = nearest_node(slave, master_nodes, nodes)
            if master is not None:
                union(slave, master)
                merged_pairs += 1

    rep_to_new = {}
    old_to_new = {}
    new_nodes = {}
    for old_id in sorted(nodes):
        rep = find(old_id)
        if rep not in rep_to_new:
            rep_to_new[rep] = len(rep_to_new) + 1
            new_nodes[rep_to_new[rep]] = nodes[rep]
        old_to_new[old_id] = rep_to_new[rep]

    for key, old_id in list(node_map.items()):
        node_map[key] = old_to_new[old_id]

    for etype, elems in flat_elements.items():
        for idx, (eid, conn, part_name, elset) in enumerate(elems):
            flat_elements[etype][idx] = (eid, [old_to_new[nid] for nid in conn], part_name, elset)

    return new_nodes, node_map, flat_elements, merged_pairs


def apply_boundaries(boundary_specs, assembly_nsets, node_map):
    bcs = {}
    for ref, dof1, dof2 in boundary_specs:
        refs = assembly_nsets.get(ref, [])
        for inst, local_id in refs:
            gid = node_map.get((inst, local_id))
            if gid is None:
                continue
            bc = bcs.setdefault(gid, [0, 0, 0, 0, 0, 0])
            for dof in range(dof1, dof2 + 1):
                if 1 <= dof <= 6:
                    bc[dof - 1] = 1
    return bcs


def node_rotation_flags(nodes, flat_elements):
    flags = {nid: [1, 1, 1] for nid in nodes}
    for etype in ("S4R", "B31"):
        for _eid, conn, _part_name, _elset in flat_elements[etype]:
            for nid in conn:
                flags[nid] = [0, 0, 0]
    return flags


def material_key(kind, part, elset, materials):
    if kind == "h8":
        sec = part["solid_sections"].get(elset, {})
        mat = materials.get(sec.get("material"), {})
        return (mat.get("E"), mat.get("nu"), mat.get("rho", 0.0))
    if kind == "plate":
        sec = part["shell_sections"].get(elset, {})
        mat = materials.get(sec.get("material"), {})
        return (mat.get("E"), mat.get("nu"), sec.get("thickness", 1.0), mat.get("rho", 0.0))
    if kind == "beam":
        sec = part["beam_sections"].get(elset, {})
        mat = materials.get(sec.get("material"), {})
        return (
            mat.get("E"), mat.get("nu"), sec.get("area", 1.0), sec.get("iy", 1.0),
            sec.get("iz", 1.0), sec.get("j", 1.0), mat.get("rho", 0.0)
        )
    if kind == "truss":
        sec = part["solid_sections"].get(elset, {})
        mat = materials.get(sec.get("material"), {})
        return (mat.get("E"), sec.get("area", 1.0), mat.get("rho", 0.0))
    raise ValueError(kind)


def build_material_map(elements, parts, materials, kind):
    mats = []
    mapping = {}
    for _eid, _conn, part_name, elset in elements:
        key = material_key(kind, parts[part_name], elset, materials)
        if None in key:
            continue
        if key not in mapping:
            mapping[key] = len(mats) + 1
            mats.append(key)
    return mats, mapping


def write_group(f, element_type, elements, mats, mat_map, parts, materials, kind):
    f.write(f"{element_type} {len(elements)} {len(mats)}\n")
    for idx, mat in enumerate(mats, start=1):
        f.write(f"{idx} " + " ".join(str(v) for v in mat) + "\n")

    for idx, (_eid, conn, part_name, elset) in enumerate(elements, start=1):
        key = material_key(kind, parts[part_name], elset, materials)
        mset = mat_map[key]
        if element_type == 5:
            ref = parts[part_name]["beam_sections"].get(elset, {}).get("ref", (0.0, 0.0, 1.0))
            f.write(f"{idx} {conn[0]} {conn[1]} {mset} {ref[0]} {ref[1]} {ref[2]}\n")
        else:
            f.write(f"{idx} " + " ".join(str(n) for n in conn) + f" {mset}\n")


def write_dat(path, parts, instances, materials, assembly_nsets, assembly_surfaces, ties,
              boundaries, loads, gravity, cli_gravity):
    nodes, node_map, flat_elements = flatten_model(parts, instances)
    nodes, node_map, flat_elements, merged_pairs = merge_tied_nodes(
        nodes, node_map, flat_elements, assembly_nsets, assembly_surfaces, ties)
    bcs = apply_boundaries(boundaries, assembly_nsets, node_map)
    rotation_flags = node_rotation_flags(nodes, flat_elements)

    h8 = flat_elements["C3D8R"] + flat_elements["C3D8"]
    plate = flat_elements["S4R"]
    beam = flat_elements["B31"]
    truss = flat_elements["T3D2"]

    group_specs = []
    for element_type, elems, kind in ((4, h8, "h8"), (6, plate, "plate"), (5, beam, "beam"), (1, truss, "truss")):
        if elems:
            mats, mat_map = build_material_map(elems, parts, materials, kind)
            group_specs.append((element_type, elems, mats, mat_map, kind))

    if gravity is None and cli_gravity is not None:
        gravity = cli_gravity
    gravity_flag = 1 if gravity is not None else 0
    if gravity is None:
        gravity = (0.0, 0.0, 0.0)

    mapped_loads = []
    for ref, dof, val in loads:
        if ref in assembly_nsets:
            for inst, local_id in assembly_nsets[ref]:
                gid = node_map.get((inst, local_id))
                if gid is not None:
                    mapped_loads.append((gid, dof, val))

    with open(path, "w", encoding="utf-8") as f:
        f.write("Converted from Abaqus inp\n")
        f.write(f"{len(nodes)} {len(group_specs)} 1 1\n")

        for nid in sorted(nodes):
            bc = [0, 0, 0, *rotation_flags[nid]]
            if nid in bcs:
                for dof, value in enumerate(bcs[nid]):
                    if value:
                        bc[dof] = value
            x, y, z = nodes[nid]
            f.write(f"{nid} " + " ".join(str(v) for v in bc) + f" {x} {y} {z}\n")

        f.write("1\n")
        f.write(f"{len(mapped_loads)} {gravity_flag} {gravity[0]} {gravity[1]} {gravity[2]}\n")
        for nid, dof, val in mapped_loads:
            f.write(f"{nid} {dof} {val}\n")

        for element_type, elems, mats, mat_map, kind in group_specs:
            write_group(f, element_type, elems, mats, mat_map, parts, materials, kind)

    return {
        "nodes": len(nodes),
        "h8": len(h8),
        "plate": len(plate),
        "beam": len(beam),
        "truss": len(truss),
        "materials": len(materials),
        "bcs": sum(sum(bc) for bc in bcs.values()),
        "ties": len(ties),
        "merged_pairs": merged_pairs,
        "load_cases": 1,
        "gravity": gravity_flag,
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
        output = args.inp[:-4] + ".dat" if args.inp.lower().endswith(".inp") else args.inp + ".dat"

    parts, instances, materials, assembly_nsets, assembly_surfaces, ties, boundaries, loads, gravity = parse_inp(args.inp)
    stats = write_dat(output, parts, instances, materials, assembly_nsets, assembly_surfaces, ties,
                      boundaries, loads, gravity, tuple(args.gravity) if args.gravity else None)

    print(f"Nodes: {stats['nodes']}")
    print(f"C3D8R/H8: {stats['h8']}")
    print(f"S4R/Plate4: {stats['plate']}")
    print(f"B31/Beam3D2: {stats['beam']}")
    print(f"T3D2/Truss3D2: {stats['truss']}")
    print(f"Materials: {stats['materials']}")
    print(f"Boundary DOFs: {stats['bcs']}")
    print(f"Tie constraints: {stats['ties']}")
    print(f"Tie node merges: {stats['merged_pairs']}")
    print(f"Load cases: {stats['load_cases']}")
    print(f"Gravity load: {stats['gravity']}")


if __name__ == "__main__":
    main()
