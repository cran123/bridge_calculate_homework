#!/usr/bin/env python
"""Generate compact Abaqus input files from the course-design STAP++ cases."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


@dataclass
class Node:
    label: int
    bc: list[int]
    xyz: list[float]


@dataclass
class Load:
    node: int
    dof: int
    value: float


@dataclass
class MpcTerm:
    node: int
    dof: int
    coefficient: float


@dataclass
class Mpc:
    terms: list[MpcTerm]


@dataclass
class Group:
    etype: int
    materials: list[list[float]]
    elements: list[list[float]]


def clean_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip() and not line.strip().startswith("#")]


def parse_dat(path: Path):
    lines = clean_lines(path)
    title = lines[0]
    nnode, ngroup, ncase, _ = map(int, lines[1].split()[:4])
    idx = 2
    nodes = []
    for _ in range(nnode):
        parts = lines[idx].split()
        nodes.append(Node(int(parts[0]), [int(x) for x in parts[1:7]], [float(x) for x in parts[7:10]]))
        idx += 1
    loads: list[Load] = []
    for _ in range(ncase):
        idx += 1
        header = lines[idx].split()
        nloads = int(header[0])
        idx += 1
        for _ in range(nloads):
            node, dof, value = lines[idx].split()[:3]
            loads.append(Load(int(node), int(dof), float(value)))
            idx += 1
        if len(header) >= 6:
            idx += int(header[5])
    mpcs: list[Mpc] = []
    if idx < len(lines):
        maybe_mpc_or_group = lines[idx].split()
        if len(maybe_mpc_or_group) == 1:
            nmpc = int(float(maybe_mpc_or_group[0]))
            idx += 1
            for _ in range(nmpc):
                nterms = int(float(lines[idx].split()[0]))
                idx += 1
                terms = []
                for _ in range(nterms):
                    node, dof, coefficient = lines[idx].split()[:3]
                    terms.append(MpcTerm(int(node), int(dof), float(coefficient)))
                    idx += 1
                mpcs.append(Mpc(terms))
    groups = []
    for _ in range(ngroup):
        etype, nelem, nmat = map(int, lines[idx].split()[:3])
        idx += 1
        mats = []
        for _ in range(nmat):
            mats.append([float(x) for x in lines[idx].split()])
            idx += 1
        elems = []
        for _ in range(nelem):
            elems.append([float(x) for x in lines[idx].split()])
            idx += 1
        groups.append(Group(etype, mats, elems))
    return title, nodes, loads, mpcs, groups


def abaqus_element_type(etype: int, mats: list[list[float]]) -> str:
    if etype == 1:
        return "T3D2"
    if etype == 4:
        return "C3D8"
    if etype == 5:
        return "B31"
    if etype == 6:
        return "S4R"
    if etype == 8:
        nu = mats[0][2] if mats and len(mats[0]) > 2 else 0.3
        return "C3D8H" if nu > 0.45 else "C3D8"
    raise ValueError(f"Unsupported element type {etype}")


def mat_name(group_index: int, set_id: int) -> str:
    return f"MAT_G{group_index}_{set_id}"


def write_inp(dat: Path, target: Path) -> None:
    title, nodes, loads, mpcs, groups = parse_dat(dat)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="ascii", newline="\n") as out:
        out.write(f"*Heading\n{title} - Abaqus reference generated from STAP++ DAT\n")
        out.write("*Preprint, echo=NO, model=NO, history=NO, contact=NO\n")
        out.write("*Node\n")
        for node in nodes:
            out.write(f"{node.label}, {node.xyz[0]:.12g}, {node.xyz[1]:.12g}, {node.xyz[2]:.12g}\n")

        global_element = 1
        section_records = []
        for gi, group in enumerate(groups, start=1):
            abq_type = abaqus_element_type(group.etype, group.materials)
            out.write(f"*Element, type={abq_type}, elset=G{gi}\n")
            for elem in group.elements:
                conn_count = {1: 2, 4: 8, 5: 2, 6: 4, 8: 8}[group.etype]
                conn = [int(x) for x in elem[1 : 1 + conn_count]]
                out.write(f"{global_element}, " + ", ".join(str(x) for x in conn) + "\n")
                mat = int(elem[1 + conn_count])
                section_records.append((gi, global_element, group.etype, mat, elem))
                global_element += 1

        for gi, group in enumerate(groups, start=1):
            for mat in group.materials:
                set_id = int(mat[0])
                e = mat[1]
                nu = 0.3
                if group.etype in (4, 6, 8) and len(mat) > 2:
                    nu = mat[2]
                if group.etype == 5 and len(mat) > 2:
                    nu = mat[2]
                out.write(f"*Material, name={mat_name(gi, set_id)}\n")
                out.write("*Elastic\n")
                out.write(f"{e:.12g}, {nu:.12g}\n")

        for gi, group in enumerate(groups, start=1):
            for mat in group.materials:
                set_id = int(mat[0])
                labels = [str(label) for g, label, _, m, _ in section_records if g == gi and m == set_id]
                if not labels:
                    continue
                elset = f"G{gi}_M{set_id}"
                out.write(f"*Elset, elset={elset}\n")
                out.write(", ".join(labels) + "\n")
                material = mat_name(gi, set_id)
                if group.etype == 1:
                    area = mat[2]
                    out.write(f"*Solid Section, elset={elset}, material={material}\n")
                    out.write(f"{area:.12g}\n")
                elif group.etype in (4, 8):
                    out.write(f"*Solid Section, elset={elset}, material={material}\n")
                elif group.etype == 6:
                    thickness = mat[3]
                    out.write(f"*Shell Section, elset={elset}, material={material}\n")
                    out.write(f"{thickness:.12g}\n")
                elif group.etype == 5:
                    area, iy, iz, j = mat[3], mat[4], mat[5], mat[6]
                    e = mat[1]
                    nu = mat[2]
                    g = e / (2.0 * (1.0 + nu))
                    out.write(f"*Beam General Section, elset={elset}, section=GENERAL\n")
                    out.write(f"{area:.12g}, {iy:.12g}, 0., {iz:.12g}, {j:.12g}\n")
                    ref = next((rec[4] for rec in section_records if rec[0] == gi and rec[3] == set_id), None)
                    if ref is not None and len(ref) >= 7:
                        out.write(f"{ref[-3]:.12g}, {ref[-2]:.12g}, {ref[-1]:.12g}\n")
                    else:
                        out.write("0., 0., 1.\n")
                    out.write(f"{e:.12g}, {g:.12g}\n")

        for mpc in mpcs:
            out.write("*Equation\n")
            out.write(f"{len(mpc.terms)}\n")
            fields = []
            for term in mpc.terms:
                fields.extend([str(term.node), str(term.dof), f"{term.coefficient:.12g}"])
            out.write(", ".join(fields) + "\n")

        out.write("*Boundary\n")
        for node in nodes:
            for dof, code in enumerate(node.bc, start=1):
                if code:
                    out.write(f"{node.label}, {dof}, {dof}, 0.\n")
        out.write("*Step, name=Load, nlgeom=NO\n")
        out.write("*Static\n1., 1., 1e-05, 1.\n")
        if loads:
            out.write("*Cload\n")
            for load in loads:
                out.write(f"{load.node}, {load.dof}, {load.value:.12g}\n")
        out.write("*Output, field\n")
        out.write("*Node Output\nU, UR\n")
        out.write("*Element Output, directions=YES\nS, SF, SM\n")
        out.write("*End Step\n")


def main() -> None:
    generated = []
    for dat in sorted((ROOT / "cases").glob("*/input/*.dat")):
        target = dat.parents[1] / "abaqus" / (dat.stem + "_abaqus.inp")
        write_inp(dat, target)
        generated.append(target)
    for path in generated:
        print(path)


if __name__ == "__main__":
    main()
