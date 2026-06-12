#!/usr/bin/env python3
"""Fast wrapper for the original Abaqus-to-STAP++ converter.

The original converter is preserved as ``tools/abaqus_to_stappp.py``.  This
wrapper imports the implementation from ``others/tools/abaqus_to_stappp.py`` and
patches the slow tie-MPC path by caching master-surface facets.  The command
line interface and output DAT format are intentionally kept compatible.
"""

from pathlib import Path
import runpy


ORIGINAL = Path(__file__).resolve().parents[1] / "others" / "tools" / "abaqus_to_stappp.py"


def install_fast_tie_builder(module):
    surface_nodes = module["surface_nodes"]
    node_element_kinds = module["node_element_kinds"]
    add_mpc_term = module["add_mpc_term"]
    compact_mpc_terms = module["compact_mpc_terms"]
    project_to_quad = module["project_to_quad"]
    project_to_line = module["project_to_line"]
    inverse_distance_weights = module["inverse_distance_weights"]
    nearest_node = module["nearest_node"]

    def facets_for_master_surface(master_key, master_nodes, flat_elements):
        master_set = set(master_nodes)
        facets = []
        h8_faces = (
            (0, 1, 2, 3),
            (4, 5, 6, 7),
            (0, 1, 5, 4),
            (1, 2, 6, 5),
            (2, 3, 7, 6),
            (3, 0, 4, 7),
        )

        for etype in ("C3D8R", "C3D8"):
            for _eid, conn, _part_name, _elset, *_metadata in flat_elements[etype]:
                for face in h8_faces:
                    face_conn = tuple(conn[i] for i in face)
                    if all(nid in master_set for nid in face_conn):
                        facets.append(("quad", face_conn))

        for _eid, conn, _part_name, _elset, *_metadata in flat_elements["S4R"]:
            if all(nid in master_set for nid in conn):
                facets.append(("quad", tuple(conn)))

        for _eid, conn, _part_name, _elset, *_metadata in flat_elements["B31"]:
            if all(nid in master_set for nid in conn):
                facets.append(("line", tuple(conn)))

        return facets

    def cached_master_surface_weights(source, master_nodes, facets, nodes, node_based=False):
        if node_based:
            master = nearest_node(source, master_nodes, nodes)
            return [(master, 1.0)] if master is not None else []

        best = None
        for kind, conn in facets:
            projected = project_to_quad(source, conn, nodes) if kind == "quad" else project_to_line(source, conn, nodes)
            if projected is None:
                continue
            distance2, weights = projected
            if best is None or distance2 < best[0]:
                best = (distance2, weights)
        if best is not None:
            return best[1]
        return inverse_distance_weights(source, master_nodes, nodes)

    def build_tie_mpcs_fast(nodes, node_map, flat_elements, assembly_nsets, assembly_surfaces,
                            assembly_surface_types, ties):
        kinds = node_element_kinds(flat_elements)
        mpcs = []
        generated = 0
        cached = {}

        for tie in ties:
            slave_nodes = surface_nodes(tie.slave_surface, assembly_nsets, assembly_surfaces, node_map)
            master_nodes = surface_nodes(tie.master_surface, assembly_nsets, assembly_surfaces, node_map)
            if not slave_nodes or not master_nodes:
                continue

            node_based_master = assembly_surface_types.get(tie.master_surface, "").upper() == "NODE"
            cache_key = (tie.master_surface, tuple(master_nodes))
            if cache_key not in cached:
                cached[cache_key] = facets_for_master_surface(cache_key, master_nodes, flat_elements)
            facets = cached[cache_key]

            for slave in slave_nodes:
                weights = cached_master_surface_weights(slave, master_nodes, facets, nodes, node_based=node_based_master)
                sx, sy, sz = nodes[slave]
                cx = sum(nodes[nid][0] * weight for nid, weight in weights)
                cy = sum(nodes[nid][1] * weight for nid, weight in weights)
                cz = sum(nodes[nid][2] * weight for nid, weight in weights)
                rx, ry, rz = sx - cx, sy - cy, sz - cz

                slave_has_rot = kinds.get(slave, {}).get("rotation", False)
                master_has_rot = bool(weights) and all(kinds.get(nid, {}).get("rotation", False) for nid, _w in weights)
                master_is_solid = bool(weights) and any(kinds.get(nid, {}).get("solid", False) for nid, _w in weights)

                for comp in range(3):
                    terms = {}
                    add_mpc_term(terms, slave, comp + 1, 1.0)
                    for master, weight in weights:
                        add_mpc_term(terms, master, comp + 1, -weight)

                    if not tie.no_rotation and slave_has_rot and master_is_solid:
                        if comp == 0:
                            add_mpc_term(terms, slave, 5, rz)
                            add_mpc_term(terms, slave, 6, -ry)
                        elif comp == 1:
                            add_mpc_term(terms, slave, 4, -rz)
                            add_mpc_term(terms, slave, 6, rx)
                        else:
                            add_mpc_term(terms, slave, 4, ry)
                            add_mpc_term(terms, slave, 5, -rx)

                    compact = compact_mpc_terms(terms)
                    if len(compact) >= 2:
                        mpcs.append(compact)
                        generated += 1

                if not tie.no_rotation and slave_has_rot and master_has_rot:
                    for dof in range(4, 7):
                        terms = {}
                        add_mpc_term(terms, slave, dof, 1.0)
                        for master, weight in weights:
                            add_mpc_term(terms, master, dof, -weight)
                        compact = compact_mpc_terms(terms)
                        if len(compact) >= 2:
                            mpcs.append(compact)
                            generated += 1

        return mpcs, generated

    module["build_tie_mpcs"] = build_tie_mpcs_fast


def main():
    module = runpy.run_path(str(ORIGINAL), run_name="abaqus_to_stappp_original")
    install_fast_tie_builder(module)
    module["main"]()


if __name__ == "__main__":
    main()
