"""Generate and run an Abaqus S4R plate case, then compare with STAP++ Plate4.

Run from the repository root with Abaqus Python:

    abaqus cae noGUI=tools/compare_plate4_s4r.py

Optional arguments after ``--``:

    abaqus cae noGUI=tools/compare_plate4_s4r.py -- --mesh 4 --no-submit

The script creates a 1 x 1 m cantilever plate with S4R elements, gravity load,
and the same material/thickness used by the STAP++ Plate4 convergence case.
It writes an Abaqus inp, optionally submits the Abaqus job, extracts key ODB
results, runs the matching STAP++ dat file if the executable exists, and writes
a CSV comparison table.
"""

from __future__ import print_function

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
WORK_DIR = os.path.join(ROOT, "results", "abaqus_plate4_s4r")

LENGTH = 1.0
WIDTH = 1.0
THICKNESS = 0.1
E = 3.0e10
NU = 0.2
DENSITY = 2500.0
GRAVITY = -9.81


def parse_args():
    argv = sys.argv
    if "--" in argv:
        argv = argv[argv.index("--") + 1:]
    else:
        argv = []

    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=int, default=4, help="number of S4R elements along each side")
    parser.add_argument("--job", default=None, help="Abaqus job name")
    parser.add_argument("--no-submit", action="store_true", help="only write the Abaqus inp")
    parser.add_argument("--no-stappp", action="store_true", help="skip running STAP++")
    parser.add_argument("--stappp-exe", default=os.path.join(ROOT, "build", "Debug", "stap++.exe"))
    return parser.parse_args(argv)


def ensure_work_dir():
    if not os.path.isdir(WORK_DIR):
        os.makedirs(WORK_DIR)


def create_abaqus_model(mesh_n, job_name):
    from abaqus import mdb
    from abaqusConstants import (
        DEFORMABLE_BODY,
        THREE_D,
        OFF,
        ON,
        UNIFORM,
        UNSET,
        STANDARD,
        GRAV,
        S4R,
        MIDDLE_SURFACE,
        FROM_SECTION,
        CARTESIAN,
    )
    import regionToolset
    import mesh

    model_name = "Plate4_S4R_Compare"
    if model_name in mdb.models:
        del mdb.models[model_name]

    model = mdb.Model(name=model_name)

    sketch = model.ConstrainedSketch(name="plate_sketch", sheetSize=2.0)
    sketch.rectangle(point1=(0.0, 0.0), point2=(LENGTH, WIDTH))
    part = model.Part(name="Plate", dimensionality=THREE_D, type=DEFORMABLE_BODY)
    part.BaseShell(sketch=sketch)
    del model.sketches["plate_sketch"]

    face_region = regionToolset.Region(faces=part.faces[:])

    material = model.Material(name="PlateMaterial")
    material.Elastic(table=((E, NU),))
    material.Density(table=((DENSITY,),))

    model.HomogeneousShellSection(
        name="PlateSection",
        material="PlateMaterial",
        thickness=THICKNESS,
        thicknessType=UNIFORM,
    )
    part.SectionAssignment(
        region=face_region,
        sectionName="PlateSection",
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        offsetField="",
        thicknessAssignment=FROM_SECTION,
    )

    part.seedPart(size=LENGTH / float(mesh_n), deviationFactor=0.1, minSizeFactor=0.1)
    part.setElementType(
        regions=face_region,
        elemTypes=(mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD),),
    )
    part.generateMesh()

    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)
    instance = assembly.Instance(name="Plate-1", part=part, dependent=ON)

    fixed_edges = instance.edges.getByBoundingBox(
        xMin=-1.0e-8,
        xMax=1.0e-8,
        yMin=-1.0e-8,
        yMax=WIDTH + 1.0e-8,
        zMin=-1.0e-8,
        zMax=1.0e-8,
    )
    model.EncastreBC(name="ClampX0", createStepName="Initial", region=regionToolset.Region(edges=fixed_edges))

    model.StaticStep(name="Step-1", previous="Initial", nlgeom=OFF)
    model.Gravity(name="Gravity", createStepName="Step-1", comp3=GRAVITY, distributionType=UNIFORM, field="")
    model.fieldOutputRequests["F-Output-1"].setValues(variables=("U", "S", "SF", "SM"))

    job_path = os.path.join(WORK_DIR, job_name)
    job = mdb.Job(name=job_name, model=model_name, description="S4R comparison for STAP++ Plate4")
    job.writeInput(consistencyChecking=OFF)

    inp_from = os.path.join(os.getcwd(), job_name + ".inp")
    inp_to = job_path + ".inp"
    if os.path.abspath(inp_from) != os.path.abspath(inp_to) and os.path.exists(inp_from):
        shutil.move(inp_from, inp_to)

    return job


def submit_abaqus_job(job, job_name):
    old_cwd = os.getcwd()
    os.chdir(WORK_DIR)
    try:
        job.submit(consistencyChecking=False)
        job.waitForCompletion()
    finally:
        os.chdir(old_cwd)


def extract_abaqus_results(job_name):
    from odbAccess import openOdb

    odb_path = os.path.join(WORK_DIR, job_name + ".odb")
    if not os.path.exists(odb_path):
        return {}

    odb = openOdb(odb_path, readOnly=True)
    try:
        frame = odb.steps["Step-1"].frames[-1]
        u_field = frame.fieldOutputs["U"]

        max_u_mag = 0.0
        max_abs_u3 = 0.0
        tip_u3_values = []

        for value in u_field.values:
            u1, u2, u3 = value.data
            mag = math.sqrt(u1 * u1 + u2 * u2 + u3 * u3)
            max_u_mag = max(max_u_mag, mag)
            max_abs_u3 = max(max_abs_u3, abs(u3))

            node = odb.rootAssembly.instances["PLATE-1"].nodes[value.nodeLabel - 1]
            if abs(node.coordinates[0] - LENGTH) < 1.0e-6:
                tip_u3_values.append(u3)

        result = {
            "abaqus_max_u_mag": max_u_mag,
            "abaqus_max_abs_u3": max_abs_u3,
            "abaqus_tip_mean_u3": sum(tip_u3_values) / len(tip_u3_values) if tip_u3_values else None,
        }

        if "SM" in frame.fieldOutputs:
            sm_values = frame.fieldOutputs["SM"].values
            result["abaqus_max_abs_sm1"] = max(abs(v.data[0]) for v in sm_values)
            result["abaqus_max_abs_sm2"] = max(abs(v.data[1]) for v in sm_values)

        return result
    finally:
        odb.close()


def run_stappp(mesh_n, stappp_exe, skip_stappp):
    dat_path = os.path.join(ROOT, "tests", "convergence", "plate4_cantilever_gravity_%dx%d.dat" % (mesh_n, mesh_n))
    if skip_stappp or not os.path.exists(stappp_exe) or not os.path.exists(dat_path):
        return dat_path

    subprocess.check_call([stappp_exe, dat_path], cwd=ROOT)
    return dat_path


def extract_stappp_results(dat_path):
    out_path = os.path.splitext(dat_path)[0] + ".out"
    if not os.path.exists(out_path):
        return {}

    in_displacement = False
    in_stress = False
    max_u_mag = 0.0
    max_abs_u3 = 0.0
    tip_u3_values = []
    max_abs_mx = 0.0
    max_abs_my = 0.0

    node_x = {}
    with open(dat_path, "r") as f:
        lines = f.readlines()
    num_nodes = int(lines[1].split()[0])
    for line in lines[2:2 + num_nodes]:
        parts = line.split()
        node_x[int(parts[0])] = float(parts[7])

    with open(out_path, "r") as f:
        for line in f:
            if "D I S P L A C E M E N T S" in line:
                in_displacement = True
                in_stress = False
                continue
            if "S T R E S S" in line:
                in_displacement = False
                in_stress = True
                continue
            if "S O L U T I O N" in line:
                in_displacement = False
                in_stress = False

            parts = line.split()
            if in_displacement and len(parts) == 7 and parts[0].isdigit():
                node = int(parts[0])
                ux, uy, uz = [float(v) for v in parts[1:4]]
                mag = math.sqrt(ux * ux + uy * uy + uz * uz)
                max_u_mag = max(max_u_mag, mag)
                max_abs_u3 = max(max_abs_u3, abs(uz))
                if abs(node_x.get(node, 0.0) - LENGTH) < 1.0e-8:
                    tip_u3_values.append(uz)

            if in_stress and len(parts) == 8 and parts[0].isdigit():
                mx = float(parts[4])
                my = float(parts[5])
                max_abs_mx = max(max_abs_mx, abs(mx))
                max_abs_my = max(max_abs_my, abs(my))

    return {
        "stappp_max_u_mag": max_u_mag,
        "stappp_max_abs_u3": max_abs_u3,
        "stappp_tip_mean_u3": sum(tip_u3_values) / len(tip_u3_values) if tip_u3_values else None,
        "stappp_max_abs_mx": max_abs_mx,
        "stappp_max_abs_my": max_abs_my,
    }


def rel_error(stappp, abaqus):
    if stappp is None or abaqus is None or abs(abaqus) < 1.0e-30:
        return None
    return abs(stappp - abaqus) / abs(abaqus)


def write_comparison(mesh_n, abaqus_results, stappp_results):
    csv_path = os.path.join(WORK_DIR, "plate4_s4r_compare_%dx%d.csv" % (mesh_n, mesh_n))
    rows = [
        ("max_u_mag", stappp_results.get("stappp_max_u_mag"), abaqus_results.get("abaqus_max_u_mag")),
        ("max_abs_u3", stappp_results.get("stappp_max_abs_u3"), abaqus_results.get("abaqus_max_abs_u3")),
        ("tip_mean_u3", stappp_results.get("stappp_tip_mean_u3"), abaqus_results.get("abaqus_tip_mean_u3")),
        ("max_abs_mx_or_sm1", stappp_results.get("stappp_max_abs_mx"), abaqus_results.get("abaqus_max_abs_sm1")),
        ("max_abs_my_or_sm2", stappp_results.get("stappp_max_abs_my"), abaqus_results.get("abaqus_max_abs_sm2")),
    ]

    with open(csv_path, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["quantity", "stappp", "abaqus", "relative_error"])
        for name, stappp, abaqus in rows:
            err = rel_error(stappp, abaqus)
            writer.writerow([
                name,
                "" if stappp is None else "%.12e" % stappp,
                "" if abaqus is None else "%.12e" % abaqus,
                "" if err is None else "%.6e" % err,
            ])

    return csv_path


def main():
    args = parse_args()
    ensure_work_dir()
    job_name = args.job or "plate4_s4r_%dx%d" % (args.mesh, args.mesh)

    old_cwd = os.getcwd()
    os.chdir(WORK_DIR)
    try:
        job = create_abaqus_model(args.mesh, job_name)
        inp_path = os.path.join(WORK_DIR, job_name + ".inp")
        print("Wrote Abaqus inp:", inp_path)

        if not args.no_submit:
            submit_abaqus_job(job, job_name)
            print("Finished Abaqus job:", job_name)
    finally:
        os.chdir(old_cwd)

    abaqus_results = {} if args.no_submit else extract_abaqus_results(job_name)
    dat_path = run_stappp(args.mesh, args.stappp_exe, args.no_stappp)
    stappp_results = extract_stappp_results(dat_path)
    csv_path = write_comparison(args.mesh, abaqus_results, stappp_results)

    print("STAP++ dat:", dat_path)
    print("Comparison CSV:", csv_path)
    print("Abaqus results:", abaqus_results)
    print("STAP++ results:", stappp_results)


if __name__ == "__main__":
    main()
