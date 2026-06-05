"""Create and run Abaqus reference jobs for the STAP++ course-design cases.

This script is intended for:

    abaqus cae noGUI=abaqus_scripts/run_validation_jobs.py -- --repo-root <repo>
"""

from __future__ import annotations

import argparse
import os
import sys

from abaqus import mdb
from abaqusConstants import (  # type: ignore
    DEFORMABLE_BODY,
    THREE_D,
    ON,
    OFF,
    CARTESIAN,
    UNIFORM,
    NODAL,
    MIDDLE_SURFACE,
    FROM_SECTION,
    DURING_ANALYSIS,
    DEFAULT,
    SINGLE,
    PERCENTAGE,
    ANALYSIS,
    C3D8,
    B31,
    S4R,
    T3D2,
    N1_COSINES,
    UNSET,
)
import regionToolset  # type: ignore


def parse_args() -> argparse.Namespace:
    argv = sys.argv
    if "--" in argv:
        argv = argv[argv.index("--") + 1 :]
    else:
        argv = []
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", default=os.getcwd())
    parser.add_argument("--submit", action="store_true", default=True)
    return parser.parse_args(argv)


def job_dir(repo_root: str, name: str) -> str:
    path = os.path.join(repo_root, "cases", name, "abaqus")
    if not os.path.isdir(path):
        os.makedirs(path)
    return path


def make_part_from_nodes(model, name, nodes, elements, elem_type):
    part = model.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    for label, xyz in nodes:
        part.Node(coordinates=xyz, label=label)
    for label, conn in elements:
        part.Element(nodes=conn, label=label, elemShape=None)
    return part


def assign_section(part, section_name, elem_code):
    region = regionToolset.Region(elements=part.elements)
    if elem_code == "bar":
        part.SectionAssignment(region=region, sectionName=section_name)
        part.setElementType(regions=region, elemTypes=(T3D2,))
    elif elem_code == "beam":
        part.SectionAssignment(region=region, sectionName=section_name)
        part.assignBeamSectionOrientation(region=region, method=N1_COSINES, n1=(0.0, 0.0, 1.0))
        part.setElementType(regions=region, elemTypes=(B31,))
    elif elem_code == "plate":
        part.SectionAssignment(region=region, sectionName=section_name, offset=0.0,
                               offsetType=MIDDLE_SURFACE, offsetField="",
                               thicknessAssignment=FROM_SECTION)
        part.setElementType(regions=region, elemTypes=(S4R,))
    else:
        part.SectionAssignment(region=region, sectionName=section_name)
        part.setElementType(regions=region, elemTypes=(C3D8,))


def create_bar(repo_root: str):
    model = mdb.Model(name="bar_patch")
    model.Material(name="steel").Elastic(table=((1.0e6, 0.3),))
    model.TrussSection(name="bar_section", material="steel", area=1.0)
    part = make_part_from_nodes(
        model,
        "bar",
        [(1, (0.0, 0.0, 0.0)), (2, (1.0, 0.0, 0.0))],
        [(1, (1, 2))],
        "bar",
    )
    assign_section(part, "bar_section", "bar")
    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)
    inst = assembly.Instance(name="bar-1", part=part, dependent=ON)
    model.StaticStep(name="load", previous="Initial")
    model.EncastreBC(name="fixed", createStepName="Initial", region=regionToolset.Region(nodes=inst.nodes.sequenceFromLabels((1,))))
    model.ConcentratedForce(name="tip", createStepName="load", region=regionToolset.Region(nodes=inst.nodes.sequenceFromLabels((2,))), cf1=1.0)
    return write_and_submit(repo_root, "bar_patch", model)


def create_hex(repo_root: str):
    model = mdb.Model(name="hex8_patch")
    model.Material(name="mat").Elastic(table=((1.0e3, 0.3),))
    model.HomogeneousSolidSection(name="solid", material="mat")
    nodes = [
        (1, (0, 0, 0)), (2, (1, 0, 0)), (3, (1, 1, 0)), (4, (0, 1, 0)),
        (5, (0, 0, 1)), (6, (1, 0, 1)), (7, (1, 1, 1)), (8, (0, 1, 1)),
    ]
    part = make_part_from_nodes(model, "hex", nodes, [(1, (1, 2, 3, 4, 5, 6, 7, 8))], "hex")
    assign_section(part, "solid", "hex")
    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)
    inst = assembly.Instance(name="hex-1", part=part, dependent=ON)
    model.StaticStep(name="load", previous="Initial")
    model.DisplacementBC(name="left", createStepName="Initial",
                         region=regionToolset.Region(nodes=inst.nodes.sequenceFromLabels((1, 4, 5, 8))),
                         u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    for label in (2, 3, 6, 7):
        model.ConcentratedForce(name="fx_%s" % label, createStepName="load",
                                region=regionToolset.Region(nodes=inst.nodes.sequenceFromLabels((label,))), cf1=0.25)
    return write_and_submit(repo_root, "hex8_patch", model)


def write_and_submit(repo_root: str, case_name: str, model):
    work = job_dir(repo_root, case_name)
    old = os.getcwd()
    os.chdir(work)
    try:
        job_name = case_name + "_abaqus"
        job = mdb.Job(name=job_name, model=model.name, type=ANALYSIS,
                      multiprocessingMode=DEFAULT, numCpus=1, numDomains=1,
                      memory=90, memoryUnits=PERCENTAGE, explicitPrecision=SINGLE,
                      nodalOutputPrecision=SINGLE)
        job.writeInput()
        job.submit(consistencyChecking=OFF)
        job.waitForCompletion()
    finally:
        os.chdir(old)
    return os.path.join(work, job_name + ".odb")


def main() -> None:
    args = parse_args()
    repo_root = os.path.abspath(args.repo_root)
    created = []
    created.append(create_bar(repo_root))
    created.append(create_hex(repo_root))
    extractor = os.path.join(repo_root, "abaqus_scripts", "extract_odb_results.py")
    for odb in created:
        os.system('abaqus python "%s" --odb "%s"' % (extractor, odb))
    print("Created Abaqus validation jobs:")
    for path in created:
        print(path)


if __name__ == "__main__":
    main()
