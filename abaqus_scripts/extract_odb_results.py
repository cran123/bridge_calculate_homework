import csv
import os
import sys

from odbAccess import openOdb  # type: ignore


def parse_odb_arg(argv):
    if "--odb" in argv:
        idx = argv.index("--odb")
        if idx + 1 < len(argv):
            return argv[idx + 1]
    raise SystemExit("Usage: abaqus python extract_odb_results.py --odb model.odb")


def main():
    odb_path = parse_odb_arg(sys.argv[1:])
    odb = openOdb(odb_path, readOnly=True)
    step = odb.steps[list(odb.steps.keys())[-1]]
    frame = step.frames[-1]
    base = os.path.splitext(odb_path)[0]

    with open(base + "_nodes.csv", "wb") as stream:
        writer = csv.writer(stream)
        writer.writerow(["instance", "node", "u1", "u2", "u3"])
        disp = frame.fieldOutputs["U"]
        for value in disp.values:
            data = list(value.data)
            while len(data) < 3:
                data.append(0.0)
            writer.writerow([value.instance.name, value.nodeLabel] + data[:3])

    if "S" in frame.fieldOutputs:
        with open(base + "_stress.csv", "wb") as stream:
            writer = csv.writer(stream)
            writer.writerow(["instance", "element", "s11", "s22", "s33", "s12", "s13", "s23", "mises"])
            for value in frame.fieldOutputs["S"].values:
                data = list(value.data)
                while len(data) < 6:
                    data.append(0.0)
                writer.writerow([value.instance.name, value.elementLabel] + data[:6] + [value.mises])

    odb.close()
    print(base + "_nodes.csv")


if __name__ == "__main__":
    main()
