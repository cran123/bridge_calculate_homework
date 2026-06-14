"""Extract node coordinates + displacements from ODB."""
import csv
import sys
from odbAccess import openOdb


def main():
    odb_path = sys.argv[sys.argv.index("--odb") + 1] if "--odb" in sys.argv else None
    if not odb_path:
        print("Usage: abaqus python extract_odb_coords.py --odb model.odb")
        return

    odb = openOdb(odb_path, readOnly=True)
    step = odb.steps[list(odb.steps.keys())[-1]]
    frame = step.frames[-1]
    
    # Get instance node coordinates
    assembly = odb.rootAssembly
    base = odb_path.replace(".odb", "")

    with open(base + "_nodes_full.csv", "wb") as f:
        writer = csv.writer(f)
        writer.writerow(["instance", "node", "x", "y", "z", "u1", "u2", "u3"])
        disp = frame.fieldOutputs["U"]
        for value in disp.values:
            inst = value.instance
            node_label = value.nodeLabel
            # Get coordinates
            try:
                node_obj = inst.getNodeFromLabel(node_label)
                coord = node_obj.coordinates
            except:
                coord = [0.0, 0.0, 0.0]
            data = list(value.data)
            while len(data) < 3:
                data.append(0.0)
            writer.writerow([inst.name, node_label, coord[0], coord[1], coord[2]] + data[:3])

    odb.close()
    print("Done: " + base + "_nodes_full.csv")


if __name__ == "__main__":
    main()
