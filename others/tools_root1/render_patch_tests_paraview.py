"""Render STAP++ patch-test VTU files with ParaView Python.

Run with ParaView's Python, for example:

    pvpython tools/render_patch_tests_paraview.py

The script writes PNG screenshots under ``cases/patch_test_paraview``.
"""

from pathlib import Path
from xml.sax.saxutils import escape

from paraview.simple import (  # type: ignore
    ColorBy,
    CreateRenderView,
    Delete,
    GetColorTransferFunction,
    GetLayout,
    Hide,
    HideScalarBarIfNotNeeded,
    ResetSession,
    SaveScreenshot,
    SetActiveSource,
    Show,
    Text,
    XMLUnstructuredGridReader,
    WarpByVector,
)


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "cases" / "patch_test_paraview"
VTU_DIR = OUT_DIR / "subdivided_vtu"

MODELS = [
    {
        "key": "bar",
        "prefix": "01_bar",
        "label": "bar axial patch",
        "file": VTU_DIR / "bar_patch_6_segments.vtu",
        "ux_range": (0.0, 1.0e-6),
        "scale": 200000.0,
        "camera": "xy",
        "line_width": 12.0,
        "note": "Six bar elements; UX is linear from 0 to 1e-6.",
    },
    {
        "key": "plate",
        "prefix": "02_plate",
        "label": "Plate4 membrane patch",
        "file": VTU_DIR / "plate4_patch_3x3.vtu",
        "ux_range": (0.0, 3.2e-6),
        "scale": 80000.0,
        "camera": "xy",
        "line_width": 2.0,
        "note": "3x3 Plate4 elements; UX is linear across the patch.",
    },
    {
        "key": "beam",
        "prefix": "03_beam",
        "label": "Beam3D2 axial patch",
        "file": VTU_DIR / "beam_patch_6_segments.vtu",
        "ux_range": (0.0, 1.0e-6),
        "scale": 200000.0,
        "camera": "xy",
        "line_width": 12.0,
        "note": "Six Beam3D2 elements; axial UX is linear across the beam patch.",
    },
    {
        "key": "hex",
        "prefix": "04_hex8",
        "label": "Hex8 uniaxial solid patch",
        "file": VTU_DIR / "hex8_patch_2x2x2.vtu",
        "ux_range": (0.0, 1.0e-3),
        "scale": 120.0,
        "camera": "iso",
        "line_width": 2.0,
        "note": "2x2x2 Hex8 elements; UX is linear with x.",
    },
    {
        "key": "bbarhex8",
        "prefix": "05_bbarhex8",
        "label": "BbarHex8 uniaxial solid patch",
        "file": VTU_DIR / "bbarhex8_patch_2x2x2.vtu",
        "ux_range": (0.0, 1.0e-3),
        "scale": 120.0,
        "camera": "iso",
        "line_width": 2.0,
        "note": "2x2x2 BbarHex8 elements; UX is linear with x.",
    },
]

PATCH_FORMS = [
    {
        "key": "all_node_displacement",
        "index": "type1",
        "title": "Type 1 - prescribed displacement on every node",
        "marker": "AllNodeDisplacement",
        "note": "All nodes are constrained to the target affine displacement field.",
    },
    {
        "key": "boundary_displacement",
        "index": "type2",
        "title": "Type 2 - prescribed displacement on boundary nodes",
        "marker": "BoundaryDisplacement",
        "note": "Only boundary nodes are constrained; interior nodes should recover the same affine field.",
    },
    {
        "key": "equivalent_force_boundary",
        "index": "type3",
        "title": "Type 3 - equivalent force boundary",
        "marker": "EquivalentForceNode",
        "note": "Loaded boundary nodes carry equivalent nodal forces for the target constant stress state.",
    },
]

PATCH_TEST_TYPES = [
    "1. Prescribed displacement on every node: directly impose the target affine displacement field.",
    "2. Prescribed displacement on boundary nodes: constrain only the patch boundary to the affine field.",
    "3. Equivalent force boundary: apply nodal forces equivalent to the same constant stress state.",
]


def fmt(values):
    return " ".join(f"{value:.12e}" for value in values)


def write_vtu(path, title, points, cells, cell_type, displacements, point_data, cell_data):
    path.parent.mkdir(parents=True, exist_ok=True)
    magnitudes = [
        (ux * ux + uy * uy + uz * uz) ** 0.5 for ux, uy, uz in displacements
    ]
    offsets = []
    total = 0
    for cell in cells:
        total += len(cell)
        offsets.append(total)

    lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">',
        f"  <!-- {escape(title)} -->",
        "  <UnstructuredGrid>",
        f'    <Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(cells)}">',
        "      <Points>",
        '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">',
    ]
    lines.extend(f"          {fmt(point)}" for point in points)
    lines.extend(
        [
            "        </DataArray>",
            "      </Points>",
            "      <Cells>",
            '        <DataArray type="Int32" Name="connectivity" format="ascii">',
            "          " + " ".join(str(node) for cell in cells for node in cell),
            "        </DataArray>",
            '        <DataArray type="Int32" Name="offsets" format="ascii">',
            "          " + " ".join(str(offset) for offset in offsets),
            "        </DataArray>",
            '        <DataArray type="UInt8" Name="types" format="ascii">',
            "          " + " ".join(str(cell_type) for _ in cells),
            "        </DataArray>",
            "      </Cells>",
            "      <PointData>",
            '        <DataArray type="Float64" Name="Displacement" NumberOfComponents="3" format="ascii">',
        ]
    )
    lines.extend(f"          {fmt(disp)}" for disp in displacements)
    lines.extend(
        [
            "        </DataArray>",
            '        <DataArray type="Float64" Name="DisplacementMagnitude" format="ascii">',
        ]
    )
    lines.extend(f"          {value:.12e}" for value in magnitudes)
    lines.append("        </DataArray>")
    component_data = {
        "UX": [disp[0] for disp in displacements],
        "UY": [disp[1] for disp in displacements],
        "UZ": [disp[2] for disp in displacements],
    }
    component_data.update(point_data)
    for name, values in component_data.items():
        lines.append(f'        <DataArray type="Float64" Name="{name}" format="ascii">')
        lines.extend(f"          {value:.12e}" for value in values)
        lines.append("        </DataArray>")
    lines.extend(["      </PointData>", "      <CellData>"])
    for name, values in cell_data.items():
        lines.append(f'        <DataArray type="Float64" Name="{name}" format="ascii">')
        lines.extend(f"          {value:.12e}" for value in values)
        lines.append("        </DataArray>")
    lines.extend(
        [
            "      </CellData>",
            "    </Piece>",
            "  </UnstructuredGrid>",
            "</VTKFile>",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def write_marker_vtu(path, title, points, displacements, flags):
    selected_points = []
    selected_displacements = []
    for point, displacement, flag in zip(points, displacements, flags):
        if flag > 0.5:
            selected_points.append(point)
            selected_displacements.append(displacement)
    cells = [(i,) for i in range(len(selected_points))]
    write_vtu(
        path,
        title,
        selected_points,
        cells,
        1,
        selected_displacements,
        {"Marker": [1.0] * len(selected_points)},
        {"MarkerCell": [1.0] * len(selected_points)},
    )


def boundary_flags(points, dimensions):
    active = dimensions
    flags = []
    for x, y, z in points:
        checks = []
        if "x" in active:
            checks.append(abs(x) < 1.0e-12 or abs(x - 1.0) < 1.0e-12)
        if "y" in active:
            checks.append(abs(y) < 1.0e-12 or abs(y - 1.0) < 1.0e-12)
        if "z" in active:
            checks.append(abs(z) < 1.0e-12 or abs(z - 1.0) < 1.0e-12)
        on_boundary = any(checks)
        flags.append(1.0 if on_boundary else 0.0)
    return flags


def force_flags(points, dimensions):
    flags = []
    for x, y, z in points:
        on_loaded_face = abs(x - 1.0) < 1.0e-12
        if "y" in dimensions:
            on_loaded_face = on_loaded_face and (
                abs(y) < 1.0e-12 or abs(y - 1.0) < 1.0e-12
            )
        if "z" in dimensions:
            on_loaded_face = on_loaded_face and (
                abs(z) < 1.0e-12 or abs(z - 1.0) < 1.0e-12
            )
        flags.append(1.0 if on_loaded_face else 0.0)
    return flags


def patch_point_data(points, dimensions):
    boundary = boundary_flags(points, dimensions)
    return {
        "AllNodeDisplacement": [1.0] * len(points),
        "BoundaryDisplacement": boundary,
        "EquivalentForceNode": force_flags(points, dimensions),
    }


def generate_bar_patch():
    nx = 6
    points = [(i / nx, 0.0, 0.0) for i in range(nx + 1)]
    cells = [(i, i + 1) for i in range(nx)]
    displacements = [(1.0e-6 * x, 0.0, 0.0) for x, _, _ in points]
    markers = patch_point_data(points, "x")
    write_vtu(
        VTU_DIR / "bar_patch_6_segments.vtu",
        "Six-element bar patch with an affine axial displacement field",
        points,
        cells,
        3,
        displacements,
        markers,
        {
            "PatchElementId": [float(i + 1) for i in range(nx)],
            "AxialForce": [1.0] * nx,
            "SXX": [1.0] * nx,
        },
    )
    for form in PATCH_FORMS:
        write_marker_vtu(
            VTU_DIR / f"bar_patch_6_segments_{form['key']}_markers.vtu",
            f"Bar patch markers for {form['title']}",
            points,
            displacements,
            markers[form["marker"]],
        )


def generate_beam_patch():
    nx = 6
    points = [(i / nx, 0.12, 0.0) for i in range(nx + 1)]
    cells = [(i, i + 1) for i in range(nx)]
    displacements = [(1.0e-6 * x, 0.0, 0.0) for x, _, _ in points]
    markers = patch_point_data(points, "x")
    write_vtu(
        VTU_DIR / "beam_patch_6_segments.vtu",
        "Six-element Beam3D2 axial patch with an affine axial displacement field",
        points,
        cells,
        3,
        displacements,
        markers,
        {
            "PatchElementId": [float(i + 1) for i in range(nx)],
            "AxialForce": [1.0] * nx,
            "BeamMy": [0.0] * nx,
            "BeamMz": [0.0] * nx,
            "BeamTorque": [0.0] * nx,
        },
    )
    for form in PATCH_FORMS:
        write_marker_vtu(
            VTU_DIR / f"beam_patch_6_segments_{form['key']}_markers.vtu",
            f"Beam patch markers for {form['title']}",
            points,
            displacements,
            markers[form["marker"]],
        )


def generate_plate_patch():
    nx = ny = 3

    def node(i, j):
        return j * (nx + 1) + i

    points = [(i / nx, j / ny, 0.0) for j in range(ny + 1) for i in range(nx + 1)]
    cells = [
        (node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1))
        for j in range(ny)
        for i in range(nx)
    ]
    displacements = [(3.2e-6 * x, 0.0, 0.0) for x, _, _ in points]
    markers = patch_point_data(points, "xy")
    write_vtu(
        VTU_DIR / "plate4_patch_3x3.vtu",
        "Nine-element Plate4 membrane patch with an affine displacement field",
        points,
        cells,
        9,
        displacements,
        markers,
        {
            "PatchElementId": [float(i + 1) for i in range(len(cells))],
            "SXX": [1.0e5] * len(cells),
            "SYY": [2.0e4] * len(cells),
            "SXY": [0.0] * len(cells),
        },
    )
    for form in PATCH_FORMS:
        write_marker_vtu(
            VTU_DIR / f"plate4_patch_3x3_{form['key']}_markers.vtu",
            f"Plate4 patch markers for {form['title']}",
            points,
            displacements,
            markers[form["marker"]],
        )


def generate_hex_patch():
    nx = ny = nz = 2

    def node(i, j, k):
        return k * (ny + 1) * (nx + 1) + j * (nx + 1) + i

    points = [
        (i / nx, j / ny, k / nz)
        for k in range(nz + 1)
        for j in range(ny + 1)
        for i in range(nx + 1)
    ]
    cells = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                cells.append(
                    (
                        node(i, j, k),
                        node(i + 1, j, k),
                        node(i + 1, j + 1, k),
                        node(i, j + 1, k),
                        node(i, j, k + 1),
                        node(i + 1, j, k + 1),
                        node(i + 1, j + 1, k + 1),
                        node(i, j + 1, k + 1),
                    )
                )
    displacements = [
        (1.0e-3 * x, -3.0e-4 * y, -3.0e-4 * z) for x, y, z in points
    ]
    markers = patch_point_data(points, "xyz")
    write_vtu(
        VTU_DIR / "hex8_patch_2x2x2.vtu",
        "Eight-element Hex8 patch with an affine uniaxial strain field",
        points,
        cells,
        12,
        displacements,
        markers,
        {
            "PatchElementId": [float(i + 1) for i in range(len(cells))],
            "SXX": [1.0] * len(cells),
            "SYY": [0.0] * len(cells),
            "SZZ": [0.0] * len(cells),
            "SXY": [0.0] * len(cells),
            "SYZ": [0.0] * len(cells),
            "SXZ": [0.0] * len(cells),
            "VonMises": [1.0] * len(cells),
        },
    )
    for form in PATCH_FORMS:
        write_marker_vtu(
            VTU_DIR / f"hex8_patch_2x2x2_{form['key']}_markers.vtu",
            f"Hex8 patch markers for {form['title']}",
            points,
            displacements,
            markers[form["marker"]],
        )


def generate_bbarhex8_patch():
    nx = ny = nz = 2

    def node(i, j, k):
        return k * (ny + 1) * (nx + 1) + j * (nx + 1) + i

    points = [
        (i / nx, j / ny, k / nz)
        for k in range(nz + 1)
        for j in range(ny + 1)
        for i in range(nx + 1)
    ]
    cells = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                cells.append(
                    (
                        node(i, j, k),
                        node(i + 1, j, k),
                        node(i + 1, j + 1, k),
                        node(i, j + 1, k),
                        node(i, j, k + 1),
                        node(i + 1, j, k + 1),
                        node(i + 1, j + 1, k + 1),
                        node(i, j + 1, k + 1),
                    )
                )
    displacements = [
        (1.0e-3 * x, -3.0e-4 * y, -3.0e-4 * z) for x, y, z in points
    ]
    markers = patch_point_data(points, "xyz")
    write_vtu(
        VTU_DIR / "bbarhex8_patch_2x2x2.vtu",
        "Eight-element BbarHex8 patch with an affine uniaxial strain field",
        points,
        cells,
        12,
        displacements,
        markers,
        {
            "PatchElementId": [float(i + 1) for i in range(len(cells))],
            "SXX": [1.0] * len(cells),
            "SYY": [0.0] * len(cells),
            "SZZ": [0.0] * len(cells),
            "SXY": [0.0] * len(cells),
            "SYZ": [0.0] * len(cells),
            "SXZ": [0.0] * len(cells),
            "VonMises": [1.0] * len(cells),
        },
    )
    for form in PATCH_FORMS:
        write_marker_vtu(
            VTU_DIR / f"bbarhex8_patch_2x2x2_{form['key']}_markers.vtu",
            f"BbarHex8 patch markers for {form['title']}",
            points,
            displacements,
            markers[form["marker"]],
        )


def generate_subdivided_vtus():
    generate_bar_patch()
    generate_beam_patch()
    generate_plate_patch()
    generate_hex_patch()
    generate_bbarhex8_patch()


def set_camera(view, mode):
    view.ResetCamera(False, 0.9)
    if mode == "xy":
        view.CameraPosition = [0.5, 0.5, 3.0]
        view.CameraFocalPoint = [0.5, 0.5, 0.0]
        view.CameraViewUp = [0.0, 1.0, 0.0]
    else:
        view.CameraPosition = [2.8, -3.2, 2.4]
        view.CameraFocalPoint = [0.5, 0.5, 0.5]
        view.CameraViewUp = [0.0, 0.0, 1.0]
    view.ResetCamera(False, 0.9)


def add_title(view, title, note):
    title_source = Text()
    title_source.Text = f"{title}\n{note}"
    title_display = Show(title_source, view)
    title_display.WindowLocation = "Upper Center"
    title_display.FontSize = 16
    title_display.Color = [0.05, 0.05, 0.05]
    return title_source


def marker_path(model, form):
    marker_names = {
        "bar": "bar_patch_6_segments",
        "beam": "beam_patch_6_segments",
        "plate": "plate4_patch_3x3",
        "hex": "hex8_patch_2x2x2",
        "bbarhex8": "bbarhex8_patch_2x2x2",
    }
    return VTU_DIR / f"{marker_names[model['key']]}_{form['key']}_markers.vtu"


def render_case(model, form):
    view = CreateRenderView()
    view.ViewSize = [1600, 1100]
    view.Background = [1.0, 1.0, 1.0]

    reader = XMLUnstructuredGridReader(FileName=[str(model["file"])])
    reader.PointArrayStatus = [
        "Displacement",
        "DisplacementMagnitude",
        "UX",
        "UY",
        "UZ",
        "AllNodeDisplacement",
        "BoundaryDisplacement",
        "EquivalentForceNode",
    ]
    reader.CellArrayStatus = [
        "VonMises",
        "SXX",
        "SYY",
        "SZZ",
        "SXY",
        "SYZ",
        "SXZ",
        "MX",
        "MY",
        "MXY",
        "AxialForce",
        "PatchElementId",
    ]

    base_display = Show(reader, view)
    base_display.Representation = "Wireframe"
    base_display.AmbientColor = [0.25, 0.25, 0.25]
    base_display.DiffuseColor = [0.25, 0.25, 0.25]
    base_display.LineWidth = 2.0

    warp = WarpByVector(Input=reader)
    warp.Vectors = ["POINTS", "Displacement"]
    warp.ScaleFactor = model["scale"]
    SetActiveSource(warp)

    warped_display = Show(warp, view)
    warped_display.Representation = "Surface With Edges"
    warped_display.EdgeColor = [0.08, 0.08, 0.08]
    warped_display.LineWidth = model.get("line_width", 1.5)
    try:
        warped_display.RenderLinesAsTubes = 1
    except AttributeError:
        pass
    ColorBy(warped_display, ("POINTS", "UX"))

    lut = GetColorTransferFunction("UX")
    low, high = model["ux_range"]
    lut.RescaleTransferFunction(low, high)
    lut.RGBPoints = [
        low,
        0.231373,
        0.298039,
        0.752941,
        high,
        0.705882,
        0.015686,
        0.149020,
    ]
    try:
        lut.ApplyPreset("Viridis (matplotlib)", True)
    except RuntimeError:
        pass
    lut.RescaleTransferFunction(low, high)
    warped_display.SetScalarBarVisibility(view, True)

    marker_reader = XMLUnstructuredGridReader(FileName=[str(marker_path(model, form))])
    marker_reader.PointArrayStatus = ["Displacement", "UX", "Marker"]
    marker_warp = WarpByVector(Input=marker_reader)
    marker_warp.Vectors = ["POINTS", "Displacement"]
    marker_warp.ScaleFactor = model["scale"]
    marker_display = Show(marker_warp, view)
    marker_display.Representation = "Points"
    marker_display.PointSize = 18 if model["key"] != "bar" else 22
    marker_display.DiffuseColor = [0.95, 0.05, 0.02]
    marker_display.AmbientColor = [0.95, 0.05, 0.02]
    try:
        marker_display.RenderPointsAsSpheres = 1
    except AttributeError:
        pass

    title = f"{form['title']} - {model['label']}"
    note = f"{model['note']} Red markers show {form['note'][0].lower()}{form['note'][1:]}"
    title_source = add_title(view, title, note)
    set_camera(view, model["camera"])

    screenshot = OUT_DIR / f"{model['prefix']}_{form['index']}_{form['key']}_ux.png"
    SaveScreenshot(str(screenshot), view, ImageResolution=[1600, 1100])

    HideScalarBarIfNotNeeded(lut, view)
    Hide(marker_warp, view)
    Hide(marker_reader, view)
    Hide(warp, view)
    Hide(reader, view)
    Delete(title_source)
    Delete(marker_warp)
    Delete(marker_reader)
    Delete(warp)
    Delete(reader)
    Delete(view)
    return screenshot


def write_summary(rendered):
    lines = [
        "# Patch test ParaView rendering",
        "",
        "Patch test checks whether an element formulation can reproduce an affine displacement field,",
        "therefore a constant strain and constant stress state, independent of mesh subdivision.",
        "",
        "## Three patch-test boundary forms",
        "",
        *PATCH_TEST_TYPES,
        "",
        "## Rendered STAP++ cases",
        "",
    ]
    for image in rendered:
        lines.append(f"- `{image.relative_to(ROOT)}`")
    lines.extend(
        [
            "",
            "The screenshots use explicit subdivided patch meshes so the element partition is visible:",
            "six 1D bar elements, a 3x3 Plate4 membrane patch, and a 2x2x2 Hex8 solid patch.",
            "All main screenshots color the warped patch by UX to show the affine/linear displacement field.",
            "Red markers indicate the nodes controlled by each patch-test boundary form.",
            "",
            "The generated VTU files are visualization patch fields derived from the verified",
            "single-element STAP++ patch results in this repository.",
            "",
        ]
    )
    (OUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    ResetSession()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    generate_subdivided_vtus()
    rendered = []
    for model in MODELS:
        if not model["file"].exists():
            raise FileNotFoundError(model["file"])
        for form in PATCH_FORMS:
            marker = marker_path(model, form)
            if not marker.exists():
                raise FileNotFoundError(marker)
            rendered.append(render_case(model, form))
    write_summary(rendered)
    try:
        layout = GetLayout()
    except RuntimeError:
        layout = None
    if layout:
        Delete(layout)


if __name__ == "__main__":
    main()
