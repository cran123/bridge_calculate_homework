"""
批量将 cases/ 下的 .vtu 文件导出为 PNG 截图。
用法: python tools/export_vtu_to_png.py

依赖: pip install vtk
若 VTK 不可用, 可改用 pyvista: pip install pyvista
"""

import os
import sys
import glob

# ---- 配置 ----
CASES_DIR = os.path.join(os.path.dirname(__file__), "..", "cases")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "doc", "figures")
IMAGE_WIDTH = 1200
IMAGE_HEIGHT = 900


def find_vtu_files(cases_dir: str) -> list[tuple[str, str]]:
    """扫描 cases/ 下所有 paraview/*.vtu 文件, 返回 [(case_name, vtu_path), ...]."""
    results = []
    for case_name in sorted(os.listdir(cases_dir)):
        case_path = os.path.join(cases_dir, case_name)
        if not os.path.isdir(case_path):
            continue
        pv_dir = os.path.join(case_path, "paraview")
        if not os.path.isdir(pv_dir):
            continue
        for f in sorted(os.listdir(pv_dir)):
            if f.endswith(".vtu"):
                results.append((case_name, os.path.join(pv_dir, f)))
    return results


# ============== VTK 方案 ==============
def export_with_vtk(vtu_files: list[tuple[str, str]], out_dir: str):
    """使用 vtk 模块进行离屏渲染并保存 PNG。"""
    import vtk
    from vtkmodules.vtkRenderingCore import (
        vtkRenderWindow, vtkRenderWindowInteractor, vtkRenderer,
        vtkActor, vtkPolyDataMapper, vtkDataSetMapper,
    )
    from vtkmodules.vtkRenderingAnnotation import vtkScalarBarActor
    from vtkmodules.vtkIOImage import vtkPNGWriter
    from vtkmodules.vtkRenderingOpenGL2 import vtkOSOpenGLRenderWindow

    # 检查 Displacement 向量是否存在
    def has_vector(reader, name):
        pd = reader.GetOutput()
        if not pd:
            return False
        point_data = pd.GetPointData()
        for i in range(point_data.GetNumberOfArrays()):
            if point_data.GetArrayName(i) == name:
                return True
        return False

    # 检查标量是否存在
    def has_scalar(reader, name):
        pd = reader.GetOutput()
        if not pd:
            return False
        point_data = pd.GetPointData()
        for i in range(point_data.GetNumberOfArrays()):
            if point_data.GetArrayName(i) == name:
                return True
        return False

    for case_name, vtu_path in vtu_files:
        print(f"[VTK] 处理: {case_name} -> {vtu_path}")

        # ---- 读取 ----
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtu_path)
        reader.Update()

        # ---- Warp By Vector (位移变形) ----
        warp = vtk.vtkWarpVector()
        warp.SetInputConnection(reader.GetOutputPort())
        warp.SetScaleFactor(1.0)
        if has_vector(reader, "Displacement"):
            pd = reader.GetOutput()
            pd.GetPointData().SetActiveVectors("Displacement")
            warp.SetInputArrayToProcess(
                0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "Displacement"
            )
        else:
            # 无位移向量, 跳过变形
            warp.SetScaleFactor(0.0)

        warp.Update()

        # ---- Mapper ----
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputConnection(warp.GetOutputPort())
        mapper.ScalarVisibilityOn()

        # 着色: 优先 DisplacementMagnitude, 否则尝试其他
        if has_scalar(reader, "DisplacementMagnitude"):
            warp.GetOutput().GetPointData().SetActiveScalars("DisplacementMagnitude")
            mapper.SetScalarModeToUsePointData()
        elif has_scalar(reader, "VonMises"):
            warp.GetOutput().GetPointData().SetActiveScalars("VonMises")
            mapper.SetScalarModeToUsePointData()

        mapper.SetScalarRange(warp.GetOutput().GetPointData().GetScalars().GetRange())

        # ---- Actor ----
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        # ---- Scalar Bar ----
        scalar_bar = vtkScalarBarActor()
        scalar_bar.SetLookupTable(mapper.GetLookupTable())
        scalar_bar.SetTitle("Displacement Mag")
        scalar_bar.SetNumberOfLabels(5)

        # ---- Renderer ----
        renderer = vtkRenderer()
        renderer.AddActor(actor)
        renderer.AddActor2D(scalar_bar)
        renderer.SetBackground(1.0, 1.0, 1.0)  # 白底
        renderer.ResetCamera()

        # ---- Render Window (headless) ----
        ren_win = vtk.vtkRenderWindow()
        ren_win.SetOffScreenRendering(1)
        ren_win.AddRenderer(renderer)
        ren_win.SetSize(IMAGE_WIDTH, IMAGE_HEIGHT)
        ren_win.Render()

        # ---- 导出 PNG ----
        w2i = vtk.vtkWindowToImageFilter()
        w2i.SetInput(ren_win)
        w2i.Update()

        writer = vtkPNGWriter()
        out_path = os.path.join(out_dir, f"{case_name}.png")
        writer.SetFileName(out_path)
        writer.SetInputConnection(w2i.GetOutputPort())
        writer.Write()

        print(f"  -> 已保存: {out_path}")


# ============== PyVista 方案 ==============
def export_with_pyvista(vtu_files: list[tuple[str, str]], out_dir: str):
    """使用 pyvista 进行离屏渲染 (更简洁的备选方案)。"""
    import pyvista as pv

    for case_name, vtu_path in vtu_files:
        print(f"[PyVista] 处理: {case_name} -> {vtu_path}")

        mesh = pv.read(vtu_path)

        # 设置标量
        scalars = None
        if "DisplacementMagnitude" in mesh.point_data:
            scalars = "DisplacementMagnitude"
        elif "VonMises" in mesh.point_data:
            scalars = "VonMises"

        # Warp by vector
        warped = mesh
        if "Displacement" in mesh.point_data:
            warped = mesh.warp_by_vector("Displacement", factor=1.0)

        # Plot
        plotter = pv.Plotter(off_screen=True, window_size=(IMAGE_WIDTH, IMAGE_HEIGHT))
        plotter.background_color = "white"
        plotter.add_mesh(warped, scalars=scalars, cmap="jet",
                         show_scalar_bar=True, scalar_bar_args={"title": scalars or ""})
        plotter.view_isometric()
        plotter.camera.zoom(1.2)

        out_path = os.path.join(out_dir, f"{case_name}.png")
        plotter.screenshot(out_path)
        plotter.close()
        print(f"  -> 已保存: {out_path}")


# ============== main ==============
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    vtu_files = find_vtu_files(CASES_DIR)
    if not vtu_files:
        print("未找到任何 .vtu 文件!")
        return

    print(f"找到 {len(vtu_files)} 个 .vtu 文件:")
    for name, path in vtu_files:
        print(f"  [{name}] {path}")

    # 尝试 VTK, 失败则尝试 PyVista
    try:
        import vtk
        print("\n使用 VTK 引擎导出...")
        export_with_vtk(vtu_files, OUTPUT_DIR)
    except ImportError:
        try:
            import pyvista
            print("\nVTK 未安装, 使用 PyVista 引擎导出...")
            export_with_pyvista(vtu_files, OUTPUT_DIR)
        except ImportError:
            print("\n错误: 未安装 VTK 或 PyVista!")
            print("请执行以下命令之一来安装:")
            print("  conda install -c conda-forge vtk")
            print("  pip install vtk")
            print("  pip install pyvista")

    print(f"\n全部完成! 图片保存在: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
