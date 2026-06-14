"""
桥梁模型 Abaqus 运行与对比脚本
===============================
用法:
  1. 在有 Abaqus 的机器上运行 Bridge INP:
     abaqus job=Bridge-2 input=Bridge-2.inp cpus=8 double
     abaqus job=Bridge-1 input=Bridge-1.inp cpus=4 double

  2. 从 ODB 提取节点位移:
     abaqus python abaqus_scripts/extract_odb_results.py --odb Bridge-2.odb
     abaqus python abaqus_scripts/extract_odb_results.py --odb Bridge-1.odb

  3. 运行本脚本进行对比:
     python tools/compare_bridges.py
"""

import csv
import math
import re
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DISP_ROW = re.compile(
    r"^\s*(\d+)\s+([+-]?(?:\d+\.\d*|\d*\.\d+|\d+)[Ee][+-]?\d+)\s+"
    r"([+-]?(?:\d+\.\d*|\d*\.\d+|\d+)[Ee][+-]?\d+)\s+"
    r"([+-]?(?:\d+\.\d*|\d*\.\d+|\d+)[Ee][+-]?\d+)"
)


# ===================== 配置 =====================
# 每个桥梁的 STAP++ 输出和 Abaqus CSV 位置
BRIDGE_CONFIGS = {
    "Bridge-1": {
        "stappp_out": ROOT / "data" / "data-3" / "Bridge-1.rcm.generated.out",
        "abaqus_nodes": ROOT / "data" / "data-3" / "Bridge-1_nodes.csv",
        "abaqus_inp": ROOT / "data" / "data-3" / "Bridge-1.inp",
    },
    "Bridge-2": {
        "stappp_out": ROOT / "data" / "data-3" / "Bridge-2.eigen_mpc.out",
        "abaqus_nodes": ROOT / "data" / "data-3" / "Bridge-2_nodes.csv",
        "abaqus_inp": ROOT / "data" / "data-3" / "Bridge-2.inp",
    },
    "Bridge-3": {
        "stappp_out": ROOT / "cases" / "bridge3_run" / "stappp" / "bridge3_run.out",
        "abaqus_nodes": ROOT / "cases" / "bridge3_run" / "abaqus" / "Bridge-3_nodes.csv",
        "abaqus_inp": ROOT / "data" / "data-3" / "Bridge-3.inp",
    },
}

# 桥梁描述信息
BRIDGE_INFO = {
    "Bridge-1": {
        "name": "Bridge-1 (人行桥 RCM)",
        "description": "小型人行桥，采用 RCM 节点重排序优化带宽，Skyline + LDLᵀ 求解",
        "solver": "Skyline LDLᵀ (RCM 优化)",
    },
    "Bridge-2": {
        "name": "Bridge-2 (MPC/Eigen)",
        "description": "中型桥梁，含多点约束 (MPC)，Eigen 稀疏直接求解器",
        "solver": "Eigen Sparse LDLT (MPC 消去法)",
    },
    "Bridge-3": {
        "name": "Bridge-3 (大跨悬索桥)",
        "description": "大型悬索桥，含索、塔、梁单元，MPC 耦合 + Tie 约束",
        "solver": "Eigen Sparse LDLT (refrot + norot 修正)",
    },
}


def read_stappp_displacements(path: Path) -> dict[int, tuple[float, float, float]]:
    """从 STAP++ .out 文件读取节点位移。"""
    if not path.exists():
        return {}
    values = {}
    in_displacements = False
    text = path.read_text(errors="ignore")
    for line in text.splitlines():
        if "D I S P L A C E M E N T S" in line:
            in_displacements = True
            continue
        if in_displacements and ("S T R E S S" in line or "R E A C T I O N" in line):
            break
        if not in_displacements:
            continue
        match = DISP_ROW.match(line)
        if match:
            values[int(match.group(1))] = tuple(float(match.group(i)) for i in range(2, 5))
    return values


def read_abaqus_nodes(path: Path) -> dict[str, tuple[float, float, float]]:
    """从 Abaqus _nodes.csv 读取节点位移 (按 instance:node 索引)。"""
    if not path.exists():
        return {}
    values = {}
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            key = f"{row['instance']}:{row['node']}"
            values[key] = (float(row["u1"]), float(row["u2"]), float(row["u3"]))
    return values


def read_stappp_nodes_by_coord(path: Path) -> dict[tuple[float, float, float], int]:
    """从 STAP++ .out 读取节点坐标，返回 {坐标: 节点号}。"""
    if not path.exists():
        return {}
    coords = {}
    in_nodes = False
    text = path.read_text(errors="ignore")
    for line in text.splitlines():
        if "N O D A L   P O I N T   D A T A" in line:
            in_nodes = True
            continue
        if in_nodes and ("L O A D" in line or "E L E M E N T" in line):
            break
        if not in_nodes:
            continue
        parts = line.split()
        if len(parts) >= 10:
            try:
                node = int(parts[0])
                x, y, z = float(parts[7]), float(parts[8]), float(parts[9])
                coords[(x, y, z)] = node
            except (ValueError, IndexError):
                continue
    return coords


def read_abaqus_inp_nodes(path: Path) -> dict[str, tuple[float, float, float]]:
    """从 Abaqus .inp 读取所有 Part 的节点坐标，返回 {"PART-NAME:NODE-ID": (x,y,z)}。"""
    if not path.exists():
        return {}
    
    text = path.read_text(errors="ignore")
    lines = text.splitlines()
    
    nodes: dict[str, tuple[float, float, float]] = {}
    current_part: str | None = None
    in_node_section = False
    
    for line in lines:
        stripped = line.strip()
        upper = stripped.upper()
        
        # 追踪 Part name
        if upper.startswith("*PART"):
            # *Part, name=Part-Cable50  → 提取 name=
            current_part = None
            for token in stripped.split(","):
                token = token.strip()
                if "=" in token:
                    key, val = token.split("=", 1)
                    if key.strip().upper() == "NAME":
                        current_part = val.strip()
                        break
            if current_part is None:
                # 没有 name=，用 Part 后的第一个 token
                name = stripped.replace("*Part", "").replace("*PART", "").strip().lstrip(",").strip()
                if name:
                    current_part = name
            in_node_section = False
            continue
        
        if upper.startswith("*END PART") or upper.startswith("*ENDPART"):
            current_part = None
            in_node_section = False
            continue
        
        if upper.startswith("*ASSEMBLY"):
            current_part = "ASSEMBLY"
            in_node_section = False
            continue
        
        if upper.startswith("*END ASSEMBLY"):
            current_part = None
            in_node_section = False
            continue
        
        # 检测 *Node 块的开始
        if upper.startswith("*NODE"):
            if current_part:
                in_node_section = True
            continue
        
        # 遇到下一个 * 关键字 → 结束当前节点块
        if in_node_section and stripped.startswith("*"):
            in_node_section = False
            continue
        
        # 读取节点行
        if in_node_section and current_part:
            parts = stripped.split(",")
            if len(parts) >= 4:
                try:
                    local_id = int(parts[0].strip())
                    x = float(parts[1].strip())
                    y = float(parts[2].strip())
                    z = float(parts[3].strip())
                    key = f"{current_part}:{local_id}"
                    nodes[key] = (x, y, z)
                except (ValueError, IndexError):
                    continue
    
    return nodes


def match_by_coordinate(
    stap_coords: dict[tuple[float, float, float], int],
    abaqus_nodes: dict[str, tuple[float, float, float]],
) -> dict[int, str]:
    """通过坐标匹配 STAP++ 全局节点号与 Abaqus instance:node。"""
    # 构建 Abaqus 坐标 → instance:node 映射（四舍五入到 3 位小数）
    abaqus_coord_map: dict[tuple[float, float, float], str] = {}
    for key, coord in abaqus_nodes.items():
        rounded = (round(coord[0], 3), round(coord[1], 3), round(coord[2], 3))
        abaqus_coord_map[rounded] = key

    matched: dict[int, str] = {}
    for coord, stap_node in stap_coords.items():
        rounded = (round(coord[0], 3), round(coord[1], 3), round(coord[2], 3))
        if rounded in abaqus_coord_map:
            matched[stap_node] = abaqus_coord_map[rounded]

    return matched


def compute_rms(errors: list[float]) -> float:
    if not errors:
        return 0.0
    return math.sqrt(sum(e * e for e in errors) / len(errors))


def compare_bridge(name: str, config: dict) -> dict | None:
    """对比单个桥梁的 STAP++ 和 Abaqus 位移结果。"""
    print(f"\n{'='*60}")
    print(f"  对比: {BRIDGE_INFO[name]['name']}")
    print(f"{'='*60}")

    stap_path = config["stappp_out"]
    abaqus_csv = config["abaqus_nodes"]
    abaqus_inp = config["abaqus_inp"]

    if not stap_path.exists():
        print(f"  ⚠ STAP++ 输出缺失: {stap_path}")
        return None
    if not abaqus_csv.exists():
        print(f"  ⚠ Abaqus 节点 CSV 缺失: {abaqus_csv}")
        print(f"  → 请先运行 Abaqus: abaqus job={name} input={name}.inp")
        return None
    if not abaqus_inp.exists():
        print(f"  ⚠ Abaqus INP 缺失: {abaqus_inp}")
        return None

    # 读 STAP++ 位移
    stap_disp = read_stappp_displacements(stap_path)
    print(f"  STAP++ 节点数: {len(stap_disp)}")

    # 读 Abaqus INP 节点坐标 (per-instance)
    abaqus_nodes_coord = read_abaqus_inp_nodes(abaqus_inp)
    print(f"  Abaqus INP 节点数: {len(abaqus_nodes_coord)}")

    # 读 Abaqus 位移
    abaqus_disp = read_abaqus_nodes(abaqus_csv)
    print(f"  Abaqus 位移行数: {len(abaqus_disp)}")

    # 读 STAP++ 节点坐标并匹配
    stap_coords = read_stappp_nodes_by_coord(stap_path)
    print(f"  STAP++ 坐标映射: {len(stap_coords)} 节点")

    matched = match_by_coordinate(stap_coords, abaqus_nodes_coord)
    print(f"  坐标匹配节点数: {len(matched)}")

    if len(matched) == 0:
        # 回退: 按节点号直接匹配 (如果 STAP++ 和 Abaqus 节点编号一致)
        print("  坐标匹配失败，尝试节点号直接匹配...")
        # Abaqus CSV 使用 instance:node 格式，无法直接匹配
        print("  ✗ 无法匹配，需要 Abaqus 结果文件")
        return None

    # 逐节点对比
    errors_u1 = []
    errors_u2 = []
    errors_u3 = []
    errors_vec = []
    max_err = {"node": 0, "u1": 0, "u2": 0, "u3": 0, "vec": 0}

    for stap_node, abaqus_key in sorted(matched.items()):
        u_stap = stap_disp.get(stap_node)
        u_aba = abaqus_disp.get(abaqus_key)
        if u_stap is None or u_aba is None:
            continue

        du1 = abs(u_stap[0] - u_aba[0])
        du2 = abs(u_stap[1] - u_aba[1])
        du3 = abs(u_stap[2] - u_aba[2])
        vec_err = math.sqrt(du1**2 + du2**2 + du3**2)

        errors_u1.append(du1)
        errors_u2.append(du2)
        errors_u3.append(du3)
        errors_vec.append(vec_err)

        if vec_err > max_err["vec"]:
            max_err = {
                "node": stap_node,
                "abaqus_key": abaqus_key,
                "u1": du1,
                "u2": du2,
                "u3": du3,
                "vec": vec_err,
            }

    total = len(errors_vec)
    if total == 0:
        print("  ✗ 无有效对比数据")
        return None

    result = {
        "name": name,
        "display_name": BRIDGE_INFO[name]["name"],
        "description": BRIDGE_INFO[name]["description"],
        "solver": BRIDGE_INFO[name]["solver"],
        "total_nodes": total,
        "stap_nodes": len(stap_disp),
        "abaqus_nodes": len(abaqus_disp),
        "rms_u1": compute_rms(errors_u1),
        "rms_u2": compute_rms(errors_u2),
        "rms_u3": compute_rms(errors_u3),
        "rms_vec": compute_rms(errors_vec),
        "max_abs_u1": max(errors_u1) if errors_u1 else 0,
        "max_abs_u2": max(errors_u2) if errors_u2 else 0,
        "max_abs_u3": max(errors_u3) if errors_u3 else 0,
        "max_abs_vec": max(errors_vec) if errors_vec else 0,
        "max_node": max_err["node"],
    }

    # 输出结果
    print(f"\n  对比节点数: {total}")
    print(f"  位移向量 RMSE: {result['rms_vec']:.6e}")
    print(f"  分量 RMSE: U1={result['rms_u1']:.6e}, U2={result['rms_u2']:.6e}, U3={result['rms_u3']:.6e}")
    print(f"  最大绝对误差:")
    print(f"    U1={result['max_abs_u1']:.6e}, U2={result['max_abs_u2']:.6e}, U3={result['max_abs_u3']:.6e}")
    print(f"    向量={result['max_abs_vec']:.6e} @ 节点 {max_err['node']}")

    return result


def main():
    os.makedirs(ROOT / "cases", exist_ok=True)

    print("=" * 60)
    print("  STAP++ vs Abaqus 桥梁模型位移对比")
    print("=" * 60)

    results = []
    for name in ["Bridge-1", "Bridge-2", "Bridge-3"]:
        result = compare_bridge(name, BRIDGE_CONFIGS[name])
        if result:
            results.append(result)

    if not results:
        print("\n无有效对比结果。")
        print("\n请先运行 Abaqus 作业:")
        for name in ["Bridge-1", "Bridge-2"]:
            print(f"  abaqus job={name} input={name}.inp cpus=4 double")
            print(f"  abaqus python abaqus_scripts/extract_odb_results.py --odb {name}.odb")
        return

    # 写入汇总 CSV
    summary_path = ROOT / "cases" / "bridge_comparison.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "bridge", "nodes", "rms_vec", "rms_u1", "rms_u2", "rms_u3",
            "max_u1", "max_u2", "max_u3", "max_vec", "max_node"
        ])
        writer.writeheader()
        for r in results:
            writer.writerow({
                "bridge": r["name"], "nodes": r["total_nodes"],
                "rms_vec": r["rms_vec"], "rms_u1": r["rms_u1"],
                "rms_u2": r["rms_u2"], "rms_u3": r["rms_u3"],
                "max_u1": r["max_abs_u1"], "max_u2": r["max_abs_u2"],
                "max_u3": r["max_abs_u3"], "max_vec": r["max_abs_vec"],
                "max_node": r["max_node"],
            })
    print(f"\n汇总已保存: {summary_path}")

    # 输出汇总表
    print(f"\n{'='*80}")
    print("  桥梁对比汇总")
    print(f"{'='*80}")
    header = f"{'模型':<20} {'对比节点':>8} {'Vec RMSE':>12} {'Max Vec':>12}"
    print(header)
    print("-" * 80)
    for r in results:
        print(f"{r['display_name']:<20} {r['total_nodes']:>8} {r['rms_vec']:>12.6e} {r['max_abs_vec']:>12.6e}")

    print(f"\n{'='*80}")
    print("  详细分量对比")
    print(f"{'='*80}")
    detail_header = f"{'模型':<20} {'U1 RMS':>12} {'U2 RMS':>12} {'U3 RMS':>12} {'U1 Max':>12} {'U2 Max':>12} {'U3 Max':>12}"
    print(detail_header)
    print("-" * 80)
    for r in results:
        print(f"{r['display_name']:<20} {r['rms_u1']:>12.6e} {r['rms_u2']:>12.6e} {r['rms_u3']:>12.6e} "
              f"{r['max_abs_u1']:>12.6e} {r['max_abs_u2']:>12.6e} {r['max_abs_u3']:>12.6e}")


if __name__ == "__main__":
    main()
