# 团队单元合并接口约定

本文档用于 H8、Beam、Truss、Plate 四类单元合并前对齐接口。各分支合并到 `dev` 前必须按本文档检查，避免单元编号、输入格式、自重和后处理字段互相冲突。

## 1. 单元编号

程序内部统一使用 `ElementTypes` 枚举，输入文件中单元组第一项即为单元类型编号。

| 编号 | ABAQUS 类型 | 程序类名 | 说明 |
| --- | --- | --- | --- |
| 1 | T3D2 | `CBar` / `CTruss3D2` | 三维二节点杆单元，现有 `Bar` 可作为基础 |
| 4 | C3D8 / C3D8R | `CHex8` | 八节点实体单元 |
| 5 | B31 | `CBeam3D2` | 三维二节点梁单元 |
| 6 | S4R | `CPlate4` | 工程版四节点板壳单元 |

合并要求：

- 不允许重复占用编号。
- `ElementGroup::CalculateMemberSize()`、`AllocateElements()`、`AllocateMaterials()` 中必须同时注册单元和材料。
- `Outputter::OutputElementInfo()` 和 `Outputter::OutputElementStress()` 必须有对应分支。
- 遇到未实现类型时保留清晰报错。

## 2. 材料格式

所有材料行第一项统一为 `set`，并要求材料集按 `1, 2, 3...` 顺序输入。

| 单元 | 材料行格式 |
| --- | --- |
| Truss3D2 | `set E A rho` |
| H8 | `set E nu rho` |
| Beam3D2 | `set E nu A Iy Iz J rho` |
| Plate4 | `set E nu thickness rho` |

合并要求：

- `rho` 统一表示质量密度，单位按输入模型一致处理。
- `nu` 统一表示泊松比。
- Beam 的剪切模量统一由 `G = E / (2 * (1 + nu))` 得到。
- 如果 Beam 后续需要截面外缘距离用于最大应力，可新增 `cy cz`，但必须同步更新文档和转换器。

## 3. 自重接口

所有单元统一通过基类接口提供自重等效节点力：

```cpp
virtual void ElementBodyForce(double* bodyForce, const double gravity[3]);
```

约定：

- `bodyForce` 长度必须等于该单元 `ND_`。
- 函数内部先清零 `bodyForce`。
- `gravity[3]` 是加速度向量，例如 `0 0 -9.81`。
- 返回值按单元局部位置矩阵顺序排列。
- 不支持自重的单元使用基类默认零向量。

各单元第一版自重策略：

| 单元 | 自重等效荷载 |
| --- | --- |
| Truss3D2 | `rho * A * L * gravity / 2` 分配到两个节点平动自由度 |
| H8 | 使用 Gauss 积分计算 `rho * N^T * gravity dV` |
| Beam3D2 | `rho * A * L * gravity / 2` 分配到两个节点平动自由度 |
| Plate4 | `rho * thickness * area * gravity / 4` 分配到四个节点平动自由度 |

荷载工况输入格式：

```text
loadcase_id
nloads gravity_flag gx gy gz
node dof value
...
```

兼容旧格式：

```text
loadcase_id
nloads
node dof value
...
```

## 4. VTU 字段

ParaView 后处理统一输出 `.vtu`，每个荷载工况一个文件：

```text
model_lc1.vtu
model_lc2.vtu
```

VTK 单元类型：

| 程序类型 | VTK 类型 |
| --- | --- |
| Truss3D2 | `VTK_LINE = 3` |
| Beam3D2 | `VTK_LINE = 3` |
| Plate4 | `VTK_QUAD = 9` |
| H8 | `VTK_HEXAHEDRON = 12` |

PointData 统一字段：

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `Displacement` | 3 分量 | `ux uy uz` |
| `Rotation` | 3 分量 | `rx ry rz`，实体/杆可为 0 |
| `DisplacementMagnitude` | 标量 | 位移模长 |

CellData 统一字段：

| 字段 | 说明 |
| --- | --- |
| `ElementType` | 程序单元编号 |
| `VonMises` | 实体和板壳优先输出，其他单元可为 0 |
| `SXX SYY SZZ SXY SYZ SXZ` | 实体应力；非实体缺省为 0 |
| `SXX SYY SXY` | 板壳膜应力；非板壳缺省为 0 |
| `MX MY MXY` | 板壳弯矩；非板壳缺省为 0 |
| `AxialForce` | 杆/梁轴力；其他单元缺省为 0 |
| `BeamMy BeamMz BeamTorque` | 梁内力；其他单元缺省为 0 |

合并要求：

- 字段名不要随意改，报告和 ParaView 截图按这些字段引用。
- 某单元没有该字段物理意义时输出 `0`，不要省略字段。
- 后续若新增字段，应保持旧字段不变。

## 5. 混合模型 DAT 输入格式

统一输入结构：

```text
title
NUMNP NUMEG NLCASE MODEX
node bcx bcy bcz bcrx bcry bcrz x y z
...
loadcase_id
nloads gravity_flag gx gy gz
node dof value
...
element_type NUME NUMMAT
material lines...
element lines...
...
```

节点自由度顺序：

```text
ux uy uz rx ry rz
```

边界码：

| 值 | 含义 |
| --- | --- |
| `0` | 自由，参与求解 |
| `1` | 约束，不参与求解 |

单元组顺序建议：

```text
4 H8
6 Plate4
5 Beam3D2
1 Truss3D2
```

各单元连接格式：

```text
H8:
ele n1 n2 n3 n4 n5 n6 n7 n8 mat

Plate4:
ele n1 n2 n3 n4 mat

Beam3D2:
ele n1 n2 mat vx vy vz

Truss3D2:
ele n1 n2 mat
```

合并要求：

- 全局节点统一 6 自由度。
- 实体和杆只取平动自由度。
- 梁和板壳使用平动与转角自由度。
- 各单元必须自己覆盖 `GenerateLocationMatrix()`，只取本单元实际使用的自由度。

## 6. Bridge-1 转换与整体求解

`Bridge-1.inp` 到 `.dat` 的转换器统一放在：

```text
tools/abaqus_to_stappp.py
```

第一版必须识别：

```text
*Node
*Element, type=C3D8R
*Element, type=S4R
*Element, type=B31
*Element, type=T3D2
*Material
*Solid Section
*Shell Section
*Beam Section
*Boundary
*Cload
```

转换策略：

- ABAQUS `C3D8R` 输出为程序 `H8`，编号 `4`。
- ABAQUS `S4R` 输出为程序 `Plate4`，编号 `6`。
- ABAQUS `B31` 输出为程序 `Beam3D2`，编号 `5`。
- ABAQUS `T3D2` 输出为程序 `Truss3D2`，编号 `1`。
- 如果 inp 中没有显式重力，转换器通过命令行参数写入 `gravity_flag gx gy gz`。

转换后必须打印统计：

```text
Nodes: ...
C3D8R/H8: ...
S4R/Plate4: ...
B31/Beam3D2: ...
T3D2/Truss3D2: ...
Materials: ...
Boundary DOFs: ...
Load cases: ...
```

整体求解验收：

- `Bridge-1.generated.dat` 能完成数据检查。
- 至少一个自重荷载工况能进入求解。
- 输出 `.out` 包含节点位移和各类单元结果。
- 输出 `.vtu` 能在 ParaView 打开。
- 报告中对比 ABAQUS 的最大位移、关键节点位移和关键单元应力/内力。

## 7. 合并前检查清单

每个单元分支合并前至少确认：

- 单元编号和材料格式符合本文档。
- `GenerateLocationMatrix()` 只包含本单元实际自由度。
- `ElementStiffness()` 输出 STAPpp 所需的上三角打包格式。
- `ElementBodyForce()` 已实现或明确返回零。
- `ElementStress()` 输出不导致 `Outputter` 崩溃。
- 至少有一个单元测试 `.dat`。
- 至少跑通一次 `cmake --build build`。
- 不提交 `*.out`、`*.vtu`、`build/`。
