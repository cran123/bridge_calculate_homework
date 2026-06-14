#!/usr/bin/env python3

import argparse
import math


def build_d_matrix(E, nu):
    factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
    a = (1.0 - nu) * factor
    b = nu * factor
    c = 0.5 * (1.0 - 2.0 * nu) * factor
    D = [[0.0 for _ in range(6)] for _ in range(6)]
    D[0][0] = a; D[0][1] = b; D[0][2] = b
    D[1][0] = b; D[1][1] = a; D[1][2] = b
    D[2][0] = b; D[2][1] = b; D[2][2] = a
    D[3][3] = c
    D[4][4] = c
    D[5][5] = c
    return D


def compute_shape(r, s, t):
    rp = 1.0 + r
    rm = 1.0 - r
    sp = 1.0 + s
    sm = 1.0 - s
    tp = 1.0 + t
    tm = 1.0 - t

    N = [0.0] * 8
    dNdr = [0.0] * 8
    dNds = [0.0] * 8
    dNdt = [0.0] * 8

    N[0] = 0.125 * rm * sm * tm
    N[1] = 0.125 * rp * sm * tm
    N[2] = 0.125 * rp * sp * tm
    N[3] = 0.125 * rm * sp * tm
    N[4] = 0.125 * rm * sm * tp
    N[5] = 0.125 * rp * sm * tp
    N[6] = 0.125 * rp * sp * tp
    N[7] = 0.125 * rm * sp * tp

    dNdr[0] = -0.125 * sm * tm
    dNdr[1] =  0.125 * sm * tm
    dNdr[2] =  0.125 * sp * tm
    dNdr[3] = -0.125 * sp * tm
    dNdr[4] = -0.125 * sm * tp
    dNdr[5] =  0.125 * sm * tp
    dNdr[6] =  0.125 * sp * tp
    dNdr[7] = -0.125 * sp * tp

    dNds[0] = -0.125 * rm * tm
    dNds[1] = -0.125 * rp * tm
    dNds[2] =  0.125 * rp * tm
    dNds[3] =  0.125 * rm * tm
    dNds[4] = -0.125 * rm * tp
    dNds[5] = -0.125 * rp * tp
    dNds[6] =  0.125 * rp * tp
    dNds[7] =  0.125 * rm * tp

    dNdt[0] = -0.125 * rm * sm
    dNdt[1] = -0.125 * rp * sm
    dNdt[2] = -0.125 * rp * sp
    dNdt[3] = -0.125 * rm * sp
    dNdt[4] =  0.125 * rm * sm
    dNdt[5] =  0.125 * rp * sm
    dNdt[6] =  0.125 * rp * sp
    dNdt[7] =  0.125 * rm * sp

    return N, dNdr, dNds, dNdt


def compute_jacobian(nodes, dNdr, dNds, dNdt):
    J = [[0.0 for _ in range(3)] for _ in range(3)]
    for i in range(8):
        J[0][0] += dNdr[i] * nodes[i][0]
        J[0][1] += dNds[i] * nodes[i][0]
        J[0][2] += dNdt[i] * nodes[i][0]

        J[1][0] += dNdr[i] * nodes[i][1]
        J[1][1] += dNds[i] * nodes[i][1]
        J[1][2] += dNdt[i] * nodes[i][1]

        J[2][0] += dNdr[i] * nodes[i][2]
        J[2][1] += dNds[i] * nodes[i][2]
        J[2][2] += dNdt[i] * nodes[i][2]

    detJ = (J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
            - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
            + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]))

    if abs(detJ) <= 1e-12:
        return None, None

    invDet = 1.0 / detJ
    invJ = [[0.0 for _ in range(3)] for _ in range(3)]
    invJ[0][0] =  (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * invDet
    invJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * invDet
    invJ[0][2] =  (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * invDet

    invJ[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]) * invDet
    invJ[1][1] =  (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * invDet
    invJ[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]) * invDet

    invJ[2][0] =  (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * invDet
    invJ[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]) * invDet
    invJ[2][2] =  (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * invDet

    return invJ, detJ


def build_b_matrix(dNdx, dNdy, dNdz):
    B = [[0.0 for _ in range(24)] for _ in range(6)]
    for i in range(8):
        col = i * 3
        B[0][col] = dNdx[i]
        B[1][col + 1] = dNdy[i]
        B[2][col + 2] = dNdz[i]

        B[3][col] = dNdy[i]
        B[3][col + 1] = dNdx[i]

        B[4][col + 1] = dNdz[i]
        B[4][col + 2] = dNdy[i]

        B[5][col] = dNdz[i]
        B[5][col + 2] = dNdx[i]
    return B


def mat_mul(A, B):
    m = len(A)
    n = len(B[0])
    k = len(B)
    out = [[0.0 for _ in range(n)] for _ in range(m)]
    for i in range(m):
        for j in range(n):
            s = 0.0
            for t in range(k):
                s += A[i][t] * B[t][j]
            out[i][j] = s
    return out


def mat_vec(A, x):
    m = len(A)
    n = len(A[0])
    out = [0.0 for _ in range(m)]
    for i in range(m):
        s = 0.0
        for j in range(n):
            s += A[i][j] * x[j]
        out[i] = s
    return out


def build_k_matrix(nodes, E, nu):
    D = build_d_matrix(E, nu)
    K = [[0.0 for _ in range(24)] for _ in range(24)]

    gp = [-1.0 / math.sqrt(3.0), 1.0 / math.sqrt(3.0)]
    for r in gp:
        for s in gp:
            for t in gp:
                N, dNdr, dNds, dNdt = compute_shape(r, s, t)
                invJ, detJ = compute_jacobian(nodes, dNdr, dNds, dNdt)
                if invJ is None:
                    continue
                dNdx = [invJ[0][0] * dNdr[i] + invJ[0][1] * dNds[i] + invJ[0][2] * dNdt[i] for i in range(8)]
                dNdy = [invJ[1][0] * dNdr[i] + invJ[1][1] * dNds[i] + invJ[1][2] * dNdt[i] for i in range(8)]
                dNdz = [invJ[2][0] * dNdr[i] + invJ[2][1] * dNds[i] + invJ[2][2] * dNdt[i] for i in range(8)]

                B = build_b_matrix(dNdx, dNdy, dNdz)
                DB = mat_mul(D, B)
                for i in range(24):
                    for j in range(24):
                        s = 0.0
                        for k in range(6):
                            s += B[k][i] * DB[k][j]
                        K[i][j] += s * detJ
    return K


def main():
    parser = argparse.ArgumentParser(description="Generate consistent nodal forces for Hex8 using f = K u")
    parser.add_argument("--E", type=float, default=1.0e3)
    parser.add_argument("--nu", type=float, default=0.3)
    parser.add_argument("--eps", type=float, default=0.01, help="Prescribed uniform strain in x: ux = eps * x")
    parser.add_argument("--out", type=str, default="data/verify_h8_consistent.dat")
    args = parser.parse_args()

    nodes = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
    ]

    K = build_k_matrix(nodes, args.E, args.nu)

    u = [0.0] * 24
    for i, (x, y, z) in enumerate(nodes):
        u[3 * i] = args.eps * x
        u[3 * i + 1] = 0.0
        u[3 * i + 2] = 0.0

    f = mat_vec(K, u)

    # Boundary conditions: fix node1 (ux,uy,uz); fix node2 (uy,uz); fix node4 (uz)
    bcs = {
        1: (1, 1, 1),
        2: (0, 1, 1),
        3: (0, 0, 0),
        4: (1, 0, 1),
        5: (0, 0, 0),
        6: (0, 0, 0),
        7: (0, 0, 0),
        8: (0, 0, 0),
    }

    loads = []
    for i in range(8):
        nid = i + 1
        bc = bcs.get(nid, (0, 0, 0))
        for d in range(3):
            if bc[d] == 0:
                val = f[i * 3 + d]
                if abs(val) > 1e-12:
                    loads.append((nid, d + 1, val))

    with open(args.out, "w", encoding="utf-8") as out:
        out.write("Hex8_Validation_Consistent\n")
        out.write("8 1 1 1\n")
        for i, (x, y, z) in enumerate(nodes, start=1):
            bc = bcs.get(i, (0, 0, 0))
            out.write(f"{i} {bc[0]} {bc[1]} {bc[2]} {x} {y} {z}\n")
        out.write("1\n")
        out.write(f"{len(loads)} 0 0 0 0\n")
        for nid, dof, val in loads:
            out.write(f"{nid} {dof} {val}\n")
        out.write("4 1 1\n")
        out.write(f"1 {args.E} {args.nu} 0.0\n")
        out.write("1 1 2 3 4 5 6 7 8 1\n")

    D = build_d_matrix(args.E, args.nu)
    strain = [args.eps, 0.0, 0.0, 0.0, 0.0, 0.0]
    stress = mat_vec(D, strain)

    print("Generated:", args.out)
    print("Expected stress (SXX, SYY, SZZ, SXY, SYZ, SXZ):")
    print(" ".join(f"{v:.6e}" for v in stress))


if __name__ == "__main__":
    main()
