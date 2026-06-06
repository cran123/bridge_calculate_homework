/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>

using namespace std;

static void BuildDMatrix(double E, double nu, double D[6][6])
{
    double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double a = (1.0 - nu) * factor;
    double b = nu * factor;
    double c = 0.5 * (1.0 - 2.0 * nu) * factor;

    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            D[i][j] = 0.0;

    D[0][0] = a; D[0][1] = b; D[0][2] = b;
    D[1][0] = b; D[1][1] = a; D[1][2] = b;
    D[2][0] = b; D[2][1] = b; D[2][2] = a;

    D[3][3] = c;
    D[4][4] = c;
    D[5][5] = c;
}

static void ComputeShape(double r, double s, double t, double N[8], double dNdr[8], double dNds[8], double dNdt[8])
{
    double rp = 1.0 + r;
    double rm = 1.0 - r;
    double sp = 1.0 + s;
    double sm = 1.0 - s;
    double tp = 1.0 + t;
    double tm = 1.0 - t;

    N[0] = 0.125 * rm * sm * tm;
    N[1] = 0.125 * rp * sm * tm;
    N[2] = 0.125 * rp * sp * tm;
    N[3] = 0.125 * rm * sp * tm;
    N[4] = 0.125 * rm * sm * tp;
    N[5] = 0.125 * rp * sm * tp;
    N[6] = 0.125 * rp * sp * tp;
    N[7] = 0.125 * rm * sp * tp;

    dNdr[0] = -0.125 * sm * tm;
    dNdr[1] =  0.125 * sm * tm;
    dNdr[2] =  0.125 * sp * tm;
    dNdr[3] = -0.125 * sp * tm;
    dNdr[4] = -0.125 * sm * tp;
    dNdr[5] =  0.125 * sm * tp;
    dNdr[6] =  0.125 * sp * tp;
    dNdr[7] = -0.125 * sp * tp;

    dNds[0] = -0.125 * rm * tm;
    dNds[1] = -0.125 * rp * tm;
    dNds[2] =  0.125 * rp * tm;
    dNds[3] =  0.125 * rm * tm;
    dNds[4] = -0.125 * rm * tp;
    dNds[5] = -0.125 * rp * tp;
    dNds[6] =  0.125 * rp * tp;
    dNds[7] =  0.125 * rm * tp;

    dNdt[0] = -0.125 * rm * sm;
    dNdt[1] = -0.125 * rp * sm;
    dNdt[2] = -0.125 * rp * sp;
    dNdt[3] = -0.125 * rm * sp;
    dNdt[4] =  0.125 * rm * sm;
    dNdt[5] =  0.125 * rp * sm;
    dNdt[6] =  0.125 * rp * sp;
    dNdt[7] =  0.125 * rm * sp;
}

static bool ComputeJacobian(const double nodes[8][3],
    const double dNdr[8], const double dNds[8], const double dNdt[8],
    double J[3][3], double invJ[3][3], double& detJ)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            J[i][j] = 0.0;

    for (int i = 0; i < 8; i++)
    {
        J[0][0] += dNdr[i] * nodes[i][0];
        J[0][1] += dNds[i] * nodes[i][0];
        J[0][2] += dNdt[i] * nodes[i][0];

        J[1][0] += dNdr[i] * nodes[i][1];
        J[1][1] += dNds[i] * nodes[i][1];
        J[1][2] += dNdt[i] * nodes[i][1];

        J[2][0] += dNdr[i] * nodes[i][2];
        J[2][1] += dNds[i] * nodes[i][2];
        J[2][2] += dNdt[i] * nodes[i][2];
    }

    detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
         - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
         + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

    if (fabs(detJ) <= 1e-12)
        return false;

    double invDet = 1.0 / detJ;

    invJ[0][0] =  (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * invDet;
    invJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * invDet;
    invJ[0][2] =  (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * invDet;

    invJ[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]) * invDet;
    invJ[1][1] =  (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * invDet;
    invJ[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]) * invDet;

    invJ[2][0] =  (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * invDet;
    invJ[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]) * invDet;
    invJ[2][2] =  (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * invDet;

    return true;
}

int main(int argc, char* argv[])
{
    double E = 1.0e3;
    double nu = 0.3;
    double eps = 0.01;
    string outFile = "..\\data\\verify_h8_consistent.dat";

    if (argc > 1)
        outFile = argv[1];

    const double nodes[8][3] = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {1.0, 1.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 0.0, 1.0},
        {1.0, 1.0, 1.0},
        {0.0, 1.0, 1.0}
    };

    double D[6][6];
    BuildDMatrix(E, nu, D);

    double K[24][24] = {0.0};

    double gp[2] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
    for (int ir = 0; ir < 2; ir++)
    for (int is = 0; is < 2; is++)
    for (int it = 0; it < 2; it++)
    {
        double N[8], dNdr[8], dNds[8], dNdt[8];
        ComputeShape(gp[ir], gp[is], gp[it], N, dNdr, dNds, dNdt);

        double J[3][3], invJ[3][3], detJ;
        if (!ComputeJacobian(nodes, dNdr, dNds, dNdt, J, invJ, detJ))
            continue;

        double dNdx[8], dNdy[8], dNdz[8];
        for (int i = 0; i < 8; i++)
        {
            dNdx[i] = invJ[0][0] * dNdr[i] + invJ[0][1] * dNds[i] + invJ[0][2] * dNdt[i];
            dNdy[i] = invJ[1][0] * dNdr[i] + invJ[1][1] * dNds[i] + invJ[1][2] * dNdt[i];
            dNdz[i] = invJ[2][0] * dNdr[i] + invJ[2][1] * dNds[i] + invJ[2][2] * dNdt[i];
        }

        double B[6][24] = {0.0};
        for (int i = 0; i < 8; i++)
        {
            int col = i * 3;
            B[0][col] = dNdx[i];
            B[1][col + 1] = dNdy[i];
            B[2][col + 2] = dNdz[i];

            B[3][col] = dNdy[i];
            B[3][col + 1] = dNdx[i];

            B[4][col + 1] = dNdz[i];
            B[4][col + 2] = dNdy[i];

            B[5][col] = dNdz[i];
            B[5][col + 2] = dNdx[i];
        }

        double DB[6][24] = {0.0};
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 24; j++)
                for (int k = 0; k < 6; k++)
                    DB[i][j] += D[i][k] * B[k][j];

        for (int i = 0; i < 24; i++)
            for (int j = 0; j < 24; j++)
            {
                double sum = 0.0;
                for (int k = 0; k < 6; k++)
                    sum += B[k][i] * DB[k][j];
                K[i][j] += sum * detJ;
            }
    }

    double u[24] = {0.0};
    for (int i = 0; i < 8; i++)
    {
        u[i * 3] = eps * nodes[i][0];
        u[i * 3 + 1] = -nu * eps * nodes[i][1];
        u[i * 3 + 2] = -nu * eps * nodes[i][2];
    }

    double f[24] = {0.0};
    for (int i = 0; i < 24; i++)
    {
        for (int j = 0; j < 24; j++)
            f[i] += K[i][j] * u[j];
    }

    int bc[8][3] = {
        {1, 1, 1},
        {0, 0, 0},
        {0, 0, 0},
        {1, 0, 1},
        {1, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {1, 0, 0}
    };

    vector<tuple<int, int, double>> loads;
    for (int i = 0; i < 8; i++)
    {
        for (int d = 0; d < 3; d++)
        {
            if (bc[i][d] == 0 && fabs(f[i * 3 + d]) > 1e-12)
                loads.emplace_back(i + 1, d + 1, f[i * 3 + d]);
        }
    }

    ofstream out(outFile);
    if (!out)
    {
        cerr << "*** Error *** Failed to write " << outFile << endl;
        return 1;
    }

    out << "Hex8_Validation_Consistent\n";
    out << "8 1 1 1\n";
    out << setprecision(6) << scientific;
    for (int i = 0; i < 8; i++)
    {
        out << (i + 1) << " " << bc[i][0] << " " << bc[i][1] << " " << bc[i][2] << " "
            << nodes[i][0] << " " << nodes[i][1] << " " << nodes[i][2] << "\n";
    }
    out << "1\n";
    out << loads.size() << " 0 0 0 0\n";
    for (auto& load : loads)
        out << get<0>(load) << " " << get<1>(load) << " " << get<2>(load) << "\n";
    out << "4 1 1\n";
    out << "1 " << E << " " << nu << " 0.0\n";
    out << "1 1 2 3 4 5 6 7 8 1\n";

    double strain[6] = {eps, -nu * eps, -nu * eps, 0.0, 0.0, 0.0};
    double stress[6] = {0.0};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            stress[i] += D[i][j] * strain[j];

    cout << "Generated: " << outFile << "\n";
    cout << "Expected stress (SXX SYY SZZ SXY SYZ SXZ):";
    for (int i = 0; i < 6; i++)
        cout << " " << stress[i];
    cout << endl;

    return 0;
}
