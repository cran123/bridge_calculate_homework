/*****************************************************************************/
/*  STAP++ : B-bar Hex8 element — B-bar method for reduced locking           */
/*  Reference: OpenSees BbarBrick.cpp                                        */
/*****************************************************************************/

#include "BbarHex8.h"
#include "Material.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

namespace
{
    void ComputeShape(double r, double s, double t,
        double N[8], double dNdr[8], double dNds[8], double dNdt[8])
    {
        double rp = 1.0 + r, rm = 1.0 - r;
        double sp = 1.0 + s, sm = 1.0 - s;
        double tp = 1.0 + t, tm = 1.0 - t;

        N[0] = 0.125 * rm * sm * tm;
        N[1] = 0.125 * rp * sm * tm;
        N[2] = 0.125 * rp * sp * tm;
        N[3] = 0.125 * rm * sp * tm;
        N[4] = 0.125 * rm * sm * tp;
        N[5] = 0.125 * rp * sm * tp;
        N[6] = 0.125 * rp * sp * tp;
        N[7] = 0.125 * rm * sp * tp;

        dNdr[0] = -0.125 * sm * tm;  dNdr[1] =  0.125 * sm * tm;
        dNdr[2] =  0.125 * sp * tm;  dNdr[3] = -0.125 * sp * tm;
        dNdr[4] = -0.125 * sm * tp;  dNdr[5] =  0.125 * sm * tp;
        dNdr[6] =  0.125 * sp * tp;  dNdr[7] = -0.125 * sp * tp;

        dNds[0] = -0.125 * rm * tm;  dNds[1] = -0.125 * rp * tm;
        dNds[2] =  0.125 * rp * tm;  dNds[3] =  0.125 * rm * tm;
        dNds[4] = -0.125 * rm * tp;  dNds[5] = -0.125 * rp * tp;
        dNds[6] =  0.125 * rp * tp;  dNds[7] =  0.125 * rm * tp;

        dNdt[0] = -0.125 * rm * sm;  dNdt[1] = -0.125 * rp * sm;
        dNdt[2] = -0.125 * rp * sp;  dNdt[3] = -0.125 * rm * sp;
        dNdt[4] =  0.125 * rm * sm;  dNdt[5] =  0.125 * rp * sm;
        dNdt[6] =  0.125 * rp * sp;  dNdt[7] =  0.125 * rm * sp;
    }

    bool ComputeJacobian(const CNode* const nodes[8],
        const double dNdr[8], const double dNds[8], const double dNdt[8],
        double J[3][3], double invJ[3][3], double& detJ)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                J[i][j] = 0.0;

        for (int i = 0; i < 8; i++)
        {
            J[0][0] += dNdr[i] * nodes[i]->XYZ[0];
            J[0][1] += dNds[i] * nodes[i]->XYZ[0];
            J[0][2] += dNdt[i] * nodes[i]->XYZ[0];
            J[1][0] += dNdr[i] * nodes[i]->XYZ[1];
            J[1][1] += dNds[i] * nodes[i]->XYZ[1];
            J[1][2] += dNdt[i] * nodes[i]->XYZ[1];
            J[2][0] += dNdr[i] * nodes[i]->XYZ[2];
            J[2][1] += dNds[i] * nodes[i]->XYZ[2];
            J[2][2] += dNdt[i] * nodes[i]->XYZ[2];
        }

        detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
             - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
             + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

        if (fabs(detJ) <= 1e-12) return false;

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

    void BuildDMatrix(double E, double nu, double D[6][6])
    {
        double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        double a = (1.0 - nu) * factor;
        double b = nu * factor;
        double c = 0.5 * (1.0 - 2.0 * nu) * factor;
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++) D[i][j] = 0.0;
        D[0][0] = a; D[0][1] = b; D[0][2] = b;
        D[1][0] = b; D[1][1] = a; D[1][2] = b;
        D[2][0] = b; D[2][1] = b; D[2][2] = a;
        D[3][3] = c; D[4][4] = c; D[5][5] = c;
    }

    // Compute B-bar matrix for a single node (6×3)
    // shpBar[0]=avg dNdx, shpBar[1]=avg dNdy, shpBar[2]=avg dNdz, shpBar[3]=avg N
    void ComputeBbar(int node,
        const double dNdx[8], const double dNdy[8], const double dNdz[8],
        const double shpBar[4][8], double Bbar[6][3])
    {
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 3; j++) Bbar[i][j] = 0.0;

        const double one3 = 1.0 / 3.0;

        // Deviatoric B (using Gauss point derivatives)
        double Bdev[3][3];
        Bdev[0][0] =  2.0 * dNdx[node]; Bdev[0][1] = -dNdy[node];      Bdev[0][2] = -dNdz[node];
        Bdev[1][0] = -dNdx[node];       Bdev[1][1] =  2.0 * dNdy[node]; Bdev[1][2] = -dNdz[node];
        Bdev[2][0] = -dNdx[node];       Bdev[2][1] = -dNdy[node];       Bdev[2][2] =  2.0 * dNdz[node];

        // Volumetric B (using element-average derivatives)
        double Bvol[3][3];
        Bvol[0][0] = shpBar[0][node]; Bvol[0][1] = shpBar[1][node]; Bvol[0][2] = shpBar[2][node];
        Bvol[1][0] = shpBar[0][node]; Bvol[1][1] = shpBar[1][node]; Bvol[1][2] = shpBar[2][node];
        Bvol[2][0] = shpBar[0][node]; Bvol[2][1] = shpBar[1][node]; Bvol[2][2] = shpBar[2][node];

        // Extensional: Bbar = (Bdev + Bvol) / 3
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Bbar[i][j] = one3 * (Bdev[i][j] + Bvol[i][j]);

        // Shear terms (standard B using Gauss point derivatives)
        Bbar[3][0] = dNdy[node]; Bbar[3][1] = dNdx[node];
        Bbar[4][1] = dNdz[node]; Bbar[4][2] = dNdy[node];
        Bbar[5][0] = dNdz[node]; Bbar[5][2] = dNdx[node];
    }
}

CBbarHex8::CBbarHex8()
{
    NEN_ = 8;
    nodes_ = new CNode*[NEN_];
    ND_ = 24;
    LocationMatrix_ = new unsigned int[ND_];
    ElementMaterial_ = nullptr;
}

CBbarHex8::~CBbarHex8() {}

bool CBbarHex8::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet, N[8];
    Input >> N[0] >> N[1] >> N[2] >> N[3] >> N[4] >> N[5] >> N[6] >> N[7] >> MSet;
    ElementMaterial_ = dynamic_cast<CHex8Material*>(MaterialSets) + MSet - 1;
    for (unsigned int i = 0; i < 8; i++)
        nodes_[i] = &NodeList[N[i] - 1];
    return true;
}

void CBbarHex8::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
           << setw(9) << nodes_[1]->NodeNumber
           << setw(9) << nodes_[2]->NodeNumber
           << setw(9) << nodes_[3]->NodeNumber
           << setw(9) << nodes_[4]->NodeNumber
           << setw(9) << nodes_[5]->NodeNumber
           << setw(9) << nodes_[6]->NodeNumber
           << setw(9) << nodes_[7]->NodeNumber
           << setw(12) << ElementMaterial_->nset << endl;
}

void CBbarHex8::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    double K[24][24] = {};
    for (int i = 0; i < 24; i++)
        for (int j = 0; j < 24; j++) K[i][j] = 0.0;

    CHex8Material* mat = dynamic_cast<CHex8Material*>(ElementMaterial_);
    double D[6][6];
    BuildDMatrix(mat->E, mat->Nu, D);

    const double points[2] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};

    // Storage for 8 Gauss points
    double saved_dNdx[8][8], saved_dNdy[8][8], saved_dNdz[8][8];
    double saved_dvol[8];
    int gp = 0;

    // shpBar: [0]=avg dNdx, [1]=avg dNdy, [2]=avg dNdz, [3]=avg N
    double shpBar[4][8] = {};
    double volume = 0.0;

    // ---- First pass: save shape functions, accumulate volume & shpBar ----
    for (int ir = 0; ir < 2; ir++)
    for (int is = 0; is < 2; is++)
    for (int it = 0; it < 2; it++)
    {
        double r = points[ir], s = points[is], t = points[it];

        double N[8], dNdr[8], dNds[8], dNdt[8];
        ComputeShape(r, s, t, N, dNdr, dNds, dNdt);

        double J[3][3], invJ[3][3], detJ;
        const CNode* np[8] = {nodes_[0],nodes_[1],nodes_[2],nodes_[3],
                              nodes_[4],nodes_[5],nodes_[6],nodes_[7]};
        if (!ComputeJacobian(np, dNdr, dNds, dNdt, J, invJ, detJ)) continue;

        double dNdx[8], dNdy[8], dNdz[8];
        for (int i = 0; i < 8; i++)
        {
            dNdx[i] = invJ[0][0]*dNdr[i] + invJ[0][1]*dNds[i] + invJ[0][2]*dNdt[i];
            dNdy[i] = invJ[1][0]*dNdr[i] + invJ[1][1]*dNds[i] + invJ[1][2]*dNdt[i];
            dNdz[i] = invJ[2][0]*dNdr[i] + invJ[2][1]*dNds[i] + invJ[2][2]*dNdt[i];
        }

        double dvol = detJ * 1.0; // full integration weight = 1

        for (int i = 0; i < 8; i++)
        {
            saved_dNdx[gp][i] = dNdx[i];
            saved_dNdy[gp][i] = dNdy[i];
            saved_dNdz[gp][i] = dNdz[i];
        }
        saved_dvol[gp] = dvol;

        volume += dvol;
        for (int i = 0; i < 8; i++)
        {
            shpBar[0][i] += dvol * dNdx[i];
            shpBar[1][i] += dvol * dNdy[i];
            shpBar[2][i] += dvol * dNdz[i];
            shpBar[3][i] += dvol * N[i];
        }
        gp++;
    }

    // Normalize shpBar
    if (volume > 0.0)
        for (int p = 0; p < 4; p++)
            for (int i = 0; i < 8; i++)
                shpBar[p][i] /= volume;

    // ---- Second pass: integrate stiffness using Bbar ----
    for (int ig = 0; ig < 8; ig++)
    {
        double dvol = saved_dvol[ig];
        const double* dNdx = saved_dNdx[ig];
        const double* dNdy = saved_dNdy[ig];
        const double* dNdz = saved_dNdz[ig];

        for (int j = 0; j < 8; j++)
        {
            double Bj[6][3];
            ComputeBbar(j, dNdx, dNdy, dNdz, shpBar, Bj);

            // BjTD = Bj^T * D  (3×6)
            double BjTD[3][6] = {};
            for (int p = 0; p < 3; p++)
                for (int q = 0; q < 6; q++)
                    for (int r = 0; r < 6; r++)
                        BjTD[p][q] += Bj[r][p] * D[r][q];

            int jj = 3 * j;
            for (int k = 0; k < 8; k++)
            {
                double Bk[6][3];
                ComputeBbar(k, dNdx, dNdy, dNdz, shpBar, Bk);

                int kk = 3 * k;
                for (int p = 0; p < 3; p++)
                    for (int q = 0; q < 3; q++)
                    {
                        double sum = 0.0;
                        for (int r = 0; r < 6; r++)
                            sum += BjTD[p][r] * Bk[r][q];
                        K[jj + p][kk + q] += sum * dvol;
                    }
            }
        }
    }

    // Pack upper triangle (column-major, increasing row)
    int idx = 0;
    for (int col = 0; col < 24; col++)
        for (int row = 0; row <= col; row++)
            Matrix[idx++] = K[row][col];
}

void CBbarHex8::ElementStress(double* stress, double* Displacement)
{
    for (int i = 0; i < 6; i++) stress[i] = 0.0;

    CHex8Material* mat = dynamic_cast<CHex8Material*>(ElementMaterial_);
    double D[6][6];
    BuildDMatrix(mat->E, mat->Nu, D);

    // Extract element displacement
    double ue[24] = {};
    for (int i = 0; i < 24; i++)
        if (LocationMatrix_[i]) ue[i] = Displacement[LocationMatrix_[i] - 1];

    const double points[2] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};

    // First pass: compute shpBar
    double saved_dNdx[8][8], saved_dNdy[8][8], saved_dNdz[8][8];
    double saved_dvol[8];
    int gp = 0;
    double shpBar[4][8] = {};
    double volume = 0.0;

    for (int ir = 0; ir < 2; ir++)
    for (int is = 0; is < 2; is++)
    for (int it = 0; it < 2; it++)
    {
        double r = points[ir], s = points[is], t = points[it];
        double N[8], dNdr[8], dNds[8], dNdt[8];
        ComputeShape(r, s, t, N, dNdr, dNds, dNdt);

        double J[3][3], invJ[3][3], detJ;
        const CNode* np[8] = {nodes_[0],nodes_[1],nodes_[2],nodes_[3],
                              nodes_[4],nodes_[5],nodes_[6],nodes_[7]};
        if (!ComputeJacobian(np, dNdr, dNds, dNdt, J, invJ, detJ)) continue;

        double dNdx[8], dNdy[8], dNdz[8];
        for (int i = 0; i < 8; i++)
        {
            dNdx[i] = invJ[0][0]*dNdr[i] + invJ[0][1]*dNds[i] + invJ[0][2]*dNdt[i];
            dNdy[i] = invJ[1][0]*dNdr[i] + invJ[1][1]*dNds[i] + invJ[1][2]*dNdt[i];
            dNdz[i] = invJ[2][0]*dNdr[i] + invJ[2][1]*dNds[i] + invJ[2][2]*dNdt[i];
        }

        double dvol = detJ;
        for (int i = 0; i < 8; i++)
        {
            saved_dNdx[gp][i] = dNdx[i];
            saved_dNdy[gp][i] = dNdy[i];
            saved_dNdz[gp][i] = dNdz[i];
        }
        saved_dvol[gp] = dvol;
        volume += dvol;
        for (int i = 0; i < 8; i++)
        {
            shpBar[0][i] += dvol * dNdx[i];
            shpBar[1][i] += dvol * dNdy[i];
            shpBar[2][i] += dvol * dNdz[i];
            shpBar[3][i] += dvol * N[i];
        }
        gp++;
    }

    if (volume > 0.0)
        for (int p = 0; p < 4; p++)
            for (int i = 0; i < 8; i++)
                shpBar[p][i] /= volume;

    // Second pass: compute stress using Bbar
    double totalWeight = 0.0;
    for (int ig = 0; ig < 8; ig++)
    {
        const double* dNdx = saved_dNdx[ig];
        const double* dNdy = saved_dNdy[ig];
        const double* dNdz = saved_dNdz[ig];

        // Compute strain: ε = Σ Bbar_j * u_j
        double strain[6] = {};
        for (int j = 0; j < 8; j++)
        {
            double Bj[6][3];
            ComputeBbar(j, dNdx, dNdy, dNdz, shpBar, Bj);
            int jj = 3 * j;
            for (int r = 0; r < 6; r++)
                for (int c = 0; c < 3; c++)
                    strain[r] += Bj[r][c] * ue[jj + c];
        }

        double sigma[6] = {};
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                sigma[i] += D[i][j] * strain[j];

        double w = saved_dvol[ig];
        for (int i = 0; i < 6; i++) stress[i] += sigma[i] * w;
        totalWeight += w;
    }

    if (totalWeight > 0.0)
        for (int i = 0; i < 6; i++) stress[i] /= totalWeight;
}

void CBbarHex8::ElementBodyForce(double* bodyForce, const double gravity[3])
{
    clear(bodyForce, ND_);

    CHex8Material* mat = dynamic_cast<CHex8Material*>(ElementMaterial_);

    const double points[2] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};

    for (int ir = 0; ir < 2; ir++)
    for (int is = 0; is < 2; is++)
    for (int it = 0; it < 2; it++)
    {
        double r = points[ir], s = points[is], t = points[it];
        double N[8], dNdr[8], dNds[8], dNdt[8];
        ComputeShape(r, s, t, N, dNdr, dNds, dNdt);

        double J[3][3], invJ[3][3], detJ;
        const CNode* np[8] = {nodes_[0],nodes_[1],nodes_[2],nodes_[3],
                              nodes_[4],nodes_[5],nodes_[6],nodes_[7]};
        if (!ComputeJacobian(np, dNdr, dNds, dNdt, J, invJ, detJ)) continue;

        double scale = mat->Density * detJ * 1.0;
        for (int i = 0; i < 8; i++)
        {
            int col = 3 * i;
            bodyForce[col]     += N[i] * scale * gravity[0];
            bodyForce[col + 1] += N[i] * scale * gravity[1];
            bodyForce[col + 2] += N[i] * scale * gravity[2];
        }
    }
}
