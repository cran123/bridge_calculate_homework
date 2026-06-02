/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Hex8.h"
#include "Material.h"

#include <cmath>
#include <iomanip>

using namespace std;

namespace
{
    const bool kUseReducedIntegration = false; // Set true for C3D8R-style integration

    void ComputeShape(double r, double s, double t, double N[8], double dNdr[8], double dNds[8], double dNdt[8])
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

    void BuildDMatrix(double E, double nu, double D[6][6])
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
}

//	Constructor
CHex8::CHex8()
{
    NEN_ = 8;
    nodes_ = new CNode*[NEN_];

    ND_ = 24;
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

//	Desconstructor
CHex8::~CHex8()
{
}

//	Read element data from stream Input
bool CHex8::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;
    unsigned int N[8];

    Input >> N[0] >> N[1] >> N[2] >> N[3] >> N[4] >> N[5] >> N[6] >> N[7] >> MSet;

    ElementMaterial_ = dynamic_cast<CHex8Material*>(MaterialSets) + MSet - 1;

    for (unsigned int i = 0; i < 8; i++)
        nodes_[i] = &NodeList[N[i] - 1];

    return true;
}

//	Write element data to stream
void CHex8::Write(COutputter& output)
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

// Generate location matrix
void CHex8::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//	Calculate element stiffness matrix
void CHex8::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    double K[24][24];
    for (int i = 0; i < 24; i++)
        for (int j = 0; j < 24; j++)
            K[i][j] = 0.0;

    CHex8Material* material_ = dynamic_cast<CHex8Material*>(ElementMaterial_);

    double D[6][6];
    BuildDMatrix(material_->E, material_->Nu, D);

    double points[2] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};

    int nPoint = kUseReducedIntegration ? 1 : 2;

    for (int ir = 0; ir < nPoint; ir++)
    for (int is = 0; is < nPoint; is++)
    for (int it = 0; it < nPoint; it++)
    {
        double r = kUseReducedIntegration ? 0.0 : points[ir];
        double s = kUseReducedIntegration ? 0.0 : points[is];
        double t = kUseReducedIntegration ? 0.0 : points[it];
        double w = kUseReducedIntegration ? 8.0 : 1.0;

        double N[8], dNdr[8], dNds[8], dNdt[8];
        ComputeShape(r, s, t, N, dNdr, dNds, dNdt);

        double J[3][3], invJ[3][3], detJ;
        const CNode* nodePtr[8] = {nodes_[0], nodes_[1], nodes_[2], nodes_[3], nodes_[4], nodes_[5], nodes_[6], nodes_[7]};
        if (!ComputeJacobian(nodePtr, dNdr, dNds, dNdt, J, invJ, detJ))
            continue;

        double dNdx[8], dNdy[8], dNdz[8];
        for (int i = 0; i < 8; i++)
        {
            dNdx[i] = invJ[0][0] * dNdr[i] + invJ[0][1] * dNds[i] + invJ[0][2] * dNdt[i];
            dNdy[i] = invJ[1][0] * dNdr[i] + invJ[1][1] * dNds[i] + invJ[1][2] * dNdt[i];
            dNdz[i] = invJ[2][0] * dNdr[i] + invJ[2][1] * dNds[i] + invJ[2][2] * dNdt[i];
        }

        double B[6][24];
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 24; j++)
                B[i][j] = 0.0;

        for (int i = 0; i < 8; i++)
        {
            int col = i * 3;
            B[0][col]     = dNdx[i];
            B[1][col + 1] = dNdy[i];
            B[2][col + 2] = dNdz[i];

            B[3][col]     = dNdy[i];
            B[3][col + 1] = dNdx[i];

            B[4][col + 1] = dNdz[i];
            B[4][col + 2] = dNdy[i];

            B[5][col]     = dNdz[i];
            B[5][col + 2] = dNdx[i];
        }

        double DB[6][24];
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 24; j++)
            {
                DB[i][j] = 0.0;
                for (int k = 0; k < 6; k++)
                    DB[i][j] += D[i][k] * B[k][j];
            }

        double scale = detJ * w;
        for (int i = 0; i < 24; i++)
            for (int j = 0; j < 24; j++)
            {
                double sum = 0.0;
                for (int k = 0; k < 6; k++)
                    sum += B[k][i] * DB[k][j];
                K[i][j] += sum * scale;
            }
    }

    int idx = 0;
    for (int col = 0; col < 24; col++)
        for (int row = 0; row <= col; row++)
            Matrix[idx++] = K[row][col];
}

//	Calculate element stress (average at integration points)
void CHex8::ElementStress(double* stress, double* Displacement)
{
    for (int i = 0; i < 6; i++)
        stress[i] = 0.0;

    CHex8Material* material_ = dynamic_cast<CHex8Material*>(ElementMaterial_);

    double D[6][6];
    BuildDMatrix(material_->E, material_->Nu, D);

    double ue[24];
    for (int i = 0; i < 24; i++)
        ue[i] = 0.0;

    for (int i = 0; i < 24; i++)
        if (LocationMatrix_[i])
            ue[i] = Displacement[LocationMatrix_[i] - 1];

    double points[2] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
    int nPoint = kUseReducedIntegration ? 1 : 2;

    double totalWeight = 0.0;

    for (int ir = 0; ir < nPoint; ir++)
    for (int is = 0; is < nPoint; is++)
    for (int it = 0; it < nPoint; it++)
    {
        double r = kUseReducedIntegration ? 0.0 : points[ir];
        double s = kUseReducedIntegration ? 0.0 : points[is];
        double t = kUseReducedIntegration ? 0.0 : points[it];
        double w = kUseReducedIntegration ? 8.0 : 1.0;

        double N[8], dNdr[8], dNds[8], dNdt[8];
        ComputeShape(r, s, t, N, dNdr, dNds, dNdt);

        double J[3][3], invJ[3][3], detJ;
        const CNode* nodePtr[8] = {nodes_[0], nodes_[1], nodes_[2], nodes_[3], nodes_[4], nodes_[5], nodes_[6], nodes_[7]};
        if (!ComputeJacobian(nodePtr, dNdr, dNds, dNdt, J, invJ, detJ))
            continue;

        double dNdx[8], dNdy[8], dNdz[8];
        for (int i = 0; i < 8; i++)
        {
            dNdx[i] = invJ[0][0] * dNdr[i] + invJ[0][1] * dNds[i] + invJ[0][2] * dNdt[i];
            dNdy[i] = invJ[1][0] * dNdr[i] + invJ[1][1] * dNds[i] + invJ[1][2] * dNdt[i];
            dNdz[i] = invJ[2][0] * dNdr[i] + invJ[2][1] * dNds[i] + invJ[2][2] * dNdt[i];
        }

        double B[6][24];
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 24; j++)
                B[i][j] = 0.0;

        for (int i = 0; i < 8; i++)
        {
            int col = i * 3;
            B[0][col]     = dNdx[i];
            B[1][col + 1] = dNdy[i];
            B[2][col + 2] = dNdz[i];

            B[3][col]     = dNdy[i];
            B[3][col + 1] = dNdx[i];

            B[4][col + 1] = dNdz[i];
            B[4][col + 2] = dNdy[i];

            B[5][col]     = dNdz[i];
            B[5][col + 2] = dNdx[i];
        }

        double strain[6];
        for (int i = 0; i < 6; i++)
        {
            strain[i] = 0.0;
            for (int j = 0; j < 24; j++)
                strain[i] += B[i][j] * ue[j];
        }

        double sigma[6];
        for (int i = 0; i < 6; i++)
        {
            sigma[i] = 0.0;
            for (int j = 0; j < 6; j++)
                sigma[i] += D[i][j] * strain[j];
        }

        double weight = detJ * w;
        for (int i = 0; i < 6; i++)
            stress[i] += sigma[i] * weight;

        totalWeight += weight;
    }

    if (totalWeight > 0.0)
        for (int i = 0; i < 6; i++)
            stress[i] /= totalWeight;
}

//	Calculate element body force vector (self-weight)
void CHex8::ElementBodyForce(double* bodyForce, const double gravity[3])
{
    clear(bodyForce, ND_);

    CHex8Material* material_ = dynamic_cast<CHex8Material*>(ElementMaterial_);

    double points[2] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
    int nPoint = kUseReducedIntegration ? 1 : 2;

    for (int ir = 0; ir < nPoint; ir++)
    for (int is = 0; is < nPoint; is++)
    for (int it = 0; it < nPoint; it++)
    {
        double r = kUseReducedIntegration ? 0.0 : points[ir];
        double s = kUseReducedIntegration ? 0.0 : points[is];
        double t = kUseReducedIntegration ? 0.0 : points[it];
        double w = kUseReducedIntegration ? 8.0 : 1.0;

        double N[8], dNdr[8], dNds[8], dNdt[8];
        ComputeShape(r, s, t, N, dNdr, dNds, dNdt);

        double J[3][3], invJ[3][3], detJ;
        const CNode* nodePtr[8] = {nodes_[0], nodes_[1], nodes_[2], nodes_[3], nodes_[4], nodes_[5], nodes_[6], nodes_[7]};
        if (!ComputeJacobian(nodePtr, dNdr, dNds, dNdt, J, invJ, detJ))
            continue;

        double scale = material_->Density * detJ * w;
        for (int i = 0; i < 8; i++)
        {
            double Ni = N[i];
            int col = i * 3;
            bodyForce[col]     += Ni * scale * gravity[0];
            bodyForce[col + 1] += Ni * scale * gravity[1];
            bodyForce[col + 2] += Ni * scale * gravity[2];
        }
    }
}
