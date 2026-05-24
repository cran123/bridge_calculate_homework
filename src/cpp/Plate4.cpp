/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Plate4.h"

#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

namespace
{
	const double drillingStabilization = 1.0e-6;
	const double shearCorrection = 5.0 / 6.0;

	void ShapeFunction(double xi, double eta, double N[4], double dNdxi[4], double dNdeta[4])
	{
		N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
		N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
		N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
		N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);

		dNdxi[0] = -0.25 * (1.0 - eta);
		dNdxi[1] =  0.25 * (1.0 - eta);
		dNdxi[2] =  0.25 * (1.0 + eta);
		dNdxi[3] = -0.25 * (1.0 + eta);

		dNdeta[0] = -0.25 * (1.0 - xi);
		dNdeta[1] = -0.25 * (1.0 + xi);
		dNdeta[2] =  0.25 * (1.0 + xi);
		dNdeta[3] =  0.25 * (1.0 - xi);
	}

	bool ShapeDerivativeXY(CNode** nodes, double xi, double eta, double N[4], double dNdx[4], double dNdy[4], double& detJ)
	{
		double dNdxi[4];
		double dNdeta[4];
		ShapeFunction(xi, eta, N, dNdxi, dNdeta);

		double j11 = 0.0;
		double j12 = 0.0;
		double j21 = 0.0;
		double j22 = 0.0;

		for (unsigned int i = 0; i < 4; i++)
		{
			j11 += dNdxi[i] * nodes[i]->XYZ[0];
			j12 += dNdxi[i] * nodes[i]->XYZ[1];
			j21 += dNdeta[i] * nodes[i]->XYZ[0];
			j22 += dNdeta[i] * nodes[i]->XYZ[1];
		}

		detJ = j11 * j22 - j12 * j21;
		if (fabs(detJ) < 1.0e-14)
			return false;

		double invJ11 =  j22 / detJ;
		double invJ12 = -j12 / detJ;
		double invJ21 = -j21 / detJ;
		double invJ22 =  j11 / detJ;

		for (unsigned int i = 0; i < 4; i++)
		{
			dNdx[i] = invJ11 * dNdxi[i] + invJ12 * dNdeta[i];
			dNdy[i] = invJ21 * dNdxi[i] + invJ22 * dNdeta[i];
		}

		return true;
	}

	void AddPacked(double* Matrix, const double full[24][24])
	{
		unsigned int index = 0;
		for (unsigned int j = 0; j < 24; j++)
		{
			Matrix[index++] = full[j][j];
			for (int i = static_cast<int>(j) - 1; i >= 0; i--)
				Matrix[index++] = full[i][j];
		}
	}

	void AddSymmetric(double full[24][24], unsigned int row, unsigned int col, double value)
	{
		full[row][col] += value;
		if (row != col)
			full[col][row] += value;
	}
}

//	Read plate material data from stream Input
bool CPlateMaterial::Read(ifstream& Input)
{
	Input >> nset;
	Input >> E >> Nu >> Thickness >> Density;

	return true;
}

//	Write plate material data to stream
void CPlateMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Nu << setw(16) << Thickness
		   << setw(16) << Density << endl;
}

//	Constructor
CPlate4::CPlate4()
{
	NEN_ = 4;
	nodes_ = new CNode*[NEN_];

	ND_ = 24;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Destructor
CPlate4::~CPlate4()
{
}

//	Read element data from stream Input
bool CPlate4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;
	unsigned int N1, N2, N3, N4;

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;

	ElementMaterial_ = dynamic_cast<CPlateMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//	Write element data to stream
void CPlate4::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(9) << nodes_[2]->NodeNumber
		   << setw(9) << nodes_[3]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

// Generate location matrix
void CPlate4::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
		for (unsigned int D = 0; D < CNode::NDF; D++)
			LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//	Calculate element stiffness matrix
void CPlate4::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	double k[24][24] = {};

	CPlateMaterial* material = dynamic_cast<CPlateMaterial*>(ElementMaterial_);
	const double E = material->E;
	const double nu = material->Nu;
	const double t = material->Thickness;
	const double G = E / (2.0 * (1.0 + nu));
	const double membraneFactor = E * t / (1.0 - nu * nu);
	const double bendingFactor = E * t * t * t / (12.0 * (1.0 - nu * nu));
	const double shearFactor = shearCorrection * G * t;

	const double Dm[3][3] = {
		{membraneFactor, membraneFactor * nu, 0.0},
		{membraneFactor * nu, membraneFactor, 0.0},
		{0.0, 0.0, membraneFactor * (1.0 - nu) * 0.5}
	};

	const double Db[3][3] = {
		{bendingFactor, bendingFactor * nu, 0.0},
		{bendingFactor * nu, bendingFactor, 0.0},
		{0.0, 0.0, bendingFactor * (1.0 - nu) * 0.5}
	};

	const double gauss = 1.0 / sqrt(3.0);
	const double points[2] = {-gauss, gauss};

	for (unsigned int gpEta = 0; gpEta < 2; gpEta++)
	{
		for (unsigned int gpXi = 0; gpXi < 2; gpXi++)
		{
			double N[4], dNdx[4], dNdy[4], detJ;
			if (!ShapeDerivativeXY(nodes_, points[gpXi], points[gpEta], N, dNdx, dNdy, detJ))
			{
				cerr << "*** Error *** Invalid Plate4 element Jacobian." << endl;
				return;
			}

			double Bm[3][24] = {};
			double Bb[3][24] = {};

			for (unsigned int a = 0; a < 4; a++)
			{
				const unsigned int base = 6 * a;

				Bm[0][base + 0] = dNdx[a];
				Bm[1][base + 1] = dNdy[a];
				Bm[2][base + 0] = dNdy[a];
				Bm[2][base + 1] = dNdx[a];

				Bb[0][base + 4] = dNdx[a];
				Bb[1][base + 3] = dNdy[a];
				Bb[2][base + 3] = dNdx[a];
				Bb[2][base + 4] = dNdy[a];
			}

			for (unsigned int i = 0; i < 24; i++)
			{
				for (unsigned int j = i; j < 24; j++)
				{
					double value = 0.0;
					for (unsigned int m = 0; m < 3; m++)
						for (unsigned int n = 0; n < 3; n++)
							value += Bm[m][i] * Dm[m][n] * Bm[n][j]
								   + Bb[m][i] * Db[m][n] * Bb[n][j];

					AddSymmetric(k, i, j, value * detJ);
				}
			}
		}
	}

	double N[4], dNdx[4], dNdy[4], detJ;
	if (!ShapeDerivativeXY(nodes_, 0.0, 0.0, N, dNdx, dNdy, detJ))
	{
		cerr << "*** Error *** Invalid Plate4 element Jacobian." << endl;
		return;
	}

	double Bs[2][24] = {};
	for (unsigned int a = 0; a < 4; a++)
	{
		const unsigned int base = 6 * a;

		Bs[0][base + 2] = dNdx[a];
		Bs[0][base + 4] = N[a];
		Bs[1][base + 2] = dNdy[a];
		Bs[1][base + 3] = N[a];
	}

	for (unsigned int i = 0; i < 24; i++)
	{
		for (unsigned int j = i; j < 24; j++)
		{
			double value = shearFactor * (Bs[0][i] * Bs[0][j] + Bs[1][i] * Bs[1][j]);
			AddSymmetric(k, i, j, value * detJ * 4.0);
		}
	}

	double maxDiagonal = 0.0;
	for (unsigned int i = 0; i < 24; i++)
		if (fabs(k[i][i]) > maxDiagonal)
			maxDiagonal = fabs(k[i][i]);

	if (maxDiagonal > 0.0)
	{
		const double drilling = drillingStabilization * maxDiagonal;
		for (unsigned int a = 0; a < 4; a++)
			k[6 * a + 5][6 * a + 5] += drilling;
	}

	AddPacked(Matrix, k);
}

//	Calculate element stress
void CPlate4::ElementStress(double* stress, double* Displacement)
{
	clear(stress, 7);
}
