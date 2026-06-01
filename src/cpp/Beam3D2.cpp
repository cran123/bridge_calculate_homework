/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*****************************************************************************/

#include "Beam3D2.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace std;

namespace
{
	const double EPS = 1.0e-12;

	void ClearMatrix(double matrix[12][12])
	{
		for (unsigned int i = 0; i < 12; i++)
			for (unsigned int j = 0; j < 12; j++)
				matrix[i][j] = 0.0;
	}

	void SetSymmetric(double matrix[12][12], unsigned int i, unsigned int j, double value)
	{
		matrix[i][j] = value;
		matrix[j][i] = value;
	}

	double Dot(const double a[3], const double b[3])
	{
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}

	void Cross(const double a[3], const double b[3], double c[3])
	{
		c[0] = a[1] * b[2] - a[2] * b[1];
		c[1] = a[2] * b[0] - a[0] * b[2];
		c[2] = a[0] * b[1] - a[1] * b[0];
	}

	double Norm(double a[3])
	{
		return sqrt(Dot(a, a));
	}

	bool Normalize(double a[3])
	{
		double norm = Norm(a);
		if (norm < EPS)
			return false;

		for (unsigned int i = 0; i < 3; i++)
			a[i] /= norm;

		return true;
	}

	void FillBendingZ(double matrix[12][12], double EI, double L)
	{
		unsigned int dof[4] = {1, 5, 7, 11};
		double L2 = L * L;
		double L3 = L2 * L;
		double k[4][4] = {
			{ 12.0 * EI / L3,  6.0 * EI / L2, -12.0 * EI / L3,  6.0 * EI / L2 },
			{  6.0 * EI / L2,  4.0 * EI / L,  -6.0 * EI / L2,  2.0 * EI / L },
			{-12.0 * EI / L3, -6.0 * EI / L2,  12.0 * EI / L3, -6.0 * EI / L2 },
			{  6.0 * EI / L2,  2.0 * EI / L,  -6.0 * EI / L2,  4.0 * EI / L }
		};

		for (unsigned int i = 0; i < 4; i++)
			for (unsigned int j = i; j < 4; j++)
				SetSymmetric(matrix, dof[i], dof[j], k[i][j]);
	}

	void FillBendingY(double matrix[12][12], double EI, double L)
	{
		unsigned int dof[4] = {2, 4, 8, 10};
		double L2 = L * L;
		double L3 = L2 * L;
		double k[4][4] = {
			{ 12.0 * EI / L3, -6.0 * EI / L2, -12.0 * EI / L3, -6.0 * EI / L2 },
			{ -6.0 * EI / L2,  4.0 * EI / L,   6.0 * EI / L2,  2.0 * EI / L },
			{-12.0 * EI / L3,  6.0 * EI / L2,  12.0 * EI / L3,  6.0 * EI / L2 },
			{ -6.0 * EI / L2,  2.0 * EI / L,   6.0 * EI / L2,  4.0 * EI / L }
		};

		for (unsigned int i = 0; i < 4; i++)
			for (unsigned int j = i; j < 4; j++)
				SetSymmetric(matrix, dof[i], dof[j], k[i][j]);
	}

	void PackUpperTriangle(double full[12][12], double* packed)
	{
		for (unsigned int j = 0; j < 12; j++)
			for (unsigned int i = 0; i <= j; i++)
				packed[j * (j + 1) / 2 + j - i] = full[i][j];
	}
}

CBeam3D2::CBeam3D2()
{
	NEN_ = 2;
	nodes_ = new CNode*[NEN_];

	ND_ = 12;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
	Reference_[0] = 0.0;
	Reference_[1] = 0.0;
	Reference_[2] = 1.0;
}

CBeam3D2::~CBeam3D2()
{
}

bool CBeam3D2::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;
	unsigned int N1, N2;

	Input >> N1 >> N2 >> MSet >> Reference_[0] >> Reference_[1] >> Reference_[2];
	ElementMaterial_ = dynamic_cast<CBeamMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

void CBeam3D2::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(12) << ElementMaterial_->nset
		   << setw(16) << Reference_[0]
		   << setw(16) << Reference_[1]
		   << setw(16) << Reference_[2] << endl;
}

void CBeam3D2::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
		for (unsigned int D = 0; D < CNode::NDF; D++)
			LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

double CBeam3D2::Length() const
{
	double dx[3];
	double length2 = 0.0;
	for (unsigned int i = 0; i < 3; i++)
	{
		dx[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		length2 += dx[i] * dx[i];
	}

	return sqrt(length2);
}

bool CBeam3D2::LocalAxes(double axes[3][3]) const
{
	for (unsigned int i = 0; i < 3; i++)
		axes[0][i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	if (!Normalize(axes[0]))
		return false;

	double projection = Dot(Reference_, axes[0]);
	for (unsigned int i = 0; i < 3; i++)
		axes[1][i] = Reference_[i] - projection * axes[0][i];

	if (!Normalize(axes[1]))
	{
		double fallback[3] = {0.0, 0.0, 1.0};
		if (fabs(Dot(fallback, axes[0])) > 0.9)
		{
			fallback[0] = 0.0;
			fallback[1] = 1.0;
			fallback[2] = 0.0;
		}

		projection = Dot(fallback, axes[0]);
		for (unsigned int i = 0; i < 3; i++)
			axes[1][i] = fallback[i] - projection * axes[0][i];

		if (!Normalize(axes[1]))
			return false;
	}

	Cross(axes[0], axes[1], axes[2]);
	if (!Normalize(axes[2]))
		return false;

	Cross(axes[2], axes[0], axes[1]);
	return Normalize(axes[1]);
}

void CBeam3D2::LocalStiffness(double stiffness[12][12]) const
{
	ClearMatrix(stiffness);

	CBeamMaterial* material = dynamic_cast<CBeamMaterial*>(ElementMaterial_);
	double L = Length();
	if (!material || L < EPS)
	{
		cerr << "*** Error *** Invalid beam element material or zero length." << endl;
		exit(5);
	}

	double EA_L = material->E * material->Area / L;
	SetSymmetric(stiffness, 0, 0, EA_L);
	SetSymmetric(stiffness, 0, 6, -EA_L);
	SetSymmetric(stiffness, 6, 6, EA_L);

	double GJ_L = material->G() * material->J / L;
	SetSymmetric(stiffness, 3, 3, GJ_L);
	SetSymmetric(stiffness, 3, 9, -GJ_L);
	SetSymmetric(stiffness, 9, 9, GJ_L);

	FillBendingZ(stiffness, material->E * material->Iz, L);
	FillBendingY(stiffness, material->E * material->Iy, L);
}

void CBeam3D2::Transformation(double transform[12][12]) const
{
	ClearMatrix(transform);

	double axes[3][3];
	if (!LocalAxes(axes))
	{
		cerr << "*** Error *** Cannot define local axes for beam element." << endl;
		exit(5);
	}

	for (unsigned int block = 0; block < 4; block++)
	{
		unsigned int offset = block * 3;
		for (unsigned int i = 0; i < 3; i++)
			for (unsigned int j = 0; j < 3; j++)
				transform[offset + i][offset + j] = axes[i][j];
	}
}

void CBeam3D2::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	double local[12][12];
	double transform[12][12];
	double global[12][12];
	ClearMatrix(global);

	LocalStiffness(local);
	Transformation(transform);

	for (unsigned int i = 0; i < 12; i++)
		for (unsigned int j = 0; j < 12; j++)
			for (unsigned int a = 0; a < 12; a++)
				for (unsigned int b = 0; b < 12; b++)
					global[i][j] += transform[a][i] * local[a][b] * transform[b][j];

	PackUpperTriangle(global, Matrix);
}

void CBeam3D2::ElementEndForce(double* force, double* Displacement)
{
	double local[12][12];
	double transform[12][12];
	double globalDisplacement[12];
	double localDisplacement[12];

	LocalStiffness(local);
	Transformation(transform);

	for (unsigned int i = 0; i < 12; i++)
	{
		unsigned int equation = LocationMatrix_[i];
		globalDisplacement[i] = equation ? Displacement[equation - 1] : 0.0;
		localDisplacement[i] = 0.0;
		force[i] = 0.0;
	}

	for (unsigned int i = 0; i < 12; i++)
		for (unsigned int j = 0; j < 12; j++)
			localDisplacement[i] += transform[i][j] * globalDisplacement[j];

	for (unsigned int i = 0; i < 12; i++)
		for (unsigned int j = 0; j < 12; j++)
			force[i] += local[i][j] * localDisplacement[j];
}

void CBeam3D2::ElementStress(double* stress, double* Displacement)
{
	double force[12];
	ElementEndForce(force, Displacement);

	stress[0] = -force[0];	// positive tension
	stress[1] = force[4];	// local My at node I
	stress[2] = force[5];	// local Mz at node I
	stress[3] = -force[3];	// positive torsion
	stress[4] = force[1];	// local shear Vy at node I
	stress[5] = force[2];	// local shear Vz at node I
}

void CBeam3D2::ElementBodyForce(double* bodyForce, const double gravity[3])
{
	clear(bodyForce, ND_);

	CBeamMaterial* material = dynamic_cast<CBeamMaterial*>(ElementMaterial_);
	if (!material)
		return;

	double nodalMass = material->rho * material->Area * Length() / 2.0;
	for (unsigned int node = 0; node < 2; node++)
		for (unsigned int dof = 0; dof < 3; dof++)
			bodyForce[node * 6 + dof] = nodalMass * gravity[dof];
}
