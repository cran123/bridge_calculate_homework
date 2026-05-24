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

#include <iomanip>

using namespace std;

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
}

//	Calculate element stress
void CPlate4::ElementStress(double* stress, double* Displacement)
{
	clear(stress, 7);
}
