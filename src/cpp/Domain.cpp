/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Material.h"

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>

using namespace std;

//	Clear an array
template <class type> void clear( type* a, unsigned int N )
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

template void clear<double>(double* a, unsigned int N);

CDomain* CDomain::_instance = nullptr;

//	Constructor
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0;

	NUMNP = 0;
	NodeList = nullptr;
	
	NUMEG = 0;
	EleGrpList = nullptr;
	
	NLCASE = 0;
	NLOAD = nullptr;
	LoadCases = nullptr;
	
	NEQ = 0;

	Force = nullptr;
	StiffnessMatrix = nullptr;
	DisplacementPenalty_ = 0.0;
	MpcPenalty_ = 0.0;
	HasBufferedElementHeader_ = false;
	BufferedElementType_ = ElementTypes::UNDEFINED;
	BufferedNUME_ = 0;
	BufferedNUMMAT_ = 0;
}

//	Desconstructor
CDomain::~CDomain()
{
	delete [] NodeList;

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
	delete StiffnessMatrix;
}

//	Return pointer to the instance of the Domain class
CDomain* CDomain::GetInstance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile)
{
	Input.open(FileName);

	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	COutputter* Output = COutputter::GetInstance(OutFile);

//	Read the heading line
	Input.getline(Title, 256);
	Output->OutputHeading();

//	Read the control line
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	Read nodal point data
	if (ReadNodalPoints())
        Output->OutputNodeInfo();
    else
        return false;

//	Update equation number
	CalculateEquationNumber();
	Output->OutputEquationNumber();

//	Read load data
	if (ReadLoadCases())
        Output->OutputLoadInfo();
    else
        return false;

//	Read optional multi-point constraints
	if (!ReadMultiPointConstraints())
		return false;

//	Read element data
	if (ReadElements())
        Output->OutputElementInfo();
    else
        return false;

	return true;
}

namespace
{
	bool ReadNextDataLine(ifstream& input, string& line)
	{
		while (getline(input, line))
		{
			size_t first = line.find_first_not_of(" \t\r\n");
			if (first == string::npos)
				continue;
			if (line[first] == '#')
				continue;
			return true;
		}
		return false;
	}
}

//	Read nodal point data
bool CDomain::ReadNodalPoints()
{

//	Read nodal point data lines
	NodeList = new CNode[NUMNP];

//	Loop over for all nodal points
	for (unsigned int np = 0; np < NUMNP; np++)
    {
		if (!NodeList[np].Read(Input))
			return false;
    
        if (NodeList[np].NodeNumber != np + 1)
        {
            cerr << "*** Error *** Nodes must be inputted in order !" << endl
            << "   Expected node number : " << np + 1 << endl
            << "   Provided node number : " << NodeList[np].NodeNumber << endl;
        
            return false;
        }
    }

	return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber()
{
	NEQ = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode[dof]) 
				NodeList[np].bcode[dof] = 0;
			else
			{
				NEQ++;
				NodeList[np].bcode[dof] = NEQ;
			}
		}
	}
}

//	Read load case data
bool CDomain::ReadLoadCases()
{
//	Read load data lines
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases

//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
    {
        unsigned int LL;
        Input >> LL;
        
        if (LL != lcase + 1)
        {
            cerr << "*** Error *** Load case must be inputted in order !" << endl
            << "   Expected load case : " << lcase + 1 << endl
            << "   Provided load case : " << LL << endl;
            
            return false;
        }

        LoadCases[lcase].Read(Input);
    }

	return true;
}

bool CDomain::ReadMultiPointConstraints()
{
	string line;
	if (!ReadNextDataLine(Input, line))
		return true;

	istringstream header(line);
	vector<double> values;
	double value;
	while (header >> value)
		values.push_back(value);

	if (values.empty())
		return true;

	if (values.size() >= 3)
	{
		BufferedElementType_ = static_cast<ElementTypes>(static_cast<int>(values[0]));
		BufferedNUME_ = static_cast<unsigned int>(values[1]);
		BufferedNUMMAT_ = static_cast<unsigned int>(values[2]);
		HasBufferedElementHeader_ = true;
		return true;
	}

	unsigned int nmpc = static_cast<unsigned int>(values[0]);
	MpcConstraints_.clear();
	MpcConstraints_.reserve(nmpc);

	for (unsigned int i = 0; i < nmpc; i++)
	{
		unsigned int nterms;
		Input >> nterms;

		CMpcConstraint constraint;
		constraint.terms.reserve(nterms);
		constraint.equations.reserve(nterms);

		for (unsigned int j = 0; j < nterms; j++)
		{
			CMpcTerm term;
			Input >> term.node >> term.dof >> term.coefficient;
			constraint.terms.push_back(term);

			unsigned int eq = 0;
			if (term.node >= 1 && term.node <= NUMNP && term.dof >= 1 && term.dof <= CNode::NDF)
				eq = NodeList[term.node - 1].bcode[term.dof - 1];
			constraint.equations.push_back(eq);
		}

		MpcConstraints_.push_back(constraint);
	}

	return true;
}

// Read element data
bool CDomain::ReadElements()
{
    EleGrpList = new CElementGroup[NUMEG];

//	Loop over for all element group
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		if (EleGrp == 0 && HasBufferedElementHeader_)
		{
			if (!EleGrpList[EleGrp].Read(Input, BufferedElementType_, BufferedNUME_, BufferedNUMMAT_))
				return false;
		}
		else if (!EleGrpList[EleGrp].Read(Input))
		{
			return false;
		}
	}
    
    return true;
}

//	Calculate column heights
void CDomain::CalculateColumnHeights()
{
#ifdef _DEBUG_
    COutputter* Output = COutputter::GetInstance();
    *Output << setw(9) << "Ele = " << setw(22) << "Location Matrix" << endl;
#endif

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
        {
            CElement& Element = ElementGrp[Ele];

            // Generate location matrix
            Element.GenerateLocationMatrix();
            
#ifdef _DEBUG_
            unsigned int* LocationMatrix = Element.GetLocationMatrix();
            
            *Output << setw(9) << Ele+1;
            for (int i=0; i<Element.GetND(); i++)
                *Output << setw(5) << LocationMatrix[i];
            *Output << endl;
#endif

            StiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
        }
    }

    StiffnessMatrix->CalculateMaximumHalfBandwidth();
    
#ifdef _DEBUG_
    *Output << endl;
	Output->PrintColumnHeights();
#endif

}

//    Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//    and calculate the column heights and address of diagonal elements
void CDomain::AllocateMatrices()
{
    //    Allocate for global force/displacement vector
    Force = new double[NEQ];
    
    //  Create the banded stiffness matrix
    StiffnessMatrix = new CSkylineMatrix<double>(NEQ);
    
    //    Calculate column heights
    CalculateColumnHeights();
    
    //    Calculate address of diagonal elements in banded matrix
    StiffnessMatrix->CalculateDiagnoalAddress();
    
    //    Allocate for banded global stiffness matrix
    StiffnessMatrix->Allocate();
    
    COutputter* Output = COutputter::GetInstance();
    Output->OutputTotalSystemData();
}

#ifdef STAPPP_USE_EIGEN
void CDomain::AllocateForceVector()
{
	delete[] Force;
	Force = new double[NEQ];
}

void CDomain::AssembleSparseStiffnessMatrix(Eigen::SparseMatrix<double>& SparseMatrix)
{
	vector<vector<int> > rowsByColumn(NEQ);

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		CElementGroup& ElementGrp = EleGrpList[EleGrp];
		unsigned int NUME = ElementGrp.GetNUME();

		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			CElement& Element = ElementGrp[Ele];
			Element.GenerateLocationMatrix();
			unsigned int* LocationMatrix = Element.GetLocationMatrix();
			unsigned int ND = Element.GetND();

			for (unsigned int j = 0; j < ND; j++)
			{
				unsigned int Lj = LocationMatrix[j];
				if (!Lj)
					continue;
				vector<int>& rows = rowsByColumn[Lj - 1];
				for (unsigned int i = 0; i < ND; i++)
				{
					unsigned int Li = LocationMatrix[i];
					if (Li)
						rows.push_back(static_cast<int>(Li - 1));
				}
			}
		}
	}

	vector<int> reserve(NEQ, 0);
	for (unsigned int col = 0; col < NEQ; col++)
	{
		vector<int>& rows = rowsByColumn[col];
		sort(rows.begin(), rows.end());
		rows.erase(unique(rows.begin(), rows.end()), rows.end());
		reserve[col] = static_cast<int>(rows.size());
	}

	SparseMatrix.resize(NEQ, NEQ);
	SparseMatrix.reserve(reserve);
	for (unsigned int col = 0; col < NEQ; col++)
	{
		for (unsigned int k = 0; k < rowsByColumn[col].size(); k++)
			SparseMatrix.insert(rowsByColumn[col][k], col) = 0.0;
	}
	SparseMatrix.makeCompressed();

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		CElementGroup& ElementGrp = EleGrpList[EleGrp];
		unsigned int NUME = ElementGrp.GetNUME();
		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			CElement& Element = ElementGrp[Ele];
			Element.ElementStiffness(Matrix);
			unsigned int* LocationMatrix = Element.GetLocationMatrix();
			unsigned int ND = Element.GetND();

			for (unsigned int j = 0; j < ND; j++)
			{
				unsigned int Lj = LocationMatrix[j];
				if (!Lj)
					continue;

				unsigned int DiagjElement = (j + 1) * j / 2;
				for (unsigned int i = 0; i <= j; i++)
				{
					unsigned int Li = LocationMatrix[i];
					if (!Li)
						continue;

					double value = Matrix[DiagjElement + j - i];
					if (value == 0.0)
						continue;

					SparseMatrix.coeffRef(Li - 1, Lj - 1) += value;
					if (Li != Lj)
						SparseMatrix.coeffRef(Lj - 1, Li - 1) += value;
				}
			}
		}

		delete[] Matrix;
	}

	SparseMatrix.prune(0.0);
	SparseMatrix.makeCompressed();
}

void CDomain::ApplyDisplacementPenalty(Eigen::SparseMatrix<double>& SparseMatrix)
{
	double maxDiag = 0.0;
	for (unsigned int i = 0; i < NEQ; i++)
	{
		double v = fabs(SparseMatrix.coeff(i, i));
		if (v > maxDiag)
			maxDiag = v;
	}

	if (maxDiag <= 0.0)
		maxDiag = 1.0;

	DisplacementPenalty_ = maxDiag * 1.0e10;

	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
	{
		CLoadCaseData* LoadData = &LoadCases[lcase];
		for (unsigned int i = 0; i < LoadData->ndisp; i++)
		{
			unsigned int eq = NodeList[LoadData->dispNode[i] - 1].bcode[LoadData->dispDof[i] - 1];
			if (eq)
				SparseMatrix.coeffRef(eq - 1, eq - 1) += DisplacementPenalty_;
		}
	}

	SparseMatrix.makeCompressed();
}
#endif

//	Assemble the banded gloabl stiffness matrix
void CDomain::AssembleStiffnessMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
            Element.ElementStiffness(Matrix);
            StiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
        }

		delete[] Matrix;
		Matrix = nullptr;
	}

#ifdef _DEBUG_
	COutputter* Output = COutputter::GetInstance();
	Output->PrintStiffnessMatrix();
#endif

}

void CDomain::ApplyMultiPointConstraintPenalty()
{
	if (!StiffnessMatrix || MpcConstraints_.empty())
		return;

	double maxDiag = 0.0;
	for (unsigned int i = 1; i <= NEQ; i++)
	{
		double v = fabs((*StiffnessMatrix)(i, i));
		if (v > maxDiag)
			maxDiag = v;
	}

	if (maxDiag <= 0.0)
		maxDiag = 1.0;

	MpcPenalty_ = maxDiag * 1.0e8;

	for (unsigned int mpc = 0; mpc < MpcConstraints_.size(); mpc++)
	{
		CMpcConstraint& constraint = MpcConstraints_[mpc];
		for (unsigned int j = 0; j < constraint.terms.size(); j++)
		{
			unsigned int eqj = constraint.equations[j];
			if (!eqj)
				continue;

			double cj = constraint.terms[j].coefficient;
			for (unsigned int i = 0; i <= j; i++)
			{
				unsigned int eqi = constraint.equations[i];
				if (!eqi)
					continue;

				double ci = constraint.terms[i].coefficient;
				(*StiffnessMatrix)(eqi, eqj) += MpcPenalty_ * ci * cj;
			}
		}
	}
}

//	Apply penalty terms for displacement boundary conditions
void CDomain::ApplyDisplacementPenalty()
{
	if (!StiffnessMatrix)
		return;

	double maxDiag = 0.0;
	for (unsigned int i = 1; i <= NEQ; i++)
	{
		double v = fabs((*StiffnessMatrix)(i, i));
		if (v > maxDiag)
			maxDiag = v;
	}

	if (maxDiag <= 0.0)
		maxDiag = 1.0;

	DisplacementPenalty_ = maxDiag * 1.0e10;

	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
	{
		CLoadCaseData* LoadData = &LoadCases[lcase];
		for (unsigned int i = 0; i < LoadData->ndisp; i++)
		{
			unsigned int eq = NodeList[LoadData->dispNode[i] - 1].bcode[LoadData->dispDof[i] - 1];
			if (eq)
				(*StiffnessMatrix)(eq, eq) += DisplacementPenalty_;
		}
	}
}

//	Assemble the global nodal force vector for load case LoadCase
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];

    clear(Force, NEQ);

//	Loop over for all concentrated loads in load case LoadCase
	for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
        
        if(dof) // The DOF is activated
            Force[dof - 1] += LoadData->load[lnum];
	}

	if (LoadData->hasGravity)
	{
		for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
			unsigned int NUME = ElementGrp.GetNUME();
			double* BodyForce = new double[ElementGrp[0].GetND()];

			for (unsigned int Ele = 0; Ele < NUME; Ele++)
			{
				CElement& Element = ElementGrp[Ele];
				Element.ElementBodyForce(BodyForce, LoadData->gravity);
				unsigned int* LocationMatrix = Element.GetLocationMatrix();

				for (unsigned int i = 0; i < Element.GetND(); i++)
				{
					unsigned int dof = LocationMatrix[i];
					if (dof)
						Force[dof - 1] += BodyForce[i];
				}
			}

			delete[] BodyForce;
		}
	}

	if (DisplacementPenalty_ > 0.0)
	{
		for (unsigned int i = 0; i < LoadData->ndisp; i++)
		{
			unsigned int eq = NodeList[LoadData->dispNode[i] - 1].bcode[LoadData->dispDof[i] - 1];
			if (eq)
				Force[eq - 1] += DisplacementPenalty_ * LoadData->dispValue[i];
		}
	}

	return true;
}

