/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Node.h"
#include "ElementGroup.h"
#include "Outputter.h"
#include "Solver.h"
#include "LoadCaseData.h"
#include "SkylineMatrix.h"
#include "MpcConstraint.h"

#ifdef STAPPP_USE_EIGEN
#include <Eigen/Sparse>
#endif

#include <vector>

using namespace std;

//!	Clear an array
template <class type> void clear( type* a, unsigned int N );

//!	Domain class : Define the problem domain
/*!	Only a single instance of Domain class can be created */
class CDomain
{
private:

//!	The instance of the Domain class
	static CDomain* _instance;

//!	Input file stream for reading data from input data file
	ifstream Input;

//!	Heading information for use in labeling the outpu
	char Title[256]; 

//!	Solution MODEX
/*!		0 : Data check only;
		1 : Execution */
	unsigned int MODEX;

//!	Total number of nodal points
	unsigned int NUMNP;

//!	List of all nodes in the domain
	CNode* NodeList;

//!	Total number of element groups.
/*! An element group consists of a convenient collection of elements with same type */
	unsigned int NUMEG;

//! Element group list
    CElementGroup* EleGrpList;
    
//!	Number of load cases
	unsigned int NLCASE;

//!	List of all load cases
	CLoadCaseData* LoadCases;

//!	Number of concentrated loads applied in each load case
	unsigned int* NLOAD;

//!	Total number of equations in the system
	unsigned int NEQ;

//!	Banded stiffness matrix
/*! A one-dimensional array storing only the elements below the	skyline of the 
    global stiffness matrix. */
    CSkylineMatrix<double>* StiffnessMatrix;

//!	Global nodal force/displacement vector
	double* Force;

//!	Penalty factor for displacement boundary conditions
	double DisplacementPenalty_;

//!	Penalty factor for multi-point constraints
	double MpcPenalty_;

//!	Multi-point constraints
	vector<CMpcConstraint> MpcConstraints_;

//!	Element group header already read when parsing legacy input
	bool HasBufferedElementHeader_;
	ElementTypes BufferedElementType_;
	unsigned int BufferedNUME_;
	unsigned int BufferedNUMMAT_;

private:

//!	Constructor
	CDomain();

//!	Desconstructor
	~CDomain();

public:

//!	Return pointer to the instance of the Domain class
	static CDomain* GetInstance();

//!	Read domain data from the input data file
	bool ReadData(string FileName, string OutFile);

//!	Read nodal point data
	bool ReadNodalPoints();

//!	Read load case data
	bool ReadLoadCases();

//!	Read optional multi-point constraints
	bool ReadMultiPointConstraints();

//!	Read element data
	bool ReadElements();

//!	Calculate global equation numbers corresponding to every degree of freedom of each node
	void CalculateEquationNumber();

//!	Calculate column heights
	void CalculateColumnHeights();

//! Allocate storage for matrices
/*!	Allocate Force, ColumnHeights, DiagonalAddress and StiffnessMatrix and 
    calculate the column heights and address of diagonal elements */
	void AllocateMatrices();

#ifdef STAPPP_USE_EIGEN
//! Allocate only global force/displacement vector for direct sparse assembly
	void AllocateForceVector();

//! Assemble directly into an Eigen sparse matrix, bypassing skyline storage
	void AssembleSparseStiffnessMatrix(Eigen::SparseMatrix<double>& SparseMatrix);

//! Apply penalty terms for displacement boundary conditions to a sparse matrix
	void ApplyDisplacementPenalty(Eigen::SparseMatrix<double>& SparseMatrix);
#endif

//!	Assemble the banded gloabl stiffness matrix
	void AssembleStiffnessMatrix();

//!	Assemble the global nodal force vector for load case LoadCase
	bool AssembleForce(unsigned int LoadCase); 

//!	Apply penalty terms for displacement boundary conditions
	void ApplyDisplacementPenalty();

//!	Apply penalty terms for multi-point constraints
	void ApplyMultiPointConstraintPenalty();

//!	Return solution mode
	inline unsigned int GetMODEX() { return MODEX; }

//!	Return the title of problem
	inline string GetTitle() { return Title; }

//!	Return the total number of equations
	inline unsigned int GetNEQ() { return NEQ; }

//!	Return the total number of nodal points
	inline unsigned int GetNUMNP() { return NUMNP; }

//!	Return the node list
	inline CNode* GetNodeList() { return NodeList; }

//!	Return total number of element groups
	inline unsigned int GetNUMEG() { return NUMEG; }

//! Return element group list
    inline CElementGroup* GetEleGrpList() { return EleGrpList; }

//!	Return pointer to the global nodal force vector
	inline double* GetForce() { return Force; }

//!	Return pointer to the global nodal displacement vector
	inline double* GetDisplacement() { return Force; }

//!	Return the total number of load cases
	inline unsigned int GetNLCASE() { return NLCASE; }

//!	Return the number of concentrated loads applied in each load case
	inline unsigned int* GetNLOAD() { return NLOAD; }

//!	Return the list of load cases
	inline CLoadCaseData* GetLoadCases() { return LoadCases; }

//!	Return pointer to the banded stiffness matrix
	inline CSkylineMatrix<double>* GetStiffnessMatrix() { return StiffnessMatrix; }

//!	Return multi-point constraints
	inline const vector<CMpcConstraint>& GetMpcConstraints() const { return MpcConstraints_; }

};
