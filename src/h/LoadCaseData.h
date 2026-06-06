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

#include "Outputter.h"

using namespace std;

//! Class LoadData is used to store load data
class CLoadCaseData
{
public:

	unsigned int nloads;	//!< Number of concentrated loads in this load case
	unsigned int* node;		//!< Node number to which this load is applied
	unsigned int* dof;		//!< Degree of freedom number for this load component
	double* load;			//!< Magnitude of load
	bool hasGravity;		//!< Whether gravity load is applied
	double gravity[3];		//!< Gravity acceleration vector
	unsigned int ndisp;		//!< Number of prescribed displacements
	unsigned int* dispNode;		//!< Node number for displacement
	unsigned int* dispDof;		//!< Degree of freedom number for displacement
	double* dispValue;		//!< Prescribed displacement value

public:

	CLoadCaseData() : nloads(0), node(NULL), dof(NULL), load(NULL), hasGravity(false), gravity{0.0, 0.0, 0.0},
		ndisp(0), dispNode(NULL), dispDof(NULL), dispValue(NULL) {};
	~CLoadCaseData();

//!	Set nloads, and new array node, dof and load
	void Allocate(unsigned int num);

//!	Set ndisp, and new array dispNode, dispDof and dispValue
	void AllocateDisp(unsigned int num);

//!	Read load case data from stream Input
	bool Read(ifstream& Input);

//!	Write load case data to stream
	void Write(COutputter& output);
};
