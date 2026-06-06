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

#include "Element.h"

using namespace std;

//! Material class for four-node engineering S4R plate/shell element
class CPlateMaterial : public CMaterial
{
public:

	double Nu;          //!< Poisson's ratio
	double Thickness;   //!< Plate thickness
	double Density;     //!< Mass density

public:

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to stream
	virtual void Write(COutputter& output);
};

//! Four-node engineering S4R plate/shell element for XY-plane plates
class CPlate4 : public CElement
{
public:

//!	Constructor
	CPlate4();

//!	Destructor
	~CPlate4();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//! Generate location matrix
	virtual void GenerateLocationMatrix();

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

//!	Calculate equivalent nodal body force
	virtual void ElementBodyForce(double* bodyForce, const double gravity[3]);
};
