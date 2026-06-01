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

//! Hexahedral 8-node solid element
class CHex8 : public CElement
{
public:

//!	Constructor
	CHex8();

//!	Desconstructor
	~CHex8();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress (Sxx, Syy, Szz, Sxy, Syz, Sxz)
	virtual void ElementStress(double* stress, double* Displacement);

//!	Calculate element body force vector (self-weight)
	virtual void ElementBodyForce(double* bodyForce, const double gravity[3]);
};
