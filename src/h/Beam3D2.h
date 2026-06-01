/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! Two-node 3D Euler-Bernoulli beam element, compatible with ABAQUS B31
class CBeam3D2 : public CElement
{
private:
	double Reference_[3];	//!< Reference vector defining the local y direction

	double Length() const;
	bool LocalAxes(double axes[3][3]) const;
	void LocalStiffness(double stiffness[12][12]) const;
	void Transformation(double transform[12][12]) const;

public:

//!	Constructor
	CBeam3D2();

//!	Desconstructor
	~CBeam3D2();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//! Generate location matrix
	virtual void GenerateLocationMatrix();

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress and section forces
	virtual void ElementStress(double* stress, double* Displacement);

//!	Calculate equivalent nodal body force
	virtual void ElementBodyForce(double* bodyForce, const double gravity[3]);

//!	Calculate local element end force vector
	void ElementEndForce(double* force, double* Displacement);
};
