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

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:

	unsigned int nset;	//!< Number of set
	
	double E;  //!< Young's modulus

public:

//! Virtual deconstructor
    virtual ~CMaterial() {};

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input) = 0;

//!	Write material data to Stream
    virtual void Write(COutputter& output) = 0;

};

//!	Material class for bar element
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element
	double rho;		//!< Mass density

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};

//!	Material and section class for a 3D beam element
class CBeamMaterial : public CMaterial
{
public:
	double nu;		//!< Poisson's ratio
	double Area;	//!< Sectional area
	double Iy;		//!< Second moment of area about local y axis
	double Iz;		//!< Second moment of area about local z axis
	double J;		//!< Torsion constant
	double rho;		//!< Mass density

public:

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);

//!	Return shear modulus
	double G() const { return E / (2.0 * (1.0 + nu)); }
};

//! Material class for Hex8 element
class CHex8Material : public CMaterial
{
public:
	double Nu;		//!< Poisson's ratio
	double Density;	//!< Material density for self-weight

public:
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};
