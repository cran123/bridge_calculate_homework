/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*                                                                           */
/*     B-bar Hex8 element (B-bar method for reduced volumetric locking)      */
/*                                                                           */
/*****************************************************************************/

#pragma once

#include "Element.h"
#include "Material.h"
#include "Node.h"

//! B-bar Hex8 element class
class CBbarHex8 : public CElement
{
public:
    CBbarHex8();
    ~CBbarHex8();

    bool Read(std::ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);
    void Write(COutputter& output);
    void GenerateLocationMatrix();
    void ElementStiffness(double* Matrix);
    void ElementStress(double* stress, double* Displacement);
    void ElementBodyForce(double* bodyForce, const double gravity[3]);
};
