/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "Node.h"

CNode::CNode(double X, double Y, double Z)
{
    XYZ[0] = X;		// Coordinates of the node
    XYZ[1] = Y;
    XYZ[2] = Z;
    
	for (unsigned int i = 0; i < NDF; i++)
		bcode[i] = 0;	// Boundary codes
};

//	Read element data from stream Input
bool CNode::Read(ifstream& Input)
{
	Input >> NodeNumber;	// node number

	string line;
	getline(Input, line);

	istringstream stream(line);
	vector<double> values;
	double value;
	while (stream >> value)
		values.push_back(value);

	if (values.size() == 9)
	{
		for (unsigned int i = 0; i < NDF; i++)
			bcode[i] = static_cast<unsigned int>(values[i]);

		XYZ[0] = values[6];
		XYZ[1] = values[7];
		XYZ[2] = values[8];
	}
	else if (values.size() == 6)
	{
		// Backward compatibility for legacy STAPpp 3-DOF node lines.
		for (unsigned int i = 0; i < 3; i++)
			bcode[i] = static_cast<unsigned int>(values[i]);

		for (unsigned int i = 3; i < NDF; i++)
			bcode[i] = 1;

		XYZ[0] = values[3];
		XYZ[1] = values[4];
		XYZ[2] = values[5];
	}
	else
	{
		cerr << "*** Error *** Invalid node line for node " << NodeNumber
			 << ". Expected 6 or 9 values after node number." << endl;
		return false;
	}

	return true;
}

//	Output nodal point data to stream
void CNode::Write(COutputter& output)
{
	output << setw(9) << NodeNumber;
	for (unsigned int i = 0; i < NDF; i++)
		output << setw(5) << bcode[i];

	output << setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
}

//	Output equation numbers of nodal point to stream
void CNode::WriteEquationNo(COutputter& output)
{
	output << setw(9) << NodeNumber << "       ";

	for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
	{
		output << setw(5) << bcode[dof];
	}

	output << endl;
}

//	Write nodal displacement
void CNode::WriteNodalDisplacement(COutputter& output, double* Displacement)
{
	output << setw(5) << NodeNumber << "        ";

	for (unsigned int j = 0; j < NDF; j++)
	{
		if (bcode[j] == 0)
		{
			output << setw(18) << 0.0;
		}
		else
		{
			output << setw(18) << Displacement[bcode[j] - 1];
		}
	}

	output << endl;
}
