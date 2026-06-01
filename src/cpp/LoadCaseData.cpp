/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "LoadCaseData.h"

#include <iomanip>
#include <iostream>
#include <sstream>

using namespace std;

CLoadCaseData :: ~CLoadCaseData()
{
	delete [] node;
	delete [] dof;
	delete [] load;
	delete [] dispNode;
	delete [] dispDof;
	delete [] dispValue;
}

void CLoadCaseData :: Allocate(unsigned int num)
{
	nloads = num;
	node = new unsigned int[nloads];
	dof = new unsigned int[nloads];
	load = new double[nloads];
}; 

void CLoadCaseData :: AllocateDisp(unsigned int num)
{
	ndisp = num;
	dispNode = new unsigned int[ndisp];
	dispDof = new unsigned int[ndisp];
	dispValue = new double[ndisp];
};

//	Read load case data from stream Input
bool CLoadCaseData :: Read(ifstream& Input)
{
//	Load case number (LL) and number of concentrated loads in this load case(NL)
	string line;
	if (!std::getline(Input >> std::ws, line))
		return false;

	std::istringstream iss(line);
	unsigned int NL = 0;
	if (!(iss >> NL))
		return false;

	gravityFlag = 1;
	gravity[0] = 0.0;
	gravity[1] = 0.0;
	gravity[2] = -9.81;

	ndisp = 0;

	unsigned int flag = 0;
	double gx = 0.0, gy = 0.0, gz = 0.0;
	if (iss >> flag >> gx >> gy >> gz)
	{
		gravityFlag = flag;
		gravity[0] = gx;
		gravity[1] = gy;
		gravity[2] = gz;

		unsigned int nd = 0;
		if (iss >> nd)
			ndisp = nd;
	}

	Allocate(NL);

	for (unsigned int i = 0; i < NL; i++)
		Input >> node[i] >> dof[i] >> load[i];

	if (ndisp > 0)
	{
		AllocateDisp(ndisp);
		for (unsigned int i = 0; i < ndisp; i++)
			Input >> dispNode[i] >> dispDof[i] >> dispValue[i];
	}

	return true;
}

//	Write load case data to stream
void CLoadCaseData::Write(COutputter& output)
{
	for (unsigned int i = 0; i < nloads; i++)
		output << setw(7) << node[i] << setw(13) << dof[i]  << setw(19) << load[i] << endl;

	if (ndisp > 0)
	{
		output << endl;
		output << "    NODE       DOF      DISPLACEMENT" << endl;
		for (unsigned int i = 0; i < ndisp; i++)
			output << setw(7) << dispNode[i] << setw(10) << dispDof[i] << setw(19) << dispValue[i] << endl;
	}
}
