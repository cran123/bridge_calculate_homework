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
#include <string>

using namespace std;

CLoadCaseData :: ~CLoadCaseData()
{
	delete [] node;
	delete [] dof;
	delete [] load;
}

void CLoadCaseData :: Allocate(unsigned int num)
{
	nloads = num;
	node = new unsigned int[nloads];
	dof = new unsigned int[nloads];
	load = new double[nloads];
}; 

//	Read load case data from stream Input
bool CLoadCaseData :: Read(ifstream& Input)
{
//	Load case number (LL) and number of concentrated loads in this load case(NL)
	
	unsigned int NL;
	string line;

	Input >> NL;
	getline(Input, line);

	istringstream lineStream(line);
	unsigned int gravityFlag = 0;
	double gx = 0.0;
	double gy = 0.0;
	double gz = 0.0;
	if (lineStream >> gravityFlag >> gx >> gy >> gz)
	{
		hasGravity = gravityFlag != 0;
		gravity[0] = gx;
		gravity[1] = gy;
		gravity[2] = gz;
	}
	else
	{
		hasGravity = false;
		gravity[0] = gravity[1] = gravity[2] = 0.0;
	}

	Allocate(NL);

	for (unsigned int i = 0; i < NL; i++)
		Input >> node[i] >> dof[i] >> load[i];

	return true;
}

//	Write load case data to stream
void CLoadCaseData::Write(COutputter& output)
{
	if (hasGravity)
	{
		output << " GRAVITY LOAD VECTOR:"
			   << setw(16) << gravity[0]
			   << setw(16) << gravity[1]
			   << setw(16) << gravity[2] << endl;
	}

	for (unsigned int i = 0; i < nloads; i++)
		output << setw(7) << node[i] << setw(13) << dof[i]  << setw(19) << load[i] << endl;
}
