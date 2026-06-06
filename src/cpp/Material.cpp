/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	string line;
	getline(Input, line);
	istringstream stream(line);

	vector<double> values;
	double value;
	while (stream >> value)
		values.push_back(value);

	if (values.size() < 2)
	{
		cerr << "*** Error *** Invalid truss material set " << nset << endl;
		return false;
	}

	E = values[0];
	Area = values[1];
	rho = values.size() > 2 ? values[2] : 0.0;

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << setw(16) << rho << endl;
}

//	Read beam material data from stream Input
bool CBeamMaterial::Read(ifstream& Input)
{
	Input >> nset;

	string line;
	getline(Input, line);
	istringstream stream(line);

	vector<double> values;
	double value;
	while (stream >> value)
		values.push_back(value);

	if (values.size() < 7)
	{
		cerr << "*** Error *** Invalid beam material set " << nset
			 << ". Expected: set E nu A Iy Iz J rho" << endl;
		return false;
	}

	E = values[0];
	nu = values[1];
	Area = values[2];
	Iy = values[3];
	Iz = values[4];
	J = values[5];
	rho = values[6];

	return true;
}

//	Write beam material data to Stream
void CBeamMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << setw(16) << Area
		   << setw(16) << Iy << setw(16) << Iz << setw(16) << J
		   << setw(16) << rho << endl;
}

//	Read Hex8 material data from stream Input
bool CHex8Material::Read(ifstream& Input)
{
	Input >> nset;

	string line;
	getline(Input, line);
	istringstream stream(line);

	vector<double> values;
	double value;
	while (stream >> value)
		values.push_back(value);

	if (values.size() < 2)
	{
		cerr << "*** Error *** Invalid Hex8 material set " << nset << endl;
		return false;
	}

	E = values[0];
	Nu = values[1];
	Density = values.size() > 2 ? values[2] : 0.0;

	return true;
}

//	Write Hex8 material data to Stream
void CHex8Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Nu << setw(16) << Density << endl;
}
