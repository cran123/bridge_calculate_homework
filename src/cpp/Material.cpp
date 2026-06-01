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

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	string line;
	if (!std::getline(Input >> std::ws, line))
		return false;

	std::istringstream iss(line);
	if (!(iss >> E >> Area))
		return false;

	if (!(iss >> Density))
		Density = 0.0;

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << setw(16) << Density << endl;
}

//	Read material data from stream Input
bool CHex8Material::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	string line;
	if (!std::getline(Input >> std::ws, line))
		return false;

	std::istringstream iss(line);
	if (!(iss >> E >> Nu))
		return false;

	if (!(iss >> Density))
		Density = 0.0;

	return true;
}

//	Write material data to Stream
void CHex8Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Nu << setw(16) << Density << endl;
}
