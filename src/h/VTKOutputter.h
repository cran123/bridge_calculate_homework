/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*****************************************************************************/

#pragma once

#include <string>

class CVTKOutputter
{
public:
	static bool WriteVTU(const std::string& FileName, unsigned int LoadCase);
};
