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

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

//! Outputer class is used to output results
class COutputter
{
private:

//!	File stream for output
	ofstream OutputFile;

//! Mirror report output to console
    bool ConsoleOutput_;

//! Write detailed model/result tables
    bool DetailedOutput_;

//! Write input/model echo tables
    bool ModelOutput_;

//! Write result tables
    bool ResultOutput_;

//! Write VTK result files
    bool VtkOutput_;

//!	Designed as a single instance class
	static COutputter* _instance;

//! Constructor
    COutputter(string FileName);

public:

//!	Return pointer to the output file stream
	inline ofstream* GetOutputFile() { return &OutputFile; }

//!	Return whether detailed output tables are enabled
    inline bool DetailedOutputEnabled() const { return DetailedOutput_; }

//!	Return whether input/model echo tables are enabled
    inline bool ModelOutputEnabled() const { return ModelOutput_; }

//!	Return whether result tables are enabled
    inline bool ResultOutputEnabled() const { return ResultOutput_; }

//!	Return whether VTK result files are enabled
    inline bool VtkOutputEnabled() const { return VtkOutput_; }

//!	Return the single instance of the class
	static COutputter* GetInstance(string FileName = " ");

//!	Output current time and date
	void PrintTime(const struct tm * ptm, COutputter& output);

//!	Output logo and heading 
	void OutputHeading();

//!	Output nodal point data
	void OutputNodeInfo();

//!	Output equation numbers
	void OutputEquationNumber();

//!	Output element data
	void OutputElementInfo();

//!	Output bar element data
	void OutputBarElements(unsigned int EleGrp);

//!	Output beam element data
	void OutputBeamElements(unsigned int EleGrp);

//!	Output hex8 element data
	void OutputHex8Elements(unsigned int EleGrp);

//!	Output plate element data
	void OutputPlateElements(unsigned int EleGrp);

//!	Output load data 
	void OutputLoadInfo(); 

//!	Output displacement data
	void OutputNodalDisplacement();

//!	Output element stresses 
	void OutputElementStress();

//!	Print total system data
	void OutputTotalSystemData();

//! Overload the operator <<
	template <typename T>
	COutputter& operator<<(const T& item) 
	{
		if (ConsoleOutput_)
			std::cout << item;
		OutputFile << item;
		return *this;
	}

	typedef std::basic_ostream<char, std::char_traits<char> > CharOstream;
	COutputter& operator<<(CharOstream& (*op)(CharOstream&)) 
	{
		if (ConsoleOutput_)
			op(std::cout);
		op(OutputFile);
		return *this;
	}

#ifdef _DEBUG_

//!	Print banded and full stiffness matrix for debuging
	void PrintStiffnessMatrix();

//!	Print address of diagonal elements for debuging
	void PrintDiagonalAddress();

//!	Print column heights for debuging
	void PrintColumnHeights();

//!	Print displacement vector for debuging
	void PrintDisplacement();

#endif

};
