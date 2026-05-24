/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*****************************************************************************/

#include "VTKOutputter.h"

#include "Domain.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

namespace
{
	double DofValue(CNode& node, unsigned int dof, double* displacement)
	{
		unsigned int equation = node.bcode[dof];
		if (!equation)
			return 0.0;

		return displacement[equation - 1];
	}

	int VTKCellType(ElementTypes type)
	{
		switch (type)
		{
			case ElementTypes::Bar:
				return 3;
			case ElementTypes::Plate:
				return 9;
			default:
				return 0;
		}
	}
}

bool CVTKOutputter::WriteVTU(const string& FileName, unsigned int LoadCase)
{
	CDomain* FEMData = CDomain::GetInstance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();

	ofstream output(FileName);
	if (!output)
	{
		cerr << "*** Error *** Cannot open VTU output file " << FileName << endl;
		return false;
	}

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int totalCells = 0;
	unsigned int totalConnectivity = 0;

	for (unsigned int eg = 0; eg < NUMEG; eg++)
	{
		CElementGroup& group = FEMData->GetEleGrpList()[eg];
		totalCells += group.GetNUME();
		for (unsigned int e = 0; e < group.GetNUME(); e++)
			totalConnectivity += group[e].GetNEN();
	}

	output << "<?xml version=\"1.0\"?>\n";
	output << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	output << "  <UnstructuredGrid>\n";
	output << "    <Piece NumberOfPoints=\"" << NUMNP << "\" NumberOfCells=\"" << totalCells << "\">\n";

	output << "      <Points>\n";
	output << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	output << setiosflags(ios::scientific) << setprecision(12);
	for (unsigned int i = 0; i < NUMNP; i++)
		output << "          " << NodeList[i].XYZ[0] << " " << NodeList[i].XYZ[1] << " " << NodeList[i].XYZ[2] << "\n";
	output << "        </DataArray>\n";
	output << "      </Points>\n";

	output << "      <Cells>\n";
	output << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
	for (unsigned int eg = 0; eg < NUMEG; eg++)
	{
		CElementGroup& group = FEMData->GetEleGrpList()[eg];
		for (unsigned int e = 0; e < group.GetNUME(); e++)
		{
			CNode** nodes = group[e].GetNodes();
			output << "          ";
			for (unsigned int n = 0; n < group[e].GetNEN(); n++)
				output << nodes[n]->NodeNumber - 1 << " ";
			output << "\n";
		}
	}
	output << "        </DataArray>\n";

	output << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	unsigned int offset = 0;
	for (unsigned int eg = 0; eg < NUMEG; eg++)
	{
		CElementGroup& group = FEMData->GetEleGrpList()[eg];
		for (unsigned int e = 0; e < group.GetNUME(); e++)
		{
			offset += group[e].GetNEN();
			output << "          " << offset << "\n";
		}
	}
	output << "        </DataArray>\n";

	output << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	for (unsigned int eg = 0; eg < NUMEG; eg++)
	{
		CElementGroup& group = FEMData->GetEleGrpList()[eg];
		int cellType = VTKCellType(group.GetElementType());
		for (unsigned int e = 0; e < group.GetNUME(); e++)
			output << "          " << cellType << "\n";
	}
	output << "        </DataArray>\n";
	output << "      </Cells>\n";

	output << "      <PointData>\n";
	output << "        <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < NUMNP; i++)
		output << "          " << DofValue(NodeList[i], 0, Displacement) << " "
			   << DofValue(NodeList[i], 1, Displacement) << " "
			   << DofValue(NodeList[i], 2, Displacement) << "\n";
	output << "        </DataArray>\n";
	output << "        <DataArray type=\"Float64\" Name=\"Rotation\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < NUMNP; i++)
		output << "          " << DofValue(NodeList[i], 3, Displacement) << " "
			   << DofValue(NodeList[i], 4, Displacement) << " "
			   << DofValue(NodeList[i], 5, Displacement) << "\n";
	output << "        </DataArray>\n";
	output << "        <DataArray type=\"Float64\" Name=\"DisplacementMagnitude\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < NUMNP; i++)
	{
		double ux = DofValue(NodeList[i], 0, Displacement);
		double uy = DofValue(NodeList[i], 1, Displacement);
		double uz = DofValue(NodeList[i], 2, Displacement);
		output << "          " << sqrt(ux * ux + uy * uy + uz * uz) << "\n";
	}
	output << "        </DataArray>\n";
	output << "      </PointData>\n";

	output << "      <CellData>\n";
	output << "        <DataArray type=\"Int32\" Name=\"ElementType\" format=\"ascii\">\n";
	for (unsigned int eg = 0; eg < NUMEG; eg++)
	{
		CElementGroup& group = FEMData->GetEleGrpList()[eg];
		for (unsigned int e = 0; e < group.GetNUME(); e++)
			output << "          " << static_cast<int>(group.GetElementType()) << "\n";
	}
	output << "        </DataArray>\n";

	const char* names[7] = {"SXX", "SYY", "SXY", "MX", "MY", "MXY", "VonMises"};
	for (unsigned int field = 0; field < 7; field++)
	{
		output << "        <DataArray type=\"Float64\" Name=\"" << names[field] << "\" format=\"ascii\">\n";
		for (unsigned int eg = 0; eg < NUMEG; eg++)
		{
			CElementGroup& group = FEMData->GetEleGrpList()[eg];
			for (unsigned int e = 0; e < group.GetNUME(); e++)
			{
				double stress[7] = {};
				if (group.GetElementType() == ElementTypes::Plate)
					group[e].ElementStress(stress, Displacement);
				output << "          " << stress[field] << "\n";
			}
		}
		output << "        </DataArray>\n";
	}
	output << "      </CellData>\n";

	output << "    </Piece>\n";
	output << "  </UnstructuredGrid>\n";
	output << "</VTKFile>\n";

	return true;
}
