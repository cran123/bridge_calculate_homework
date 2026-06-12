/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*****************************************************************************/

#include "VTKOutputter.h"

#include "Beam3D2.h"
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
			case ElementTypes::Beam:
				return 3;	// VTK_LINE
			case ElementTypes::Plate:
				return 9;	// VTK_QUAD
			case ElementTypes::H8:
			case ElementTypes::BbarH8:
				return 12;	// VTK_HEXAHEDRON
			default:
				return 0;
		}
	}

	void WriteCellField(ofstream& output, const char* name, unsigned int totalCells, const vector<double>& values)
	{
		output << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
		for (unsigned int i = 0; i < totalCells; i++)
			output << "          " << values[i] << "\n";
		output << "        </DataArray>\n";
	}
}

bool CVTKOutputter::WriteVTU(const string& FileName, unsigned int LoadCase)
{
	(void)LoadCase;

	CDomain* FEMData = CDomain::GetInstance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();

	ofstream output(FileName.c_str());
	if (!output)
	{
		cerr << "*** Error *** Cannot open VTU output file " << FileName << endl;
		return false;
	}

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int totalCells = 0;

	for (unsigned int eg = 0; eg < NUMEG; eg++)
		totalCells += FEMData->GetEleGrpList()[eg].GetNUME();

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

	vector<double> vonMises(totalCells, 0.0);
	vector<double> sxx(totalCells, 0.0);
	vector<double> syy(totalCells, 0.0);
	vector<double> szz(totalCells, 0.0);
	vector<double> sxy(totalCells, 0.0);
	vector<double> syz(totalCells, 0.0);
	vector<double> sxz(totalCells, 0.0);
	vector<double> mx(totalCells, 0.0);
	vector<double> my(totalCells, 0.0);
	vector<double> mxy(totalCells, 0.0);
	vector<double> axialForce(totalCells, 0.0);
	vector<double> beamMy(totalCells, 0.0);
	vector<double> beamMz(totalCells, 0.0);
	vector<double> beamTorque(totalCells, 0.0);

	unsigned int cell = 0;
	for (unsigned int eg = 0; eg < NUMEG; eg++)
	{
		CElementGroup& group = FEMData->GetEleGrpList()[eg];
		for (unsigned int e = 0; e < group.GetNUME(); e++)
		{
			if (group.GetElementType() == ElementTypes::Bar)
			{
				double stress = 0.0;
				group[e].ElementStress(&stress, Displacement);
				CBarMaterial& material = *dynamic_cast<CBarMaterial*>(group[e].GetElementMaterial());
				sxx[cell] = stress;
				axialForce[cell] = stress * material.Area;
			}
			else if (group.GetElementType() == ElementTypes::Beam)
			{
				double stress[6] = {};
				group[e].ElementStress(stress, Displacement);
				axialForce[cell] = stress[0];
				beamMy[cell] = stress[1];
				beamMz[cell] = stress[2];
				beamTorque[cell] = stress[3];
				sxy[cell] = stress[4];
				sxz[cell] = stress[5];
			}
			else if (group.GetElementType() == ElementTypes::H8 || group.GetElementType() == ElementTypes::BbarH8)
			{
				double stress[6] = {};
				group[e].ElementStress(stress, Displacement);
				sxx[cell] = stress[0];
				syy[cell] = stress[1];
				szz[cell] = stress[2];
				sxy[cell] = stress[3];
				syz[cell] = stress[4];
				sxz[cell] = stress[5];
				vonMises[cell] = sqrt(0.5 * ((sxx[cell] - syy[cell]) * (sxx[cell] - syy[cell])
					+ (syy[cell] - szz[cell]) * (syy[cell] - szz[cell])
					+ (szz[cell] - sxx[cell]) * (szz[cell] - sxx[cell]))
					+ 3.0 * (sxy[cell] * sxy[cell] + syz[cell] * syz[cell] + sxz[cell] * sxz[cell]));
			}
			else if (group.GetElementType() == ElementTypes::Plate)
			{
				double stress[7] = {};
				group[e].ElementStress(stress, Displacement);
				sxx[cell] = stress[0];
				syy[cell] = stress[1];
				sxy[cell] = stress[2];
				mx[cell] = stress[3];
				my[cell] = stress[4];
				mxy[cell] = stress[5];
				vonMises[cell] = sqrt(sxx[cell] * sxx[cell] - sxx[cell] * syy[cell]
					+ syy[cell] * syy[cell] + 3.0 * sxy[cell] * sxy[cell]);
			}
			cell++;
		}
	}

	WriteCellField(output, "VonMises", totalCells, vonMises);
	WriteCellField(output, "SXX", totalCells, sxx);
	WriteCellField(output, "SYY", totalCells, syy);
	WriteCellField(output, "SZZ", totalCells, szz);
	WriteCellField(output, "SXY", totalCells, sxy);
	WriteCellField(output, "SYZ", totalCells, syz);
	WriteCellField(output, "SXZ", totalCells, sxz);
	WriteCellField(output, "MX", totalCells, mx);
	WriteCellField(output, "MY", totalCells, my);
	WriteCellField(output, "MXY", totalCells, mxy);
	WriteCellField(output, "AxialForce", totalCells, axialForce);
	WriteCellField(output, "BeamMy", totalCells, beamMy);
	WriteCellField(output, "BeamMz", totalCells, beamMz);
	WriteCellField(output, "BeamTorque", totalCells, beamTorque);

	output << "      </CellData>\n";
	output << "    </Piece>\n";
	output << "  </UnstructuredGrid>\n";
	output << "</VTKFile>\n";

	return true;
}
