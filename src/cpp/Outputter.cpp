/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <ctime>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

#include "Domain.h"
#include "Outputter.h"
#include "SkylineMatrix.h"

using namespace std;

//	Output current time and date
void COutputter::PrintTime(const struct tm* ptm, COutputter &output)
{
	const char* weekday[] = {"Sunday", "Monday", "Tuesday", "Wednesday",
							 "Thursday", "Friday", "Saturday"};
	const char* month[] = {"January", "February", "March", "April", "May", "June",
						   "July", "August", "September", "October", "November", "December"};

	output << "        (";
	output << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " on ";
	output << month[ptm->tm_mon] << " " << ptm->tm_mday << ", " << ptm->tm_year + 1900 << ", "
		   << weekday[ptm->tm_wday] << ")" << endl
		   << endl;
}

COutputter* COutputter::_instance = nullptr;

//	Constructor
COutputter::COutputter(string FileName)
{
	OutputFileName_ = FileName;
	OutputFile.open(FileName);

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	Return the single instance of the class
COutputter* COutputter::GetInstance(string FileName)
{
	if (!_instance)
		_instance = new COutputter(FileName);
    
	return _instance;
}

//	Print program logo
void COutputter::OutputHeading()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "TITLE : " << FEMData->GetTitle() << endl;

	time_t rawtime;
	struct tm* timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	PrintTime(timeinfo, *this);
}

//	Print nodal data
void COutputter::OutputNodeInfo()
{
	CDomain* FEMData = CDomain::GetInstance();

	CNode* NodeList = FEMData->GetNodeList();

	*this << "C O N T R O L   I N F O R M A T I O N" << endl
		  << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int NLCASE = FEMData->GetNLCASE();
	unsigned int MODEX = FEMData->GetMODEX();

	*this << "      NUMBER OF NODAL POINTS . . . . . . . . . . (NUMNP)  =" << setw(6) << NUMNP << endl;
	*this << "      NUMBER OF ELEMENT GROUPS . . . . . . . . . (NUMEG)  =" << setw(6) << NUMEG << endl;
	*this << "      NUMBER OF LOAD CASES . . . . . . . . . . . (NLCASE) =" << setw(6) << NLCASE << endl;
	*this << "      SOLUTION MODE  . . . . . . . . . . . . . . (MODEX)  =" << setw(6) << MODEX << endl;
	*this << "         EQ.0, DATA CHECK" << endl
		  << "         EQ.1, EXECUTION" << endl
		  << endl;

	*this << " N O D A L   P O I N T   D A T A" << endl << endl;
	*this << "    NODE       BOUNDARY                         NODAL POINT" << endl
		  << "   NUMBER  CONDITION  CODES                     COORDINATES" << endl;

	for (unsigned int np = 0; np < NUMNP; np++)
		NodeList[np].Write(*this);

	*this << endl;
}

//	Output equation numbers
void COutputter::OutputEquationNumber()
{
	CDomain* FEMData = CDomain::GetInstance();
	unsigned int NUMNP = FEMData->GetNUMNP();

	CNode* NodeList = FEMData->GetNodeList();

	*this << " EQUATION NUMBERS" << endl
		  << endl;
	*this << "   NODE NUMBER   DEGREES OF FREEDOM" << endl;
	*this << "        N           X    Y    Z" << endl;

	for (unsigned int np = 0; np < NUMNP; np++) // Loop over for all node
		NodeList[np].WriteEquationNo(*this);

	*this << endl;
}

//	Output element data
void COutputter::OutputElementInfo()
{
	//	Print element group control line

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NUMEG = FEMData->GetNUMEG();

	*this << " E L E M E N T   G R O U P   D A T A" << endl
		  << endl
		  << endl;

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		*this << " E L E M E N T   D E F I N I T I O N" << endl
			  << endl;

		ElementTypes ElementType = FEMData->GetEleGrpList()[EleGrp].GetElementType();
		unsigned int NUME = FEMData->GetEleGrpList()[EleGrp].GetNUME();

		*this << " ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . =" << setw(5)
			  << ElementType << endl;
		switch (ElementType)
		{
			case ElementTypes::Bar:
				*this << "     EQ.1, TRUSS ELEMENTS" << endl;
				break;
			case ElementTypes::H8:
				*this << "     EQ.4, HEX8 SOLID ELEMENTS" << endl;
				break;
			case ElementTypes::BbarH8:
				*this << "     EQ.8, B-BAR HEX8 SOLID ELEMENTS" << endl;
				break;
			default:
				*this << "     (UNKNOWN TYPE)" << endl;
				break;
		}
		*this << endl;

		*this << " NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . =" << setw(5) << NUME
			  << endl
			  << endl;

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				OutputBarElements(EleGrp);
				break;
			case ElementTypes::H8: // Hex8 element
				OutputHex8Elements(EleGrp);
				break;
			case ElementTypes::BbarH8: // B-bar Hex8 element
				OutputHex8Elements(EleGrp);
				break;
		    default:
		        *this << ElementType << " has not been implemented yet." << endl;
		        break;
		}
	}
}

//	Output hex8 element data
void COutputter::OutputHex8Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND SOLID MATERIAL CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S          POISSON         DENSITY" << endl
		  << " NUMBER     MODULUS           RATIO" << endl
		  << "               E               NU" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset + 1;
		ElementGroup.GetMaterial(mset).Write(*this);
	}

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;

	*this << " ELEMENT     NODE     NODE     NODE     NODE     NODE     NODE     NODE     NODE   MATERIAL" << endl
		  << " NUMBER-N      1        2        3        4        5        6        7        8      SET" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele + 1;
		ElementGroup[Ele].Write(*this);
	}

	*this << endl;
}
//	Output bar element data
void COutputter::OutputBarElements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S     CROSS-SECTIONAL        DENSITY" << endl
		  << " NUMBER     MODULUS          AREA" << endl
		  << "               E              A" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
		ElementGroup.GetMaterial(mset).Write(*this);
    }

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
    
	*this << " ELEMENT     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
    }

	*this << endl;
}

//	Print load data
void COutputter::OutputLoadInfo()
{
	CDomain* FEMData = CDomain::GetInstance();

	for (unsigned int lcase = 1; lcase <= FEMData->GetNLCASE(); lcase++)
	{
		CLoadCaseData* LoadData = &FEMData->GetLoadCases()[lcase - 1];

		*this << setiosflags(ios::scientific);
		*this << " L O A D   C A S E   D A T A" << endl
			  << endl;

		*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
		*this << "     NUMBER OF CONCENTRATED LOADS . =" << setw(6) << LoadData->nloads << endl;
		*this << "     NUMBER OF DISP. CONSTRAINTS . =" << setw(6) << LoadData->ndisp << endl;
		*this << "     GRAVITY FLAG . . . . . . . . . =" << setw(6) << LoadData->gravityFlag << endl;
		*this << "     GRAVITY VECTOR . . . . . . . . =" << setw(14) << LoadData->gravity[0]
			  << setw(14) << LoadData->gravity[1] << setw(14) << LoadData->gravity[2] << endl
			  << endl;
		*this << "    NODE       DIRECTION      LOAD" << endl
			  << "   NUMBER                   MAGNITUDE" << endl;

		LoadData->Write(*this);

		*this << endl;
	}
}

//	Print nodal displacement
void COutputter::OutputNodalDisplacement()
{
	CDomain* FEMData = CDomain::GetInstance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();

	*this << setiosflags(ios::scientific);

	*this << " D I S P L A C E M E N T S" << endl
		  << endl;
	*this << "  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT" << endl;

	for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++)
		NodeList[np].WriteNodalDisplacement(*this, Displacement);

	*this << endl;
}

//	Calculate stresses
void COutputter::OutputElementStress()
{
	CDomain* FEMData = CDomain::GetInstance();

	double* Displacement = FEMData->GetDisplacement();

	unsigned int NUMEG = FEMData->GetNUMEG();

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		*this << " S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P" << setw(5)
			  << EleGrpIndex + 1 << endl
			  << endl;

		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType();

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				*this << "  ELEMENT             FORCE            STRESS" << endl
					<< "  NUMBER" << endl;

				double stress;

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.ElementStress(&stress, Displacement);

					CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());
					*this << setw(5) << Ele + 1 << setw(22) << stress * material.Area << setw(18)
						<< stress << endl;
				}

				*this << endl;

				break;

			case ElementTypes::H8: // Hex8 element
			case ElementTypes::BbarH8: // B-bar Hex8 element
				*this << "  ELEMENT        SXX          SYY          SZZ          SXY          SYZ          SXZ        VON MISES" << endl
					<< "  NUMBER" << endl;

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					double stress[6];
					Element.ElementStress(stress, Displacement);

					double sxx = stress[0];
					double syy = stress[1];
					double szz = stress[2];
					double sxy = stress[3];
					double syz = stress[4];
					double sxz = stress[5];
					double vm = sqrt(0.5 * ((sxx - syy) * (sxx - syy) + (syy - szz) * (syy - szz) + (szz - sxx) * (szz - sxx))
						+ 3.0 * (sxy * sxy + syz * syz + sxz * sxz));

					*this << setw(5) << Ele + 1
						<< setw(12) << sxx
						<< setw(12) << syy
						<< setw(12) << szz
						<< setw(12) << sxy
						<< setw(12) << syz
						<< setw(12) << sxz
						<< setw(12) << vm << endl;
				}

				*this << endl;

				break;

			default: // Invalid element type
				cerr << "*** Error *** Elment type " << ElementType
					<< " has not been implemented.\n\n";
		}
	}
}

//	Output VTK/ParaView file for a load case
void COutputter::OutputVTK(unsigned int LoadCase)
{
	CDomain* FEMData = CDomain::GetInstance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();

	std::string base = OutputFileName_;
	std::size_t dot = base.find_last_of('.');
	if (dot != std::string::npos)
		base = base.substr(0, dot);

	std::ostringstream name;
	name << base << "_lc" << LoadCase << ".vtu";

	std::ofstream out(name.str());
	if (!out)
	{
		cerr << "*** Error *** Failed to open VTK file: " << name.str() << endl;
		return;
	}

	std::vector<int> connectivity;
	std::vector<int> offsets;
	std::vector<unsigned int> cellTypes;
	std::vector<unsigned int> elementTypes;
	std::vector<double> sxx, syy, szz, sxy, syz, sxz;
	std::vector<double> sxx2, syy2, sxy2;
	std::vector<double> mx, my, mxy;
	std::vector<double> vonMises;
	std::vector<double> axialForce;
	std::vector<double> beamMy, beamMz, beamTorque;

	int offset = 0;

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType();

		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			CElement& Element = EleGrp[Ele];
			CNode** nodes = Element.GetNodes();

			if (ElementType == ElementTypes::Bar)
			{
				connectivity.push_back(nodes[0]->NodeNumber - 1);
				connectivity.push_back(nodes[1]->NodeNumber - 1);
				offset += 2;
				offsets.push_back(offset);
				cellTypes.push_back(3); // VTK_LINE
			}
			else if (ElementType == ElementTypes::H8 || ElementType == ElementTypes::BbarH8)
			{
				for (int i = 0; i < 8; i++)
					connectivity.push_back(nodes[i]->NodeNumber - 1);
				offset += 8;
				offsets.push_back(offset);
				cellTypes.push_back(12); // VTK_HEXAHEDRON
			}
			else
			{
				continue;
			}

			elementTypes.push_back(static_cast<unsigned int>(ElementType));

			if (ElementType == ElementTypes::H8 || ElementType == ElementTypes::BbarH8)
			{
				double stress[6];
				Element.ElementStress(stress, Displacement);
				double vm = sqrt(0.5 * ((stress[0] - stress[1]) * (stress[0] - stress[1])
					+ (stress[1] - stress[2]) * (stress[1] - stress[2])
					+ (stress[2] - stress[0]) * (stress[2] - stress[0]))
					+ 3.0 * (stress[3] * stress[3] + stress[4] * stress[4] + stress[5] * stress[5]));

				sxx.push_back(stress[0]);
				syy.push_back(stress[1]);
				szz.push_back(stress[2]);
				sxy.push_back(stress[3]);
				syz.push_back(stress[4]);
				sxz.push_back(stress[5]);
				vonMises.push_back(vm);
				axialForce.push_back(0.0);
			}
			else if (ElementType == ElementTypes::Bar)
			{
				double stress;
				Element.ElementStress(&stress, Displacement);
				CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());

				sxx.push_back(0.0);
				syy.push_back(0.0);
				szz.push_back(0.0);
				sxy.push_back(0.0);
				syz.push_back(0.0);
				sxz.push_back(0.0);
				vonMises.push_back(0.0);
				axialForce.push_back(stress * material.Area);
			}
			else
			{
				sxx.push_back(0.0);
				syy.push_back(0.0);
				szz.push_back(0.0);
				sxy.push_back(0.0);
				syz.push_back(0.0);
				sxz.push_back(0.0);
				vonMises.push_back(0.0);
				axialForce.push_back(0.0);
			}

			sxx2.push_back(0.0);
			syy2.push_back(0.0);
			sxy2.push_back(0.0);
			mx.push_back(0.0);
			my.push_back(0.0);
			mxy.push_back(0.0);
			beamMy.push_back(0.0);
			beamMz.push_back(0.0);
			beamTorque.push_back(0.0);
		}
	}

	unsigned int numCells = static_cast<unsigned int>(cellTypes.size());

	out << "<?xml version=\"1.0\"?>\n";
	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "  <UnstructuredGrid>\n";
	out << "    <Piece NumberOfPoints=\"" << NUMNP << "\" NumberOfCells=\"" << numCells << "\">\n";

	out << "      <PointData>\n";
	out << "        <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < NUMNP; i++)
	{
		double ux = 0.0, uy = 0.0, uz = 0.0;
		if (NodeList[i].bcode[0]) ux = Displacement[NodeList[i].bcode[0] - 1];
		if (NodeList[i].bcode[1]) uy = Displacement[NodeList[i].bcode[1] - 1];
		if (NodeList[i].bcode[2]) uz = Displacement[NodeList[i].bcode[2] - 1];
		out << "          " << ux << " " << uy << " " << uz << "\n";
	}
	out << "        </DataArray>\n";

	out << "        <DataArray type=\"Float64\" Name=\"Rotation\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < NUMNP; i++)
		out << "          0 0 0\n";
	out << "        </DataArray>\n";

	out << "        <DataArray type=\"Float64\" Name=\"DisplacementMagnitude\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < NUMNP; i++)
	{
		double ux = 0.0, uy = 0.0, uz = 0.0;
		if (NodeList[i].bcode[0]) ux = Displacement[NodeList[i].bcode[0] - 1];
		if (NodeList[i].bcode[1]) uy = Displacement[NodeList[i].bcode[1] - 1];
		if (NodeList[i].bcode[2]) uz = Displacement[NodeList[i].bcode[2] - 1];
		out << "          " << sqrt(ux * ux + uy * uy + uz * uz) << "\n";
	}
	out << "        </DataArray>\n";
	out << "      </PointData>\n";

	out << "      <CellData>\n";
	out << "        <DataArray type=\"Int32\" Name=\"ElementType\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < numCells; i++)
		out << "          " << elementTypes[i] << "\n";
	out << "        </DataArray>\n";

	out << "        <DataArray type=\"Float64\" Name=\"VonMises\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < numCells; i++)
		out << "          " << vonMises[i] << "\n";
	out << "        </DataArray>\n";

	out << "        <DataArray type=\"Float64\" Name=\"SXX SYY SZZ SXY SYZ SXZ\" NumberOfComponents=\"6\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < numCells; i++)
		out << "          " << sxx[i] << " " << syy[i] << " " << szz[i] << " " << sxy[i] << " " << syz[i] << " " << sxz[i] << "\n";
	out << "        </DataArray>\n";

	out << "        <DataArray type=\"Float64\" Name=\"SXX SYY SXY\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < numCells; i++)
		out << "          " << sxx2[i] << " " << syy2[i] << " " << sxy2[i] << "\n";
	out << "        </DataArray>\n";

	out << "        <DataArray type=\"Float64\" Name=\"MX MY MXY\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < numCells; i++)
		out << "          " << mx[i] << " " << my[i] << " " << mxy[i] << "\n";
	out << "        </DataArray>\n";

	out << "        <DataArray type=\"Float64\" Name=\"AxialForce\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < numCells; i++)
		out << "          " << axialForce[i] << "\n";
	out << "        </DataArray>\n";

	out << "        <DataArray type=\"Float64\" Name=\"BeamMy BeamMz BeamTorque\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < numCells; i++)
		out << "          " << beamMy[i] << " " << beamMz[i] << " " << beamTorque[i] << "\n";
	out << "        </DataArray>\n";
	out << "      </CellData>\n";

	out << "      <Points>\n";
	out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < NUMNP; i++)
		out << "          " << NodeList[i].XYZ[0] << " " << NodeList[i].XYZ[1] << " " << NodeList[i].XYZ[2] << "\n";
	out << "        </DataArray>\n";
	out << "      </Points>\n";

	out << "      <Cells>\n";
	out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < connectivity.size(); i++)
		out << "          " << connectivity[i] << "\n";
	out << "        </DataArray>\n";

	out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < offsets.size(); i++)
		out << "          " << offsets[i] << "\n";
	out << "        </DataArray>\n";

	out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	for (unsigned int i = 0; i < cellTypes.size(); i++)
		out << "          " << cellTypes[i] << "\n";
	out << "        </DataArray>\n";
	
	out << "      </Cells>\n";
	out << "    </Piece>\n";
	out << "  </UnstructuredGrid>\n";
	out << "</VTKFile>\n";
}

//	Print total system data
void COutputter::OutputTotalSystemData()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "	TOTAL SYSTEM DATA" << endl
		  << endl;

	*this << "     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = " << FEMData->GetNEQ()
		  << endl
		  << "     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = " << FEMData->GetStiffnessMatrix()->size()
		  << endl
		  << "     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = " << FEMData->GetStiffnessMatrix()->GetMaximumHalfBandwidth()
		  << endl
		  << "     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = " << FEMData->GetStiffnessMatrix()->size() / FEMData->GetNEQ() << endl
		  << endl
		  << endl;
}

#ifdef _DEBUG_

//	Print column heights for debuging
void COutputter::PrintColumnHeights()
{
	*this << "*** _Debug_ *** Column Heights" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();

	for (unsigned int col = 0; col < NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << ColumnHeights[col];
	}

	*this << endl
		  << endl;
}

//	Print address of diagonal elements for debuging
void COutputter::PrintDiagonalAddress()
{
	*this << "*** _Debug_ *** Address of Diagonal Element" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	for (unsigned int col = 0; col <= NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << DiagonalAddress[col];
	}

	*this << endl
		  << endl;
}

//	Print banded and full stiffness matrix for debuging
void COutputter::PrintStiffnessMatrix()
{
	*this << "*** _Debug_ *** Banded stiffness matrix" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < DiagonalAddress[NEQ] - DiagonalAddress[0]; i++)
	{
		*this << setw(14) << (*StiffnessMatrix)(i);

		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}
	}

	*this << endl
		  << endl;

	*this << "*** _Debug_ *** Full stiffness matrix" << endl;

	for (int I = 1; I <= NEQ; I++)
	{
		for (int J = 1; J <= NEQ; J++)
		{
			int J_new = (J > I) ? J : I;
			int I_new = (J > I) ? I : J;
			int H = DiagonalAddress[J_new] - DiagonalAddress[J_new - 1];
			if (J_new - I_new - H >= 0)
			{
				*this << setw(14) << 0.0;
			}
			else
			{
				*this << setw(14) << (*StiffnessMatrix)(I_new, J_new);
			}
		}

		*this << endl;
	}

	*this << endl;
}

//	Print displacement vector for debuging
void COutputter::PrintDisplacement()
{
	*this << "*** _Debug_ *** Displacement vector" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	double* Force = FEMData->GetForce();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < NEQ; i++)
	{
		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}

		*this << setw(14) << Force[i];
	}

	*this << endl
		  << endl;
}

#endif
