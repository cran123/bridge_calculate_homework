/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Bar.h"
#include "Outputter.h"
#include "Clock.h"
#include "VTKOutputter.h"

#include <cstdlib>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 2) //  Print help message
	{
	    cout << "Usage: stap++ InputFileName\n";
		exit(1);
	}

	string filename(argv[1]);
    size_t found = filename.find_last_of('.');

    // If the input file name is provided with an extension
    if (found != std::string::npos) {
        if (filename.substr(found) == ".dat")
            filename = filename.substr(0, found);
        else {
            // The input file name must has an extension of 'dat'
            cout << "*** Error *** Invalid file extension: "
                 << filename.substr(found+1) << endl;
            exit(1);
        }
    }

    string InFile = filename + ".dat";
	string OutFile = filename + ".out";

	CDomain* FEMData = CDomain::GetInstance();

    Clock timer;
    timer.Start();

//  Read data and define the problem domain
	if (!FEMData->ReadData(InFile, OutFile))
	{
		cerr << "*** Error *** Data input failed!" << endl;
		exit(1);
	}
    
    double time_input = timer.ElapsedTime();

    COutputter* Output = COutputter::GetInstance();

#ifdef STAPPP_USE_EIGEN
    const char* sparseEnv = std::getenv("STAPPP_SPARSE_ASSEMBLY");
    bool useSparseAssembly = false;
    if (sparseEnv)
    {
        string value(sparseEnv);
        useSparseAssembly = value != "0" && value != "false" && value != "FALSE";
    }
#else
    bool useSparseAssembly = false;
#endif

    if (!FEMData->GetMODEX())
    {
        *Output << "Data check completed !" << endl << endl;
        return 0;
    }

#ifdef STAPPP_USE_EIGEN
    Eigen::SparseMatrix<double> SparseStiffnessMatrix;
#endif
    CLDLTSolver* Solver = nullptr;

//  Allocate global vectors and matrices, such as the Force, ColumnHeights,
//  DiagonalAddress and StiffnessMatrix, and calculate the column heights
//  and address of diagonal elements
    if (useSparseAssembly)
    {
#ifdef STAPPP_USE_EIGEN
        FEMData->AllocateForceVector();
        FEMData->AssembleSparseStiffnessMatrix(SparseStiffnessMatrix);
        FEMData->ApplyDisplacementPenalty(SparseStiffnessMatrix);

        *Output << "	TOTAL SPARSE SYSTEM DATA" << endl
                << endl
                << "     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = "
                << FEMData->GetNEQ() << endl
                << "     NUMBER OF SPARSE NONZEROS . . . . . . . . . . .      = "
                << SparseStiffnessMatrix.nonZeros() << endl
                << "     ESTIMATED SPARSE STORAGE BYTES . . . . . . . . .     = "
                << static_cast<unsigned long long>(SparseStiffnessMatrix.nonZeros()) *
                       (sizeof(double) + sizeof(int))
                   + static_cast<unsigned long long>(SparseStiffnessMatrix.outerSize() + 1) * sizeof(int)
                << endl
                << endl
                << endl;

        if (FEMData->GetMODEX() == 2)
        {
            *Output << "Sparse matrix check completed !" << endl << endl;
            return 0;
        }

        Solver = new CLDLTSolver(&SparseStiffnessMatrix, FEMData->GetMpcConstraints());
#else
        cerr << "*** Error *** Direct sparse assembly requires Eigen solver build." << endl;
        exit(4);
#endif
    }
    else
    {
	    FEMData->AllocateMatrices();

        if (FEMData->GetMODEX() == 2)
        {
            *Output << "Matrix check completed !" << endl << endl;
            return 0;
        }
    
//  Assemble the banded gloabl stiffness matrix
	    FEMData->AssembleStiffnessMatrix();

//  Apply penalty terms for prescribed displacements
	    FEMData->ApplyDisplacementPenalty();

        Solver = new CLDLTSolver(FEMData->GetStiffnessMatrix(), FEMData->GetMpcConstraints());
    }
    
    double time_assemble = timer.ElapsedTime();

//  Solve the linear equilibrium equations for displacements
    
//  Perform L*D*L(T) factorization of stiffness matrix
    Solver->LDLT();

    double time_factorization = timer.ElapsedTime();

#ifdef _DEBUG_
    Output->PrintStiffnessMatrix();
#endif

    double time_load_solution = 0.0;
    double time_output = 0.0;

//  Loop over for all load cases
    for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
    {
        double time_case_start = timer.ElapsedTime();

//      Assemble righ-hand-side vector (force vector)
        FEMData->AssembleForce(lcase + 1);
            
//      Reduce right-hand-side force vector and back substitute
        Solver->BackSubstitution(FEMData->GetForce());

        double time_case_solved = timer.ElapsedTime();
        time_load_solution += time_case_solved - time_case_start;

        *Output << " LOAD CASE" << setw(5) << lcase + 1 << endl << endl << endl;

#ifdef _DEBUG_
        Output->PrintDisplacement();
#endif
            
        Output->OutputNodalDisplacement();

//      Calculate and output stresses of all elements
        Output->OutputElementStress();

        string VTKFile = filename + "_lc" + to_string(lcase + 1) + ".vtu";
        CVTKOutputter::WriteVTU(VTKFile, lcase + 1);

        double time_case_output = timer.ElapsedTime();
        time_output += time_case_output - time_case_solved;
    }

    double time_solution = timer.ElapsedTime();
    
    timer.Stop();
    
    *Output << "\n S O L U T I O N   T I M E   L O G   I N   S E C \n\n"
            << "     TIME FOR INPUT PHASE = " << time_input << endl
            << "     TIME FOR CALCULATION OF STIFFNESS MATRIX = " << time_assemble - time_input << endl
            << "     TIME FOR FACTORIZATION = " << time_factorization - time_assemble << endl
            << "     TIME FOR LOAD CASE SOLUTIONS = " << time_load_solution << endl
            << "     TIME FOR RESULT OUTPUT = " << time_output << endl << endl
            << "     T O T A L   S O L U T I O N   T I M E = " << time_solution << endl << endl;

    delete Solver;

	return 0;
}
