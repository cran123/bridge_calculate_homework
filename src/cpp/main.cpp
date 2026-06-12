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
    cerr << "Input read. Results will be saved to " << OutFile << endl;
    
    double time_input = timer.ElapsedTime();

    COutputter* Output = COutputter::GetInstance();
    bool modelOutput = Output->ModelOutputEnabled();
    bool resultOutput = Output->ResultOutputEnabled();
    bool vtkOutput = Output->VtkOutputEnabled();

#ifdef STAPPP_USE_EIGEN
    const char* sparseEnv = std::getenv("STAPPP_SPARSE_ASSEMBLY");
#ifdef STAPPP_USE_PARDISO
    bool useSparseAssembly = true;
#else
    bool useSparseAssembly = false;
#endif
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

        if (modelOutput)
        {
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
        }

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
    cerr << "Stiffness matrix assembled." << endl;

//  Solve the linear equilibrium equations for displacements
    
//  Perform L*D*L(T) factorization of stiffness matrix
    Solver->LDLT();

    double time_factorization = timer.ElapsedTime();
    cerr << "Factorization completed." << endl;

#ifdef _DEBUG_
    Output->PrintStiffnessMatrix();
#endif

    double time_load_solution = 0.0;
    double time_output = 0.0;
    double time_solve_completed = time_factorization;

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
        time_solve_completed = time_case_solved;

        if (resultOutput)
            *Output << " LOAD CASE" << setw(5) << lcase + 1 << endl << endl << endl;

#ifdef _DEBUG_
        Output->PrintDisplacement();
#endif
            
        if (resultOutput)
            Output->OutputNodalDisplacement();

//      Calculate and output stresses of all elements
        if (resultOutput)
            Output->OutputElementStress();

        if (vtkOutput)
        {
            string VTKFile = filename + "_lc" + to_string(lcase + 1) + ".vtu";
            CVTKOutputter::WriteVTU(VTKFile, lcase + 1);
            cerr << "Load case " << lcase + 1 << " solved. VTK saved to " << VTKFile << endl;
        }
        else
        {
            cerr << "Load case " << lcase + 1 << " solved." << endl;
        }

        double time_case_output = timer.ElapsedTime();
        time_output += time_case_output - time_case_solved;
    }

    double time_solution = timer.ElapsedTime();
    
    timer.Stop();

    double competition_solve_time = time_solve_completed - time_input;
    
    (void)time_load_solution;
    (void)time_output;
    (void)time_solution;

    *Output << "\n P E R F O R M A N C E   M E T R I C S \n\n"
            << "     COMPETITION SOLVE WALL TIME AFTER INPUT BEFORE OUTPUT, SEC = "
            << competition_solve_time << endl << endl;

    cerr << "Competition solve wall time after input before output = "
         << competition_solve_time << " sec" << endl;
    cerr << "Done. Report saved to " << OutFile << endl;

    delete Solver;

	return 0;
}
