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

#include "SkylineMatrix.h"
#include "MpcConstraint.h"

#ifdef STAPPP_USE_EIGEN
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#endif

#include <vector>

//!	LDLT solver: A in core solver using skyline storage  and column reduction scheme
class CLDLTSolver
{
private:
    
    CSkylineMatrix<double>* K_;
    const std::vector<CMpcConstraint>& MpcConstraints_;

#ifdef STAPPP_USE_EIGEN
    const Eigen::SparseMatrix<double>* SparseInput_;
    Eigen::SparseMatrix<double> SparseK_;
    Eigen::SparseMatrix<double> ReducedK_;
    Eigen::SparseMatrix<double> Transform_;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > Solver_;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int> > LLTSolver_;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SparseLUSolver_;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper,
        Eigen::IncompleteCholesky<double> > CgSolver_;
    std::vector<int> ReducedIndex_;
    bool Factorized_;
    bool UseLLT_;
    bool UseSparseLU_;
    bool UseCG_;
#endif

public:

//!	Constructor
	CLDLTSolver(CSkylineMatrix<double>* K, const std::vector<CMpcConstraint>& MpcConstraints);

#ifdef STAPPP_USE_EIGEN
//!	Constructor for direct Eigen sparse assembly
	CLDLTSolver(const Eigen::SparseMatrix<double>* K, const std::vector<CMpcConstraint>& MpcConstraints);
#endif

//!	Perform L*D*L(T) factorization of the stiffness matrix
	void LDLT();

//!	Reduce right-hand-side load vector and back substitute
	void BackSubstitution(double* Force); 
};
