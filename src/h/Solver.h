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

#ifdef STAPPP_USE_EIGEN
#include <Eigen/Sparse>
#endif

//!	LDLT solver: A in core solver using skyline storage  and column reduction scheme
class CLDLTSolver
{
private:
    
    CSkylineMatrix<double>& K;

#ifdef STAPPP_USE_EIGEN
    Eigen::SparseMatrix<double> SparseK_;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > Solver_;
    bool Factorized_;
#endif

public:

//!	Constructor
	CLDLTSolver(CSkylineMatrix<double>* K);

//!	Perform L*D*L(T) factorization of the stiffness matrix
	void LDLT();

//!	Reduce right-hand-side load vector and back substitute
	void BackSubstitution(double* Force); 
};
