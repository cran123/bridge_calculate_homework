/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Solver.h"

#include <cmath>
#include <cfloat>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

CLDLTSolver::CLDLTSolver(CSkylineMatrix<double>* K)
    : K(*K)
#ifdef STAPPP_USE_EIGEN
    , Factorized_(false)
#endif
{
}

#ifdef STAPPP_USE_EIGEN
namespace
{
    Eigen::SparseMatrix<double> BuildEigenSparseMatrix(CSkylineMatrix<double>& K)
    {
        const unsigned int N = K.dim();
        unsigned int* ColumnHeights = K.GetColumnHeights();
        vector<Eigen::Triplet<double> > triplets;

        for (unsigned int j = 1; j <= N; j++)
        {
            unsigned int first = j - ColumnHeights[j - 1];
            for (unsigned int i = first; i <= j; i++)
            {
                double value = K(i, j);
                if (fabs(value) <= 0.0)
                    continue;

                triplets.push_back(Eigen::Triplet<double>(i - 1, j - 1, value));
                if (i != j)
                    triplets.push_back(Eigen::Triplet<double>(j - 1, i - 1, value));
            }
        }

        Eigen::SparseMatrix<double> sparse(N, N);
        sparse.setFromTriplets(triplets.begin(), triplets.end());
        sparse.makeCompressed();
        return sparse;
    }
}
#endif

// LDLT facterization
void CLDLTSolver::LDLT()
{
#ifdef STAPPP_USE_EIGEN
    SparseK_ = BuildEigenSparseMatrix(K);
    Solver_.analyzePattern(SparseK_);
    Solver_.factorize(SparseK_);
    if (Solver_.info() != Eigen::Success)
    {
        cerr << "*** Error *** Eigen SimplicialLDLT factorization failed." << endl;
        exit(4);
    }

    Factorized_ = true;
    return;
#else
	unsigned int N = K.dim();
    unsigned int* ColumnHeights = K.GetColumnHeights();   // Column Hights

	for (unsigned int j = 2; j <= N; j++)      // Loop for column 2:n (Numbering starting from 1)
	{
        // Row number of the first non-zero element in column j (Numbering starting from 1)
		unsigned int mj = j - ColumnHeights[j-1];
        
		for (unsigned int i = mj+1; i <= j-1; i++)	// Loop for mj+1:j-1 (Numbering starting from 1)
		{
            // Row number of the first nonzero element in column i (Numbering starting from 1)
			unsigned int mi = i - ColumnHeights[i-1];

			double C = 0.0;
			for (unsigned int r = max(mi, mj); r <= i-1; r++)
				C += K(r,i) * K(r,j);		// C += L_ri * U_rj

			K(i,j) -= C;	// U_ij = K_ij - C
		}

		for (unsigned int r = mj; r <= j-1; r++)	// Loop for mj:j-1 (column j)
		{
			double Lrj = K(r,j) / K(r,r);	// L_rj = U_rj / D_rr
			K(j,j) -= Lrj * K(r,j);	// D_jj = K_jj - sum(L_rj*U_rj, r=mj:j-1)
			K(r,j) = Lrj;
		}

        if (fabs(K(j,j)) <= FLT_MIN)
        {
            cerr << "*** Error *** Stiffness matrix is not positive definite !" << endl
            	 << "    Euqation no = " << j << endl
            	 << "    Pivot = " << K(j,j) << endl;
            
            exit(4);
        }
    }
#endif
};

// Solve displacement by back substitution
void CLDLTSolver::BackSubstitution(double* Force)
{
#ifdef STAPPP_USE_EIGEN
    if (!Factorized_)
        LDLT();

    const unsigned int N = K.dim();
    Eigen::VectorXd rhs(N);
    for (unsigned int i = 0; i < N; i++)
        rhs[i] = Force[i];

    Eigen::VectorXd solution = Solver_.solve(rhs);
    if (Solver_.info() != Eigen::Success)
    {
        cerr << "*** Error *** Eigen SimplicialLDLT solve failed." << endl;
        exit(4);
    }

    for (unsigned int i = 0; i < N; i++)
        Force[i] = solution[i];
    return;
#else
	unsigned int N = K.dim();
    unsigned int* ColumnHeights = K.GetColumnHeights();   // Column Hights

//	Reduce right-hand-side load vector (LV = R)
	for (unsigned int i = 2; i <= N; i++)	// Loop for i=2:N (Numering starting from 1)
	{
        unsigned int mi = i - ColumnHeights[i-1];

		for (unsigned int j = mi; j <= i-1; j++)	// Loop for j=mi:i-1
			Force[i-1] -= K(j,i) * Force[j-1];	// V_i = R_i - sum_j (L_ji V_j)
	}

//	Back substitute (Vbar = D^(-1) V, L^T a = Vbar)
	for (unsigned int i = 1; i <= N; i++)	// Loop for i=1:N
		Force[i-1] /= K(i,i);	// Vbar = D^(-1) V

	for (unsigned int j = N; j >= 2; j--)	// Loop for j=N:2
	{
        unsigned int mj = j - ColumnHeights[j-1];

		for (unsigned int i = mj; i <= j-1; i++)	// Loop for i=mj:j-1
			Force[i-1] -= K(i,j) * Force[j-1];	// a_i = Vbar_i - sum_j(L_ij Vbar_j)
	}
#endif
};
