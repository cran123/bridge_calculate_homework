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
#include <iomanip>
#include <algorithm>
#include <map>
#include <set>
#include <vector>

using namespace std;

CLDLTSolver::CLDLTSolver(CSkylineMatrix<double>* K, const vector<CMpcConstraint>& MpcConstraints)
    : K(*K)
    , MpcConstraints_(MpcConstraints)
#ifdef STAPPP_USE_EIGEN
    , Factorized_(false)
    , UseSparseLU_(false)
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

    Eigen::SparseMatrix<double> BuildConstraintTransform(unsigned int N,
        const vector<CMpcConstraint>& constraints, vector<int>& reducedIndex)
    {
        vector<int> slaveOf(N, -1);
        vector<vector<pair<unsigned int, double> > > slaveRelations(N);

        for (unsigned int c = 0; c < constraints.size(); c++)
        {
            const CMpcConstraint& constraint = constraints[c];
            int slave = -1;
            double slaveCoef = 0.0;
            for (unsigned int i = 0; i < constraint.terms.size(); i++)
            {
                unsigned int eq = constraint.equations[i];
                if (!eq)
                    continue;
                double coef = constraint.terms[i].coefficient;
                if (fabs(coef - 1.0) <= 1.0e-12)
                {
                    slave = static_cast<int>(eq - 1);
                    slaveCoef = coef;
                    break;
                }
            }

            if (slave < 0)
                continue;
            if (slaveOf[slave] >= 0)
            {
                cerr << "*** Error *** MPC slave equation appears in more than one constraint." << endl;
                exit(4);
            }

            slaveOf[slave] = static_cast<int>(c);
            for (unsigned int i = 0; i < constraint.terms.size(); i++)
            {
                unsigned int eq = constraint.equations[i];
                if (!eq || static_cast<int>(eq - 1) == slave)
                    continue;
                double coef = -constraint.terms[i].coefficient / slaveCoef;
                slaveRelations[slave].push_back(make_pair(eq - 1, coef));
            }
        }

        reducedIndex.assign(N, -1);
        unsigned int reduced = 0;
        for (unsigned int eq = 0; eq < N; eq++)
        {
            if (slaveOf[eq] < 0)
                reducedIndex[eq] = reduced++;
        }

        vector<Eigen::Triplet<double> > triplets;
        triplets.reserve(N + constraints.size() * 4);

        vector<int> state(N, 0);
        vector<vector<pair<unsigned int, double> > > expandedRelations(N);

        struct Expander
        {
            const vector<int>& slaveOf;
            const vector<vector<pair<unsigned int, double> > >& slaveRelations;
            vector<int>& state;
            vector<vector<pair<unsigned int, double> > >& expandedRelations;

            vector<pair<unsigned int, double> > Expand(unsigned int eq)
            {
                if (slaveOf[eq] < 0)
                    return vector<pair<unsigned int, double> >(1, make_pair(eq, 1.0));

                if (state[eq] == 2)
                    return expandedRelations[eq];
                if (state[eq] == 1)
                {
                    cerr << "*** Error *** Cyclic MPC dependency detected." << endl;
                    exit(4);
                }

                state[eq] = 1;
                map<unsigned int, double> merged;
                for (unsigned int i = 0; i < slaveRelations[eq].size(); i++)
                {
                    unsigned int masterEq = slaveRelations[eq][i].first;
                    double coef = slaveRelations[eq][i].second;
                    vector<pair<unsigned int, double> > terms = Expand(masterEq);
                    for (unsigned int j = 0; j < terms.size(); j++)
                        merged[terms[j].first] += coef * terms[j].second;
                }

                for (map<unsigned int, double>::iterator it = merged.begin(); it != merged.end(); ++it)
                    if (fabs(it->second) > 1.0e-14)
                        expandedRelations[eq].push_back(make_pair(it->first, it->second));

                state[eq] = 2;
                return expandedRelations[eq];
            }
        };

        Expander expander = { slaveOf, slaveRelations, state, expandedRelations };

        for (unsigned int eq = 0; eq < N; eq++)
        {
            if (slaveOf[eq] < 0)
            {
                triplets.push_back(Eigen::Triplet<double>(eq, reducedIndex[eq], 1.0));
                continue;
            }

            vector<pair<unsigned int, double> > terms = expander.Expand(eq);
            for (unsigned int i = 0; i < terms.size(); i++)
            {
                unsigned int masterEq = terms[i].first;
                double coef = terms[i].second;
                int col = reducedIndex[masterEq];
                triplets.push_back(Eigen::Triplet<double>(eq, col, coef));
            }
        }

        Eigen::SparseMatrix<double> transform(N, reduced);
        transform.setFromTriplets(triplets.begin(), triplets.end());
        transform.makeCompressed();
        return transform;
    }

    void PrintReducedMatrixDiagnostics(const Eigen::SparseMatrix<double>& matrix,
        const vector<int>& reducedIndex)
    {
        vector<double> diagonal(matrix.rows(), 0.0);
        for (int col = 0; col < matrix.outerSize(); col++)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, col); it; ++it)
            {
                if (it.row() == it.col())
                    diagonal[it.row()] = it.value();
            }
        }

        double maxAbsDiag = 0.0;
        double minDiag = diagonal.empty() ? 0.0 : diagonal[0];
        double minAbsDiag = diagonal.empty() ? 0.0 : fabs(diagonal[0]);
        int nonPositiveDiag = 0;
        for (unsigned int i = 0; i < diagonal.size(); i++)
        {
            maxAbsDiag = max(maxAbsDiag, fabs(diagonal[i]));
            minDiag = min(minDiag, diagonal[i]);
            minAbsDiag = min(minAbsDiag, fabs(diagonal[i]));
            if (diagonal[i] <= 0.0)
                nonPositiveDiag++;
        }

        const double tinyLimit = max(1.0, maxAbsDiag) * 1.0e-12;
        int tinyDiag = 0;
        for (unsigned int i = 0; i < diagonal.size(); i++)
            if (fabs(diagonal[i]) <= tinyLimit)
                tinyDiag++;

        vector<unsigned int> sourceCounts(matrix.rows(), 0);
        for (unsigned int eq = 0; eq < reducedIndex.size(); eq++)
        {
            int red = reducedIndex[eq];
            if (red >= 0)
                sourceCounts[red]++;
        }

        vector<pair<double, int> > suspicious;
        suspicious.reserve(diagonal.size());
        for (unsigned int i = 0; i < diagonal.size(); i++)
            suspicious.push_back(make_pair(fabs(diagonal[i]), static_cast<int>(i)));
        sort(suspicious.begin(), suspicious.end());

        cerr << scientific << setprecision(6);
        cerr << "    Reduced matrix dimension = " << matrix.rows()
             << ", nonzeros = " << matrix.nonZeros() << endl;
        cerr << "    Reduced diagonal min = " << minDiag
             << ", min(abs) = " << minAbsDiag
             << ", max(abs) = " << maxAbsDiag << endl;
        cerr << "    Non-positive diagonal count = " << nonPositiveDiag
             << ", tiny diagonal count(|d| <= " << tinyLimit << ") = " << tinyDiag << endl;
        cerr << "    Smallest reduced diagonal entries:" << endl;
        for (unsigned int i = 0; i < suspicious.size() && i < 12; i++)
        {
            int red = suspicious[i].second;
            cerr << "      reduced eq " << red + 1
                 << ": diag = " << diagonal[red]
                 << ", source active dofs = " << sourceCounts[red] << endl;
        }
    }
}
#endif

// LDLT facterization
void CLDLTSolver::LDLT()
{
#ifdef STAPPP_USE_EIGEN
    SparseK_ = BuildEigenSparseMatrix(K);
    Transform_ = BuildConstraintTransform(K.dim(), MpcConstraints_, ReducedIndex_);
    Eigen::SparseMatrix<double> reducedK = Transform_.transpose() * SparseK_ * Transform_;
    reducedK.makeCompressed();

    Solver_.analyzePattern(reducedK);
    Solver_.factorize(reducedK);
    if (Solver_.info() == Eigen::Success)
    {
        UseSparseLU_ = false;
        Factorized_ = true;
        return;
    }

    cerr << "*** Warning *** Eigen SimplicialLDLT factorization failed; falling back to SparseLU." << endl;
    PrintReducedMatrixDiagnostics(reducedK, ReducedIndex_);
    SparseLUSolver_.analyzePattern(reducedK);
    SparseLUSolver_.factorize(reducedK);
    if (SparseLUSolver_.info() != Eigen::Success)
    {
        cerr << "*** Error *** Eigen SparseLU factorization failed." << endl;
        exit(4);
    }

    UseSparseLU_ = true;
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

    Eigen::VectorXd reducedRhs = Transform_.transpose() * rhs;
    Eigen::VectorXd solution;
    Eigen::ComputationInfo info;
    if (UseSparseLU_)
    {
        solution = SparseLUSolver_.solve(reducedRhs);
        info = SparseLUSolver_.info();
    }
    else
    {
        solution = Solver_.solve(reducedRhs);
        info = Solver_.info();
    }
    if (info != Eigen::Success)
    {
        cerr << "*** Error *** Eigen sparse solve failed." << endl;
        exit(4);
    }

    Eigen::VectorXd fullSolution = Transform_ * solution;
    for (unsigned int i = 0; i < N; i++)
        Force[i] = fullSolution[i];
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
