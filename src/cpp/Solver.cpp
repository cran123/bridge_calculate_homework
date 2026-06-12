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
#include <cstdlib>
#include <string>

using namespace std;

CLDLTSolver::CLDLTSolver(CSkylineMatrix<double>* K, const vector<CMpcConstraint>& MpcConstraints)
    : K_(K)
    , MpcConstraints_(MpcConstraints)
#ifdef STAPPP_USE_EIGEN
    , SparseInput_(nullptr)
    , Factorized_(false)
    , UseLLT_(false)
    , UseSparseLU_(false)
    , UseCG_(false)
#ifdef STAPPP_USE_PARDISO
    , PardisoMaxfct_(1)
    , PardisoMnum_(1)
    , PardisoMtype_(-2)
    , PardisoN_(0)
    , PardisoMsglvl_(0)
    , PardisoError_(0)
    , UsePardiso_(false)
    , PardisoInitialized_(false)
#endif
#endif
{
#ifdef STAPPP_USE_PARDISO
    for (unsigned int i = 0; i < 64; i++)
        PardisoPt_[i] = nullptr;
    for (unsigned int i = 0; i < 64; i++)
        PardisoIparm_[i] = 0;
#endif
}

#ifdef STAPPP_USE_EIGEN
CLDLTSolver::CLDLTSolver(const Eigen::SparseMatrix<double>* K, const vector<CMpcConstraint>& MpcConstraints)
    : K_(nullptr)
    , MpcConstraints_(MpcConstraints)
    , SparseInput_(K)
    , Factorized_(false)
    , UseLLT_(false)
    , UseSparseLU_(false)
    , UseCG_(false)
#ifdef STAPPP_USE_PARDISO
    , PardisoMaxfct_(1)
    , PardisoMnum_(1)
    , PardisoMtype_(-2)
    , PardisoN_(0)
    , PardisoMsglvl_(0)
    , PardisoError_(0)
    , UsePardiso_(false)
    , PardisoInitialized_(false)
#endif
{
#ifdef STAPPP_USE_PARDISO
    for (unsigned int i = 0; i < 64; i++)
        PardisoPt_[i] = nullptr;
    for (unsigned int i = 0; i < 64; i++)
        PardisoIparm_[i] = 0;
#endif
}

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

    bool IsEnabledEnv(const char* value, const char* name)
    {
        if (!value)
            return false;
        string v(value);
        return v == name;
    }

    int EnvInt(const char* name, int fallback)
    {
        const char* value = getenv(name);
        if (!value)
            return fallback;
        return atoi(value);
    }

    double EnvDouble(const char* name, double fallback)
    {
        const char* value = getenv(name);
        if (!value)
            return fallback;
        return atof(value);
    }

#ifdef STAPPP_USE_PARDISO
    void FillPardisoCsr(const Eigen::SparseMatrix<double>& matrix, MKL_INT mtype,
        vector<MKL_INT>& ia, vector<MKL_INT>& ja, vector<double>& a)
    {
        typedef Eigen::SparseMatrix<double, Eigen::RowMajor> RowMajorSparseMatrix;
        RowMajorSparseMatrix rowMajor(matrix);
        rowMajor.makeCompressed();

        const MKL_INT n = static_cast<MKL_INT>(rowMajor.rows());
        const bool symmetricStorage = (mtype == 2 || mtype == -2);
        ia.resize(static_cast<size_t>(n) + 1);
        ja.clear();
        a.clear();
        ja.reserve(rowMajor.nonZeros());
        a.reserve(rowMajor.nonZeros());

        const int* outer = rowMajor.outerIndexPtr();
        const int* inner = rowMajor.innerIndexPtr();
        const double* values = rowMajor.valuePtr();

        ia[0] = 1;
        for (MKL_INT i = 0; i < n; i++)
        {
            for (int k = outer[i]; k < outer[i + 1]; k++)
            {
                if (symmetricStorage && inner[k] < i)
                    continue;
                ja.push_back(static_cast<MKL_INT>(inner[k]) + 1);
                a.push_back(values[k]);
            }
            ia[static_cast<size_t>(i) + 1] = static_cast<MKL_INT>(ja.size()) + 1;
        }
    }
#endif
}
#endif

CLDLTSolver::~CLDLTSolver()
{
#ifdef STAPPP_USE_PARDISO
    if (PardisoInitialized_)
    {
        MKL_INT phase = -1;
        MKL_INT nrhs = 0;
        double ddum = 0.0;
        MKL_INT idum = 0;
        pardiso(PardisoPt_, &PardisoMaxfct_, &PardisoMnum_, &PardisoMtype_, &phase,
            &PardisoN_, &ddum, PardisoIa_.data(), PardisoJa_.data(), &idum, &nrhs,
            PardisoIparm_, &PardisoMsglvl_, &ddum, &ddum, &PardisoError_);
    }
#endif
}

// LDLT facterization
void CLDLTSolver::LDLT()
{
#ifdef STAPPP_USE_EIGEN
    if (SparseInput_)
        SparseK_ = *SparseInput_;
    else
        SparseK_ = BuildEigenSparseMatrix(*K_);

    const unsigned int dimension = static_cast<unsigned int>(SparseK_.rows());
    Transform_ = BuildConstraintTransform(dimension, MpcConstraints_, ReducedIndex_);
    ReducedK_ = Transform_.transpose() * SparseK_ * Transform_;
    ReducedK_.makeCompressed();

    const char* solverMode = getenv("STAPPP_SOLVER");
#ifdef STAPPP_USE_PARDISO
    if (!solverMode || IsEnabledEnv(solverMode, "pardiso"))
    {
        UsePardiso_ = true;
        UseCG_ = false;
        UseLLT_ = false;
        UseSparseLU_ = false;

        PardisoN_ = static_cast<MKL_INT>(ReducedK_.rows());
        PardisoMtype_ = static_cast<MKL_INT>(EnvInt("STAPPP_PARDISO_MTYPE", -2));
        PardisoMsglvl_ = static_cast<MKL_INT>(EnvInt("STAPPP_PARDISO_MSGLVL", 0));
        FillPardisoCsr(ReducedK_, PardisoMtype_, PardisoIa_, PardisoJa_, PardisoA_);

        for (unsigned int i = 0; i < 64; i++)
        {
            PardisoPt_[i] = nullptr;
            PardisoIparm_[i] = 0;
        }
        pardisoinit(PardisoPt_, &PardisoMtype_, PardisoIparm_);
        PardisoIparm_[0] = 1;
        PardisoIparm_[1] = 2;
        PardisoIparm_[7] = 2;
        PardisoIparm_[9] = 13;
        PardisoIparm_[17] = -1;
        PardisoIparm_[18] = -1;
        PardisoIparm_[26] = 1;
        PardisoIparm_[34] = 0;
        PardisoIparm_[59] = EnvInt("STAPPP_PARDISO_OOC", 0);

        MKL_INT phase = 12;
        MKL_INT nrhs = 1;
        MKL_INT idum = 0;
        double ddum = 0.0;
        pardiso(PardisoPt_, &PardisoMaxfct_, &PardisoMnum_, &PardisoMtype_, &phase,
            &PardisoN_, PardisoA_.data(), PardisoIa_.data(), PardisoJa_.data(), &idum,
            &nrhs, PardisoIparm_, &PardisoMsglvl_, &ddum, &ddum, &PardisoError_);
        if (PardisoError_ != 0)
        {
            cerr << "*** Error *** Intel oneMKL PARDISO factorization failed. error = "
                 << PardisoError_ << endl;
            PrintReducedMatrixDiagnostics(ReducedK_, ReducedIndex_);
            exit(4);
        }

        PardisoInitialized_ = true;
        Factorized_ = true;
        cerr << "    Intel oneMKL PARDISO solver prepared."
             << " mtype=" << PardisoMtype_
             << " ooc=" << PardisoIparm_[59]
             << " nonzeros=" << PardisoA_.size() << endl;
        return;
    }
#endif
    if (IsEnabledEnv(solverMode, "pcg_ic"))
    {
        UseCG_ = true;
        UseLLT_ = false;
        UseSparseLU_ = false;
        CgSolver_.setTolerance(EnvDouble("STAPPP_CG_TOL", 1.0e-8));
        CgSolver_.setMaxIterations(EnvInt("STAPPP_CG_MAXITER", 5000));
        CgSolver_.compute(ReducedK_);
        if (CgSolver_.info() != Eigen::Success)
        {
            cerr << "*** Error *** Eigen PCG/IncompleteCholesky setup failed." << endl;
            PrintReducedMatrixDiagnostics(ReducedK_, ReducedIndex_);
            exit(4);
        }
        cerr << "    Eigen PCG/IncompleteCholesky solver prepared."
             << " max_iterations=" << CgSolver_.maxIterations()
             << " tolerance=" << CgSolver_.tolerance() << endl;
        Factorized_ = true;
        return;
    }

    UseCG_ = false;
    if (IsEnabledEnv(solverMode, "direct_llt") || IsEnabledEnv(solverMode, "llt"))
    {
        LLTSolver_.analyzePattern(ReducedK_);
        LLTSolver_.factorize(ReducedK_);
        if (LLTSolver_.info() == Eigen::Success)
        {
            cerr << "    Eigen SimplicialLLT direct solver prepared." << endl;
            UseLLT_ = true;
            UseSparseLU_ = false;
            Factorized_ = true;
            return;
        }

        cerr << "*** Warning *** Eigen SimplicialLLT factorization failed; falling back to SimplicialLDLT." << endl;
    }

    UseCG_ = false;
    UseLLT_ = false;
    Solver_.analyzePattern(ReducedK_);
    Solver_.factorize(ReducedK_);
    if (Solver_.info() == Eigen::Success)
    {
        UseSparseLU_ = false;
        Factorized_ = true;
        return;
    }

    cerr << "*** Warning *** Eigen SimplicialLDLT factorization failed; falling back to SparseLU." << endl;
    PrintReducedMatrixDiagnostics(ReducedK_, ReducedIndex_);
    SparseLUSolver_.analyzePattern(ReducedK_);
    SparseLUSolver_.factorize(ReducedK_);
    if (SparseLUSolver_.info() != Eigen::Success)
    {
        cerr << "*** Error *** Eigen SparseLU factorization failed." << endl;
        exit(4);
    }

    UseSparseLU_ = true;
    Factorized_ = true;
    return;
#else
	unsigned int N = K_->dim();
    unsigned int* ColumnHeights = K_->GetColumnHeights();   // Column Hights

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
				C += (*K_)(r,i) * (*K_)(r,j);		// C += L_ri * U_rj

			(*K_)(i,j) -= C;	// U_ij = K_ij - C
		}

		for (unsigned int r = mj; r <= j-1; r++)	// Loop for mj:j-1 (column j)
		{
			double Lrj = (*K_)(r,j) / (*K_)(r,r);	// L_rj = U_rj / D_rr
			(*K_)(j,j) -= Lrj * (*K_)(r,j);	// D_jj = K_jj - sum(L_rj*U_rj, r=mj:j-1)
			(*K_)(r,j) = Lrj;
		}

        if (fabs((*K_)(j,j)) <= FLT_MIN)
        {
            cerr << "*** Error *** Stiffness matrix is not positive definite !" << endl
            	 << "    Euqation no = " << j << endl
            	 << "    Pivot = " << (*K_)(j,j) << endl;
            
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

    const unsigned int N = static_cast<unsigned int>(SparseK_.rows());
    Eigen::VectorXd rhs(N);
    for (unsigned int i = 0; i < N; i++)
        rhs[i] = Force[i];

    Eigen::VectorXd reducedRhs = Transform_.transpose() * rhs;
    Eigen::VectorXd solution;
    Eigen::ComputationInfo info;
#ifdef STAPPP_USE_PARDISO
    if (UsePardiso_)
    {
        solution = Eigen::VectorXd::Zero(reducedRhs.size());
        MKL_INT phase = 33;
        MKL_INT nrhs = 1;
        MKL_INT idum = 0;
        pardiso(PardisoPt_, &PardisoMaxfct_, &PardisoMnum_, &PardisoMtype_, &phase,
            &PardisoN_, PardisoA_.data(), PardisoIa_.data(), PardisoJa_.data(), &idum,
            &nrhs, PardisoIparm_, &PardisoMsglvl_, reducedRhs.data(), solution.data(),
            &PardisoError_);
        if (PardisoError_ != 0)
        {
            cerr << "*** Error *** Intel oneMKL PARDISO solve failed. error = "
                 << PardisoError_ << endl;
            exit(4);
        }
        info = Eigen::Success;
    }
    else
#endif
    if (UseCG_)
    {
        solution = CgSolver_.solve(reducedRhs);
        info = CgSolver_.info();
        cerr << "    Eigen PCG iterations = " << CgSolver_.iterations()
             << ", estimated error = " << CgSolver_.error() << endl;
    }
    else if (UseLLT_)
    {
        solution = LLTSolver_.solve(reducedRhs);
        info = LLTSolver_.info();
    }
    else if (UseSparseLU_)
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

    if (getenv("STAPPP_CHECK_RESIDUAL"))
    {
        Eigen::VectorXd residual = ReducedK_ * solution - reducedRhs;
        double relative = residual.norm() / max(1.0, reducedRhs.norm());
        cerr << "    Reduced equation relative residual = " << scientific << setprecision(6)
             << relative << defaultfloat << endl;
    }

    Eigen::VectorXd fullSolution = Transform_ * solution;
    for (unsigned int i = 0; i < N; i++)
        Force[i] = fullSolution[i];
    return;
#else
	unsigned int N = K_->dim();
    unsigned int* ColumnHeights = K_->GetColumnHeights();   // Column Hights

//	Reduce right-hand-side load vector (LV = R)
	for (unsigned int i = 2; i <= N; i++)	// Loop for i=2:N (Numering starting from 1)
	{
        unsigned int mi = i - ColumnHeights[i-1];

		for (unsigned int j = mi; j <= i-1; j++)	// Loop for j=mi:i-1
			Force[i-1] -= (*K_)(j,i) * Force[j-1];	// V_i = R_i - sum_j (L_ji V_j)
	}

//	Back substitute (Vbar = D^(-1) V, L^T a = Vbar)
	for (unsigned int i = 1; i <= N; i++)	// Loop for i=1:N
		Force[i-1] /= (*K_)(i,i);	// Vbar = D^(-1) V

	for (unsigned int j = N; j >= 2; j--)	// Loop for j=N:2
	{
        unsigned int mj = j - ColumnHeights[j-1];

		for (unsigned int i = mj; i <= j-1; i++)	// Loop for i=mj:j-1
			Force[i-1] -= (*K_)(i,j) * Force[j-1];	// a_i = Vbar_i - sum_j(L_ij Vbar_j)
	}
#endif
};
