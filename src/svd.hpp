#ifndef __SVD_HPP__
#define __SVD_HPP__

#include <lapacke_utils.h>
#include "matrix.hpp"
#include "ops.hpp"
#include <complex>

extern "C" void dgesvd_(char* jobu, char* jobvt, int* m, int* n,
double* A, int* lda, double *S, double *U,
int* ldu, double *VT, int* ldvt, double work[],
int* lwork, int* info);

//extern "C" void dgesvdx_(char* jobu, char* jobvt, char* range, int* m, int* n,
//double* A, int* lda, double* VL, double* VU, int* IL, int* IU, int* NS, double *S, double *U,
//int* ldu, double *VT, int* ldvt, double work[],
//int* lwork, int iwork[], int* info);


template<typename T>
int computeSVDUsingLAPACK(matrix<T>& A, matrix<T>& U, matrix<T>& Sigma, matrix<T>& VT) {
#ifdef CheckSVD
    pout << "SVD matrix: \n";
	pout << A << "\n";
#endif

    matrix<T> Mat_RearrA = arrangeElementsVerticallyForLAPACK(A);

	char computeSpecific = 'S'; //Min(m,n) Columns of U are computed
	char computeAll = 'A';
	int rows = Mat_RearrA.rowCount();
	int cols = Mat_RearrA.colCount();
    matrix<T>S(Sigma.rowCount(), 1,0);

	int info;
	int lwork = MAX(MAX(1,3*MIN(rows,cols) + MAX(rows,cols)),5*MIN(rows,cols));
//	int lwork = MAX(1,MIN(rows,cols)*(MIN(rows,cols)+4));
	double work[lwork];
//	std::cout << "MIN(rows,cols): " << MIN(rows,cols) << "\n";
	dgesvd_(&computeAll, &computeAll, &rows, &cols, &Mat_RearrA(0,0), &rows, &S(0,0), &U(0,0), &rows, &VT(0,0), &cols, work, &lwork, &info);
//	int iworklength =12*MAX(rows,cols);
//	int iworklength =rows+ 3*cols;
//	int iwork[iworklength];
//	char V = 'V';
//	char I = 'I';
//	char chA = 'A';
//	double VL = 1;
//	double VU = 2;
//	int IL = 1;
//	int IU = 2;
//	int NS = 2;
//
//	dgesvdx_(&V, &V, &chA, &rows, &cols, &A(0,0),&rows,&VL,&VU,&IL,&IU,&NS, &S(0,0), &U(0,0), &rows, &VT(0,0), &cols, work, &lwork, iwork, &info);
	if (info != 0) {
		std::cout << "DGESVD returned info: " << info << "\n";
	}
	U = transpose(U);
	VT = transpose(VT);
  fillDiagonal(Sigma, S);

#ifdef Print
	pout << "EigenValues: \n";
	pout << S << "\n";
	pout << "U: \n";
	pout << U << "\n";
	pout << "VT: \n";
	pout << VT << "\n";
#endif

	for(unsigned int i=0; i<S.rowCount(); ++i){
		
		if(S(i,0)<1e-20){
			S(i,0) = 0.0;
		}
		pout << S(i,0) << "\n";
	}

#ifdef CheckSVD
	pout << "#Check SVD: \n";
    matrix<T> SVD = U*Sigma*VT;
    pout << SVD << "\n";
#endif

	S.destroy();
	return info;
}

extern "C" void dgeev_(char* jobvl, char* jobvr, int* n, double* A, int* lda, double* wr, double* wi, double* vl, int* ldvl, double* vr, 
int* ldvr, double work[], int* lwork, int* info);
template<typename T>
int computeEigenValuesUsingLAPACK(matrix<T>& A, matrix<T>& LR, matrix<T>& LI, matrix<T>& V_R, matrix<T>& V_I) {

	char charN = 'N';
	char charV = 'V';

    matrix<T> Mat_RearrA = arrangeElementsVerticallyForLAPACK(A);
	int rows = Mat_RearrA.rowCount();
	int cols = Mat_RearrA.colCount();
    matrix<T>Lr(LR.rowCount(), 1,0);
    matrix<T>Li(LR.rowCount(), 1,0);
    matrix<T>V_L(A.rowCount(),A.rowCount());
	int info;
	int lwork = 4*rows;
	double work[lwork];
	dgeev_(&charN, &charV, &rows, &Mat_RearrA(0,0), &rows, &Lr(0,0), &Li(0,0), &V_L(0,0), &rows, &V_L(0,0), &rows, 
	work, &lwork, &info);
	if (info != 0) {
		pout << "DGEEV returned info: " << info << "\n";
	}
    fillDiagonal(LR, Lr);
    fillDiagonal(LI, Li);
	V_L = transpose(V_L);
	int i=0;
	do {
		if((Li(i,0) == -Li(i+1,0)) && (Lr(i,0) == Lr(i+1,0))){
			for(unsigned int k=0; k<V_L.rowCount(); ++k){
				V_R(k, i)    =  V_L(k,i);
				V_I(k, i)    =  V_L(k,i+1);
				V_R(k, i+1)  =  V_L(k,i);
				V_I(k, i+1)  = -V_L(k,i+1);
			}
			i=i+2;
		} else {
			for(unsigned int k=0; k<V_L.rowCount(); ++k){
				V_R(k,i) = V_L(k,i);
				V_L(k,i) = 0.0;
			}		
			i=i+1;
		}
		
	} while(i<Lr.rowCount());
#ifdef CheckEigen
    pout << "Lr (real part of eigenValues):\n";
    pout << Lr << "\n";
    pout << "Li (imaginary part of eigenValues):\n";
    pout << Li << "\n";
    pout << "V_L (real part of eigenValues):\n";
    pout << V_L << "\n";
    pout << "Vr (real part of eigenValues):\n";
    pout << V_R << "\n";
    pout << "Vi (imag part of eigenValues):\n";
    pout << V_I << "\n";
    pout << "A*v (Real):\n";
    pout << A*V_R  << "\n";
    pout << "v*L (Real):\n";
    pout << V_R*LR - V_I*LI << "\n";
    pout << "A*v (Imag):\n";
    pout << A*V_I  << "\n";
    pout << "v*L (Imag):\n";
    pout << V_R*LI + V_I*LR << "\n";
#endif
	Lr.destroy();
	Li.destroy();
	V_L.destroy();
	return info;
}

extern "C" void zgeev_(char* jobvl, char* jobvr, int* n, std::complex<double>* A, int* lda, std::complex<double>* w, std::complex<double>* vl, int* ldvl, std::complex<double>* vr, 
int* ldvr, std::complex<double> work[], int* lwork, double rwork[], int* info);
template<typename CT, typename T>
int computeEigenValuesOfAComplexMatUsingLAPACK(matrix<CT>& A, matrix<CT>& L, matrix<CT>& VR) {

	char charN = 'N';
	char charV = 'V';

    matrix<CT> Mat_RearrA = arrangeElementsVerticallyForLAPACK(A);
	int rows = Mat_RearrA.rowCount();
	int cols = Mat_RearrA.colCount();
    matrix<CT>VR_feed(A.rowCount(),A.rowCount());
	int info;
	int lwork = 4*rows;
	CT work[lwork];
	T rwork[2*rows];
  matrix<CT>L_feed(L.rowCount(), 1,0);

	zgeev_(&charN, &charV, &rows, &Mat_RearrA(0,0), &rows, &L_feed(0,0), &VR_feed(0,0), &rows, &VR_feed(0,0), &rows, 
	work, &lwork, rwork, &info);
	if (info != 0) {
		pout << "DGEEV returned info: " << info << "\n";
	}
  fillDiagonal(L, L_feed);
	VR = transpose(VR_feed);
	int i=0;
#ifdef CheckEigen
    pout << "Eigenvalues:\n";
    pout << L << "\n";
    pout << "A*v:\n";
    pout << A*VR  << "\n";
    pout << "v*L:\n";
    pout << VR*L << "\n";
#endif
	Mat_RearrA.destroy();
	L_feed.destroy();
	VR_feed.destroy();
	return info;
}

extern "C" void zgesvd_(char* jobu, char* jobvt, int* m, int* n,
std::complex<double>* A, int* lda, double *S, std::complex<double>* U,
int* ldu, std::complex<double>* VT, int* ldvt, std::complex<double>* work,
int* lwork, double* rwork, int* info);

template<typename CT, typename T>
int computeComplexSVDUsingLAPACK(matrix<CT>& A, matrix<CT>& U, matrix<CT>& Sigma, matrix<CT>& VT) {
#ifdef CheckSVD
    pout << "Complex U matrix to take SVD of: \n";
	pout << A << "\n";
#endif

	matrix<CT> Mat_RearrA = arrangeElementsVerticallyForLAPACK(A);
#ifdef CheckSVD
    pout << "Rearranged matrix Mat_RearrA: \n";
	pout << Mat_RearrA << "\n";
#endif

	char computeSpecific = 'S'; //Min(m,n) Columns of U are computed
	char computeAll = 'A';
	int rows = Mat_RearrA.rowCount();
	int cols = Mat_RearrA.colCount();
	matrix<T>S(Sigma.rowCount(), 1,0);

	int info;
	int lwork = MAX(MAX(1,3*MIN(rows,cols) + MAX(rows,cols)),5*MIN(rows,cols));
	double rwork[5*MIN(rows,cols)];
	CT  work[lwork];
	zgesvd_(&computeAll, &computeAll, &rows, &cols, &Mat_RearrA(0,0), &rows, &S(0,0), &U(0,0), &rows, &VT(0,0), &cols, &work[0], &lwork, &rwork[0], &info);
	if (info != 0) {
		std::cout << "#ZGESVD returned info: " << info << "\n";
	}
	for(unsigned int i=0; i<S.rowCount(); ++i){
		
		if(S(i,0)<1e-20){
			S(i,0) = 0.0;
		}
		pout << S(i,0) << "\n";
	}

	U = transpose(U);
	VT = transpose(VT);
  fillDiagonalRealToComplex(Sigma, S);

#ifdef Print
	pout << "#Singular values of complex U matrix: \n";
	pout << S << "\n";
	pout << "U: \n";
	pout << U << "\n";
	pout << "VT: \n";
	pout << VT << "\n";
#endif


#ifdef CheckSVD
		pout << "#Check complex reconstructed U matrix from SVD: \n";
    matrix<CT> SVD = U*Sigma*VT;
    pout << SVD << "\n";
#endif

	S.destroy();
	return info;
}

/*
 * 	SVD test from page 332 Gilbert Strang's Linear Algebra
 * 	Run the following it in main
 *
	container<nvector<Type>> MyMat;
	MyMat.push_back(nvector<Type>(2, {-1, 0}));
	MyMat.push_back(nvector<Type>(2, {1, -1}));
	MyMat.push_back(nvector<Type>(2, {0, 1}));
	matrix<Type> Mat = construct_from_column_vectors(MyMat);
	matrix<Type> Mat_Rearr = arrangeElementsVerticallyForLAPACK(Mat);
	pout << "Mat_Rearr: \n";
	pout << Mat_Rearr << "\n";
	size_type R_rows = Mat.rowCount();
	size_type R_cols = Mat.colCount();
	matrix<Type> U_r(R_rows, R_rows);
	matrix<Type> V_T_r(R_cols, R_cols);
	matrix<Type> Sigma(MAX(R_rows, R_cols),1, 0 );
	computeSVDUsingLAPACK(Mat_Rearr, U_r, Sigma, V_T_r);

	*/


#endif //SVD
