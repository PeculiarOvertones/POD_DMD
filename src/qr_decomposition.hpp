//
// Created by Makram Kamaleddine on 1/24/16.
//

#ifndef QR_DECOMPOSITION_QR_DECOMPOSITION_HPP
#define QR_DECOMPOSITION_QR_DECOMPOSITION_HPP

//#include <lapacke_utils.h>
#include "ops.hpp"
#include <complex>

template<typename T, typename InnerProd>
std::pair<matrix<T>, matrix<T>> qr_decomposition(const container<nvector<T>> basis, InnerProd inner_prod)
{
    // alias the size_type
    using size_type = typename matrix<T>::size_type;


    // perform gram-schmidt orthogonalization
    container<nvector<T>> orthonormal_basis = orthonormalize(basis, inner_prod);
//    cout << "Orthonormalized Basis: \n";
//    for (size_type i = orthonormal.begin(); i != orthonormal.end(); ++i) {
//        cout << *i << "\n";
//    }


    // transform the basis back into a matrix
    matrix<T> q = construct_from_column_vectors(orthonormal_basis);
//	pout << "size of q: " << q.rowCount() << "  "<< q.colCount() << "\n";

    matrix<T> r(q.colCount(), q.colCount());
//	pout << "size of r: " << r.rowCount() << "  "<< r.colCount() << "\n";
    for (size_type j = 0; j < r.colCount(); ++j) {
        for (size_type i = 0; i <= std::min(j,(r.rowCount()-static_cast<size_type>(1))); ++i) {
//			std::cout << "std::min(j,r.rowCount()): " << std::min(j,(r.rowCount()-static_cast<size_type>(1))) << "\n";
            r(i, j) = inner_product(orthonormal_basis[i], basis[j], inner_prod);
        }
    }

//Free container<nvector<T>>
    for (size_type i = 0; i < orthonormal_basis.size(); ++i) {
        orthonormal_basis[i].destroy();
    }
    orthonormal_basis.clear();
    container<nvector<T>>(orthonormal_basis).swap(orthonormal_basis);

#ifdef CheckQR
    pout << "Q*R: \n";
    pout << q*r << "\n";
#endif

    return std::make_pair(q, r);
}

template<typename T, typename InnerProd>
std::pair<matrix<T>, matrix<T>> qr_decomposition(const matrix<T>& mat, InnerProd inner_prod)
{
#ifdef CheckQR
    pout << "Matrix to QR factorize: \n";
    pout << mat << "\n";
#endif

    // alias the size_type
    using size_type = typename matrix<T>::size_type;

    // get the columns of the input matrix
    container<container<T>> col_collection = mat.columnCollection();

    // transform it into a vector of nvector objects
    container<nvector<T>> basis;
    std::transform(col_collection.begin(), col_collection.end(), std::back_inserter(basis), [](container<T> ctr) {
        return nvector<T>::from_container(ctr);
    });

    // perform gram-schmidt orthogonalization
    container<nvector<T>> orthonormal_basis = improved_orthonormalize(basis, inner_prod);
//    pout << "Orthonormalized Basis: \n";
//    for (size_type i = orthonormal_basis.begin(); i != orthonormal_basis.end(); ++i) {
//        pout << i->norm() << "\n";
//    }


    // transform the basis back into a matrix
    matrix<T> q = construct_from_column_vectors(orthonormal_basis);
//	pout << "size of q: " << q.rowCount() << "  "<< q.colCount() << "\n";

    matrix<T> r(q.colCount(), q.colCount());
//	pout << "size of r: " << r.rowCount() << "  "<< r.colCount() << "\n";
    for (size_type j = 0; j < r.colCount(); ++j) {
        for (size_type i = 0; i <= std::min(j,(r.rowCount()-static_cast<size_type>(1))); ++i) {
            r(i, j) = inner_product(orthonormal_basis[i], basis[j], inner_prod);
        }
    }

//Free container<container<T>>
    for (size_type i = 0; i < col_collection.size(); ++i) {
        col_collection[i].clear();
    }
    col_collection.clear();
    container<container<T>>(col_collection).swap(col_collection);

//Free container<nvector<T>>
    for (size_type i = 0; i < orthonormal_basis.size(); ++i) {
        orthonormal_basis[i].destroy();
    }
    orthonormal_basis.clear();
    container<nvector<T>>(orthonormal_basis).swap(orthonormal_basis);


#ifdef CheckQR
    pout << "Q*R: \n";
    pout << q*r << "\n";
#endif


    return std::make_pair(q, r);
}

extern "C" void zgeqrf_(int* rows, int* cols, std::complex<double>* A, int* LDA, std::complex<double>* TAU, std::complex<double>* WORK, int* LWORK, int* INFO);
extern "C" void zungqr_(int* rows, int* cols, int* K, std::complex<double>* A, int* LDA, const std::complex<double>*TAU, std::complex<double>* WORK, int* LWORK, int* INFO);


template<typename CT>
std::pair<matrix<CT>, matrix<CT>>  computeQRFactorizationUsingLAPACKOfAComplexMatrix(matrix<CT>& A) {
  using size_type = typename matrix<CT>::size_type;

  matrix<CT> Mat_RearrA = arrangeElementsVerticallyForLAPACK(A);
  int rows = Mat_RearrA.rowCount();
  int cols = Mat_RearrA.colCount();
  CT TAU[MIN(rows,cols)];
#ifdef CheckQR
    pout << "Matrix to QR factorize: \n";
    pout << A << "\n";
    pout << "Mat_RearrA: \n";
    pout << Mat_RearrA << "\n";
#endif

  int info;
  int lwork = rows * cols; 
  CT  work[MAX(1,lwork)];

	zgeqrf_(&rows, &cols, &Mat_RearrA(0,0), &rows, &TAU[0], &work[0], &lwork, &info);
  if (info != 0) {
    std::cout << "#ZGEQRF returned info: " << info << "\n";
  }
		matrix<CT> RowMajOpArray = arrangeElementsInRowMajorForCpp(Mat_RearrA);
#ifdef CheckQR
    pout << "Mat_RearrA after zgeqrf: \n";
    pout << Mat_RearrA << "\n";
    pout << "RowMajor Mat_RearrA after zgeqrf: \n";
    pout << RowMajOpArray << "\n";
#endif

	matrix<CT> R(cols,cols);
	for(size_type j = 0; j < cols; ++j) 
	{
		for(size_type k = 0; k <= j; ++k) 
		{
			R(k,j) = RowMajOpArray(k, j);
		}
	}

    //Note: rows >=cols>= 2ndparameter <cols in this case> >=0
	zungqr_(&rows, &cols, &cols, &Mat_RearrA(0,0), &rows, &TAU[0], &work[0], &lwork, &info);
  if (info != 0) {
    std::cout << "#ZUNGQR returned info: " << info << "\n";
  }

	matrix<CT> Q = arrangeElementsInRowMajorForCpp(Mat_RearrA);

#ifdef CheckQR
    pout << "R: \n";
    pout << R << "\n";
    pout << "After zungqr: \n";
    pout << Mat_RearrA << "\n";
    pout << "Q*R: \n";
    pout << Q*R << "\n";
#endif

	Mat_RearrA.destroy();

  return std::make_pair(Q,R);
}


#endif //QR_DECOMPOSITION_QR_DECOMPOSITION_HPP
