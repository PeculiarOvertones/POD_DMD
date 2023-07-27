//
// Created by Makram Kamaleddine on 1/21/16.
//

#ifndef QR_DECOMPOSITION_OPS_HPP
#define QR_DECOMPOSITION_OPS_HPP

#include <assert.h>

#include "nvector.hpp"
#include "processor.hpp"
#include <complex>


/*
 * An almost general inner product template function: if the
 * nvector has complex numbers inside it this can be the standard
 * inner product used if we're working with complex numbers.
 */
template<typename T, typename BinaryFunc>
typename std::enable_if<std::is_floating_point<T>::value, T>::type
inner_product(const nvector<T>& v1, const nvector<T>& v2, BinaryFunc func)
{
//    if (v1.size() != v2.size()) {
//	std::cout << v1.size() << "  " << v2.size() << "\n";
//        throw std::invalid_argument("Vectors must be of the same size");
//    }
    using size_type = typename nvector<T>::size_type ;

    // naive algorithm for an inner product
    T result = static_cast<T>(0);
    for (size_type i = 0; i < v1.size(); ++i)
    {
        result += func(v1[i], v2[i]);
    }

    return result;
}

/**
 *
 * Project the vector a onto the line
 * spanned by the vector e.
 * */
template<typename T, typename BinaryFunc>
nvector<typename std::enable_if<std::is_floating_point<T>::value, T>::type>
proj(const nvector<T>& e, const nvector<T>& a, BinaryFunc func)
{
	nvector<T> result = inner_product(e, a, func) * e;
	T normalize = inner_product(e, e, func);
	if(normalize > 1e-10){
		result = result/inner_product(e, e, func);
	}
    return result;
}

/*
 * Transpose the given matrix about the main diagonal
 * (i.e a standard transpose).
 */
template<typename T>
matrix<T> transpose(const matrix<T>& mat)
{
    using size_type = typename matrix<T>::size_type;

    matrix<T> result(mat.colCount(), mat.rowCount());
    for (size_type i = 0; i < result.rowCount(); ++i)
    {
        for (size_type j = 0; j < result.colCount(); ++j)
        {
            result(i, j) = mat(j, i);
        }
    }

    return result;
}

/**
 * Given a set of vectors that constitute the columns of a
 * matrix, construct the columns of the matrix.
 */
template<typename T>
matrix<T> construct_from_column_vectors(const container<nvector<T>>& basis)
{
#ifdef DEBUG
    // check if all the vectors are the same size
    size_type size = basis.front().size();
    for (size_type i = basis.begin(); i != basis.end(); ++i)
    {
        assert(size == i->size());
    }
#endif
    using size_type = typename matrix<T>::size_type;

    matrix<T> result(basis.front().size(), basis.size());

    for (size_type j = 0; j < basis.size(); ++j) {
        for (size_type i = 0; i < basis[j].size(); ++i) {
            result(i, j) = (basis[j])[i];
        }
    }

    return result;
}


template<typename T>
matrix<T> construct_from_column_vectors(const container<nvector<T>>& basis, unsigned int beginCol, unsigned int endCol)
{
    using size_type = typename matrix<T>::size_type;

    matrix<T> result(basis.front().size(), endCol-beginCol);
    for (size_type j = beginCol; j < endCol; ++j) {
        for (size_type i = 0; i < basis[j].size(); ++i) {
            result(i, j-beginCol) = (basis[j])[i];
        }
    }
    return result;
}

template<typename T>
matrix<T> time_shifted_construct_from_column_vectors(const container<nvector<T>>& basis, unsigned int beginCol, unsigned int endCol)
{
    using size_type = typename matrix<T>::size_type;

    matrix<T> result((TimeShift+1)*basis.front().size(), endCol-beginCol);
    for (size_type j = beginCol; j < endCol; ++j) {
			for(size_type t = 0; t < (TimeShift+1); ++t) {
				size_type rowSize=0;
				if(j+t > 0){
					rowSize=basis[j+t].size();
				} 
//				pout << "rowSize: " << rowSize << "\n";
        for (size_type i = 0; i < basis[j+t].size(); ++i) {
            result((i+t*rowSize), j-beginCol) = (basis[j+t])[i];
        }
			}
    }
    return result;
}

template<typename T>
matrix<T> construct_from_row_vectors(const container<nvector<T>>& basis)
{
    using size_type = typename matrix<T>::size_type;
#ifdef DEBUG
    // check if all the vectors are the same size
    size_type size = basis.front().size();
    for (auto i = basis.begin(); i != basis.end(); ++i)
    {
        assert(size == i->size());
	}
//	pout << "all vectors are of same size " << size  << " basis.size(): " << basis.size()<< "\n";
#endif

    matrix<T> result(basis.size(), basis.front().size());

    for (size_type i = 0; i < basis.size(); ++i) {
        for (size_type j = 0; j < basis[i].size(); ++j) {
            result(i, j) = (basis[i])[j];
        }
    }

    return result;
}

/**
 * Given a basis of vectors return the orthonormal basis
 * associated with it. This is essentially the Gram-Schmidt
 * procedure.
 */
template<typename T, typename InnerProd>
container<nvector<T>> orthonormalize(const container<nvector<T>>& basis, InnerProd inner_prod)
{
#ifdef DEBUG
    // check if all the vectors are the same size
    size_type size = basis.front().size();
    for (size_type i = basis.begin(); i != basis.end(); ++i)
    {
        assert(size == i->size());
    }
#endif

    using size_type = typename container<T>::size_type ;
    // orthonormal basis to return
    container<nvector<T>> result;
    // the size of a vector in this basis
    size_type size = basis.front().size();

    // do the first vector before the loop
    result.push_back(basis.front());
//    for (size_type i = 1; i < basis.size(); ++i) {
//        nvector<T> projection(size, 0);
//        for (size_type j = 0; j < i; ++j) {
//            projection = projection - proj(result[j], basis[i], inner_prod);
//        }
//        result.push_back(basis[i] + projection);
//    }
    for (size_type i = 1; i < basis.size(); ++i) {
        nvector<T> projection(size, 0);
        nvector<T> projection_old(size, 0);
        projection_old = basis[i] - proj(result[0], basis[i], inner_prod);
        for (size_type j = 1; j < i; ++j) {
			projection = projection_old - proj(result[j], projection_old, inner_prod);
            projection_old = projection;
        }
		projection = projection_old;
        result.push_back(projection);
    }

    // pseudo-map
    std::for_each(result.begin(), result.end(), [](nvector<T>& vec) {
        vec = (1.0 / vec.norm()) * vec;
    });

    return result;
}

template<typename T, typename InnerProd>
container<nvector<T>> improved_orthonormalize(const container<nvector<T>>& basis, InnerProd inner_prod)
{
#ifdef DEBUG
    // check if all the vectors are the same size
    size_type size = basis.front().size();
    for (size_type i = basis.begin(); i != basis.end(); ++i)
    {
        assert(size == i->size());
    }
#endif

    using size_type = typename container<T>::size_type ;
    // orthonormal basis to return
    container<nvector<T>> result;
    // the size of a vector in this basis
    size_type size = basis.front().size();

	
    // do the first vector before the loop
    result.push_back(basis.front());
    for (size_type i = 1; i < basis.size(); ++i) {
        nvector<T> projection(size, 0);
        nvector<T> projection_old(size, 0);
        projection_old = basis[i] - proj(result[0], basis[i], inner_prod);
        for (size_type j = 1; j < i; ++j) {
			projection = projection_old - proj(result[j], projection_old, inner_prod);
            projection_old = projection;
        }
		projection = projection_old;
        result.push_back(projection);
    }

    std::for_each(result.begin(), result.end(), [](nvector<T>& vec) {
		if(vec.norm() > 1e-10){
        	vec = (1.0 / vec.norm()) * vec;	
		}
    });
	T check = check_orthonormality(result, inner_prod);
	
    int gramSchmidtCounter=1;

	if(fabs(check) > 1e-6){
    	pout << "Orthonormality check failed: Orthonormal sum: " << check << " for gramschmidtcounter: " << gramSchmidtCounter  << "\n";
    	container<nvector<T>> result2;
		T newcheck = 0.0;
		do {	
			std::cout << "re-orthonormalizing for proc: " << MPI_PROC_ID << "\n";
			gramSchmidtCounter++;
    		result2.push_back(result.front());
    		for (size_type i = 1; i < result.size(); ++i) {
    		    nvector<T> projection(size, 0);
    		    nvector<T> projection_old(size, 0);
    		    projection_old = result[i] - proj(result2[0], result[i], inner_prod);
    		    for (size_type j = 1; j < i; ++j) {
					projection = projection_old - proj(result2[j], projection_old, inner_prod);
    		        projection_old = projection;
    		    }
				projection = projection_old;
    		    result2.push_back(projection);
    		}
    		std::for_each(result2.begin(), result2.end(), [](nvector<T>& vec) {
				if(vec.norm() > 1e-10){
    		    	vec = (1.0 / vec.norm()) * vec;
				}
    		});
    		newcheck = check_orthonormality(result2, inner_prod);
    		pout << "Orthonormal sum: " << newcheck << " for gramschmidtcounter: " << gramSchmidtCounter << "\n";
			//empty result
			
    		for (size_type i = 0; i < result.size(); ++i) {
				result[i].destroy();
    		}
			//empty result
	        result.clear();
    	    container<nvector<T>>(result).swap(result);

			result = result2;
	
			//empty result2
    		for (size_type i = 0; i < result2.size(); ++i) {
				result2[i].destroy();
    		}
	        result2.clear();
    	    container<nvector<T>>(result2).swap(result2);

		}while(fabs(newcheck) > 1e-6);
	}


    return result;
}

template<typename T, typename InnerProd>
T check_orthonormality(const container<nvector<T>>& basis, InnerProd inner_prod)
{
#ifdef DEBUG
    // check if all the vectors are the same size
    size_type size = basis.front().size();
    for (size_type i = basis.begin(); i != basis.end(); ++i)
    {
        assert(size == i->size());
    }
#endif

    using size_type = typename nvector<T>::size_type ;
    // do a very basic check: take the inner product of the first vec
    // with the rest
    T sum = 0;
    nvector<T> front = basis.front();
    for (size_type i = 1; i < basis.size(); ++i) {
//		pout << inner_product(front, basis[i], inner_prod) << "\n";
        sum += inner_product(front, basis[i], inner_prod);
    }

    return sum;
}

template<typename T>
matrix<T> matMatMul_QoldQnew_ForParallelQR(const matrix<T>& A, const matrix<T>& B)
{
#ifdef DEBUG

    if (A.colCount() != B.rowCount()/size_type(MPI_PROC_TOTAL_NUM)) {
    	std::cout << A.colCount() << "  " << B.rowCount() << "\n";
        throw std::invalid_argument("Matrices cannot be multiplied. Check dimensions!");
    }

#endif

//    pout <<"Multiplication A: \n";
//    pout << A << "\n";
//    pout << "& B: \n";
//    pout << B << "\n";
//
    using size_type = typename matrix<T>::size_type ;
    matrix<T> result(A.rowCount(), B.colCount());
	size_type AColumns = A.colCount();
	size_type ARows =  A.rowCount();
    for (size_type i = 0; i < ARows; ++i) {
        for (size_type j = 0; j < B.colCount(); ++j) {
        	for (size_type k = 0; k < AColumns; ++k) {
				result(i,j) += A(i,k) * B(k + AColumns*size_type(MPI_PROC_ID),j);
//				pout << "ijk: " << i << " " << j << " "<< k << "A(i,k): " << A(i,k) << " B(kshifted,j): " << B(k + AColumns*size_type(MPI_PROC_ID),j)  << " result(i,j): " << result(i,j) << "\n";
			}
        }
    }
//    pout << "A * B: \n";
//    pout << result << "\n";

    return result;
}

template<typename T>
matrix<T> matMatMul(const matrix<T>& A, const matrix<T>& B)
{
#ifdef DEBUG

    if (A.colCount() != B.rowCount()) {
    	std::cout << A.colCount() << "  " << B.rowCount() << "\n";
        throw std::invalid_argument("Matrices cannot be multiplied. Check dimensions!");
    }

#endif

//    pout <<"Multiplication A: \n";
//    pout << A << "\n";
//    pout << "& B: \n";
//    pout << B << "\n";

    using size_type = typename matrix<T>::size_type ;
    matrix<T> result(A.rowCount(), B.colCount());
	size_type AColumns = A.colCount();
	size_type ARows =  A.rowCount();
    for (size_type i = 0; i < ARows; ++i) {
        for (size_type j = 0; j < B.colCount(); ++j) {
        	for (size_type k = 0; k < AColumns; ++k) {
				result(i,j) += A(i,k) * B(k,j);
			}
        }
    }
//    pout << "A * B: \n";
//    pout << result << "\n";

    return result;
}


template<typename T>
matrix<T> arrangeElementsVerticallyForLAPACK(const matrix<T>& mat)
{
    using size_type = typename matrix<T>::size_type;

    matrix<T> result(mat.rowCount(), mat.colCount());
	size_type cols = mat.colCount();
	size_type r_I;
	size_type c_I;
	size_type counter=0;
    for (size_type j = 0; j < result.colCount(); ++j)
    {
    	for (size_type i = 0; i < result.rowCount(); ++i)
    	{
			r_I = size_type(counter/cols);
			c_I = counter%cols;
            result(r_I,c_I) = mat(i, j);
			counter++;
        }
    }

    return result;
}

template<typename T>
matrix<T> arrangeElementsInRowMajorForCpp(const matrix<T>& mat)
{
    using size_type = typename matrix<T>::size_type;

    matrix<T> result(mat.rowCount(), mat.colCount());
	size_type rows = mat.rowCount();
	size_type r_I;
	size_type c_I;
	size_type counter=0;
    for (size_type j = 0; j < result.rowCount(); ++j)
    {
    	for (size_type i = 0; i < result.colCount(); ++i)
    	{
				r_I = counter%rows;
				c_I = size_type(counter/rows);
      	result(r_I,c_I) = mat(j,i);
				counter++;
      }
    }

    return result;
}

template<typename T>
matrix<T> invertDiagonalVector(matrix<T>& mat)
{
    using size_type = typename matrix<T>::size_type ;
    matrix<T> result(mat.rowCount(), mat.rowCount());
    for (size_type i = 0; i < result.rowCount(); ++i)
    {
        for (size_type j = 0; j < result.colCount(); ++j)
        {
			if(i!=j){
        		result(i,j)= 0.0;
			}
			if(i==j){
        		result(i,j)= 1.0/mat(i,j);
			}
        }
    }
	return result;
}

template<typename T>
void fillDiagonal(matrix<T>& result, matrix<T>& D)
{
    using size_type = typename matrix<T>::size_type ;
    for (size_type i = 0; i < result.rowCount(); ++i)
    {
        for (size_type j = 0; j < result.colCount(); ++j)
        {
            if(i==j){
                result(i,j)= D(i,0);
            }
        }
    }
}

template<typename CT, typename T>
void fillDiagonalRealToComplex(matrix<CT>& result, matrix<T>& D)
{
    using size_type = typename matrix<T>::size_type ;
    for (size_type i = 0; i < result.rowCount(); ++i)
    {
        for (size_type j = 0; j < result.colCount(); ++j)
        {
            if(i==j){
								CT a(D(i,0),0.0);
                result(i,j)= a;
            }
        }
    }
}

template<typename T>
void reverseArrangeDiagonal(matrix<T>& D)
{
    using size_type = typename matrix<T>::size_type ;

    size_type D_rows = D.rowCount();
    matrix<T> Dr(D_rows, D_rows);
    for(size_type i=0; i<D_rows; ++i){
        Dr(D_rows-1-i,D_rows-1-i) = D(i,i);
    }
    for(size_type i=0; i<D_rows; ++i){
        D(i,i) = Dr(i,i);
    }
}

//template<typename T>
//void fillDiagonalComplex(matrix<std::complex<T>>& result, matrix<T>& Lr, matrix<T>& Li)
//{
//    using size_type = typename matrix<T>::size_type ;
//    for (size_type r = 0; r < result.rowCount(); ++r)
//    {
//        for (size_type j = 0; j < result.colCount(); ++j)
//        {
//            if(r==j){
//                (result(r,j)).real= Lr(r,0);
//                (result(r,j)).imag = Li(r,0);
//            }
//        }
//    }
//}

template<typename CT, typename T>
void constructComplexMatrixFromRealAndImagParts(matrix<CT>& C, matrix<T>& R, matrix<T>& I, int flag) 
{
    using size_type = typename matrix<T>::size_type ;
    for (size_type a = 0; a < C.rowCount(); ++a)
    {
      for (size_type b = 0; b < C.colCount(); ++b)
      {
				if(flag==2){
					CT val(R(a,b),I(a,b));
					C(a,b) = val;
				}
				if(flag==0){
					CT val(R(a,b),0.0);
					C(a,b) = val;
				}
				if(flag==1){
					CT val(0.0,I(a,b));
					C(a,b) = val;
				}
      }
    }
}

template<typename CT>
void computeComplexConjugate(matrix<CT>& result, const matrix<CT>& mat)
{
    using size_type = typename matrix<CT>::size_type;

    for (size_type i = 0; i < result.rowCount(); ++i)
    {
        for (size_type j = 0; j < result.colCount(); ++j)
        {
						CT a(mat(j,i).real(), -mat(j,i).imag());
            result(i, j) = a;
        }
    }
}

template<typename CT,typename T>
void computeSigmaPlus(matrix<CT>& result, const matrix<CT>& mat)
{
    using size_type = typename matrix<CT>::size_type;

		T eps = 1e-20;
    for (size_type i = 0; i < result.rowCount(); ++i)
    {
        for (size_type j = 0; j < result.colCount(); ++j)
        {
						if(abs(mat(j,i).real()) > eps) {
							CT a((1.0/mat(j,i).real()), 0.0);
            	result(i, j) = a;
						}
        }
    }
}
#endif //QR_DECOMPOSITION_OPS_HPP
