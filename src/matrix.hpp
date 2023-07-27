#ifndef matrix_hpp
#define matrix_hpp

#include <algorithm>
#include <iomanip>
#include <complex>

#include "common.hpp"

// fwd declaration: in order to use enable_if
//template<typename, typename = void> struct matrix;

// use matrices only for types that are arithmetic (i.e number types)
template<typename T>
//struct matrix<T,
//        typename std::enable_if<std::is_arithmetic<T>::value || std::is_same<T, std::complex<T>>::value>::type>
//struct matrix<T,
//        typename std::enable_if<std::is_arithmetic<T>::value>::type>
struct matrix
{
    typedef unsigned long size_type;

    matrix(size_type N, size_type M)
    : vec(N * M),
      rows(N),
      columns(M)
    {

    }

    matrix(size_type N, size_type M, T filler)
    : vec(N * M),
      rows(N),
      columns(M)
    {
        std::fill(vec.begin(), vec.end(), filler);
    }

    matrix(const matrix& other)
    : vec(other.vec),
      rows(other.rows),
      columns(other.columns)
    {

    }

    matrix(size_type N, size_type M, std::initializer_list<T> lst)
    : vec(N * M),
      rows(N),
      columns(M)
    {
        // only copy if they're the same size as the vector
        if (lst.size() == N * M) {
            std::copy(lst.begin(), lst.end(), vec.begin());
        }
    }

    matrix(size_type N, size_type M, container<T> vector)
    : vec(N * M),
      rows(N),
      columns(M)
    {
        if (vector.size() == N * M) {
            std::copy(vector.begin(), vector.end(), vec.begin());
        }
    }

    matrix& operator=(const matrix& other)
    {
        vec = other.vec;
        rows = other.rows;
        columns = other.columns;

        return *this;
    }

    T& operator()(size_type i, size_type j)
    {
        return vec[i * columns + j];
    }

    const T& operator()(size_type i, size_type j) const
    {
        return vec[i * columns + j];
    }

    container<T> operator()(std::string colon, size_type j)
    {
		container<T> result;
		if(colon==":") {
			result = getColumn(j);
		}
        return result;
    }
    const container<T> operator()(std::string colon, size_type j) const
    {
		container<T> result;
		if(colon==":") {
			result = getColumn(j);
		}
        return result;
    }

    size_type rowCount() const
    {
        return rows;
    }

    size_type colCount() const
    {
        return columns;
    }

    // data interface: return internal data
    container<T> data() const
    {
        return vec;
    }

	void destroy()
	{
		vec.clear();
		container<T>(vec).swap(vec);
		rows=0;
		columns=0;
	}

    /*
     * Extract the columns from this matrix into a
     * vector of vectors so that they can be used in other
     * routines (such as orthogonalization).
     */
    container<container<T>> columnCollection() const
    {
        container<container<T>> result;

        for (size_type i = 0; i < columns; ++i)
        {
            result.push_back(getColumn(i));
        }

        return result;
    }

    /*
     * Extract the rows from this matrix into a vector of
     * vectors so that they can be used in other routines
     * (such as orthogonalization).
     */
    container<container<T>> rowCollection() const
    {
        container<container<T>> result;

        for (size_type i = 0; i < rows; ++i)
        {
            result.push_back(getRow(i));
        }

        return result;
    }
 

private:
    container<T> vec;
    size_type rows;
    size_type columns;

    container<T> getRow(size_type rowIndex) const
    {
        container<T> result;

        for (size_type i = 0; i < columns; ++i)
        {
            result.push_back(this->operator()(rowIndex, i));
        }

        return result;
    }

    container<T> getColumn(size_type columnIndex) const
    {
        container<T> result;

        for (size_type i = 0; i < rows; ++i)
        {
            result.push_back(this->operator()(i, columnIndex));
        }

        return result;
    }
};

template<typename T>
std::ostream& operator<<(std::ostream& stream, const matrix<T>& mat)
{
    using size_type = typename matrix<T>::size_type;
    stream << "[\n";
    for (size_type i = 0; i < mat.rowCount(); ++i)
    {
        for (size_type j = 0; j < mat.colCount(); ++j)
        {
            if (j != mat.colCount() - 1) {
                stream << std::setw(10) << mat(i, j) << ", ";
            } else {
                stream << mat(i, j);
            }
        }
        if (i != mat.rowCount() - 1) {
            stream << ";" << std::endl;
        }
    }
    stream << "\n]";
    return stream;
}

template<typename T>
matrix<T> operator*(const matrix<T>& A, const matrix<T>& B)
{
#ifdef DEBUG
    
    if (A.colCount() != B.rowCount()) {
        std::cout << A.colCount() << "  " << B.rowCount() << "\n";
        throw std::invalid_argument("Matrices cannot be multiplied. Check dimensions!");
    }
#endif

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
    return result;
}

template<typename T>
matrix<std::complex<T>> operator*(const matrix<std::complex<T>>& A, const matrix<T>& B)
{
#ifdef DEBUG
    
    if (A.colCount() != B.rowCount()) {
        std::cout << A.colCount() << "  " << B.rowCount() << "\n";
        throw std::invalid_argument("Matrices cannot be multiplied. Check dimensions!");
    }
#endif

    using size_type = typename matrix<T>::size_type ;
    matrix<std::complex<T>> result(A.rowCount(), B.colCount());
    size_type AColumns = A.colCount();
    size_type ARows =  A.rowCount(); 
    for (size_type i = 0; i < ARows; ++i) {
        for (size_type j = 0; j < B.colCount(); ++j) {
            for (size_type k = 0; k < AColumns; ++k) {
								std::complex<T> val(std::real(A(i,k))*B(k,j), std::imag(A(i,k))*B(k,j));
                result(i,j) += val;
            }
        }
    }
    return result;
}

template<typename T>
matrix<std::complex<T>> operator*(const matrix<T>& A, const matrix<std::complex<T>>& B)
{
#ifdef DEBUG
    
    if (A.colCount() != B.rowCount()) {
        std::cout << A.colCount() << "  " << B.rowCount() << "\n";
        throw std::invalid_argument("Matrices cannot be multiplied. Check dimensions!");
    }
#endif

    using size_type = typename matrix<T>::size_type ;
    matrix<std::complex<T>> result(A.rowCount(), B.colCount());
    size_type AColumns = A.colCount();
    size_type ARows =  A.rowCount(); 
    for (size_type i = 0; i < ARows; ++i) {
        for (size_type j = 0; j < B.colCount(); ++j) {
            for (size_type k = 0; k < AColumns; ++k) {
								std::complex<T> val(A(i,k)*std::real(B(k,j)), A(i,k)*std::imag(B(k,j)));
                result(i,j) += val;
            }
        }
    }
    return result;
}

template<typename T>
container<T> operator*(const matrix<T>& A, const container<T>& B)
{
#ifdef DEBUG
    
    if (A.colCount() != B.rowCount()) {
        std::cout << A.colCount() << "  " << B.rowCount() << "\n";
        throw std::invalid_argument("Matrix and vector cannot be multiplied. Check dimensions!");
    }
#endif

    using size_type = typename matrix<T>::size_type ;
    container<T> result(A.rowCount());
    size_type AColumns = A.colCount();
    size_type ARows =  A.rowCount(); 
    for (size_type i = 0; i < ARows; ++i) {
        for (size_type k = 0; k < AColumns; ++k) {
            result[i] += A(i,k) * B[k];
        }
    }
    return result;
}

template<typename T>
matrix<T> operator+(const matrix<T>& A, const matrix<T>& B)
{
#ifdef DEBUG
    
    if (A.colCount() != B.rowCount()) {
        std::cout << A.colCount() << "  " << B.rowCount() << "\n";
        throw std::invalid_argument("Matrices cannot be multiplied. Check dimensions!");
    }
#endif

    using size_type = typename matrix<T>::size_type ;
    matrix<T> result(A.rowCount(), B.colCount());
    size_type AColumns = A.colCount();
    size_type ARows =  A.rowCount(); 
    for (size_type i = 0; i < ARows; ++i) {
        for (size_type j = 0; j < B.colCount(); ++j) {
                result(i,j) = A(i,j) +  B(i,j);
        }
    }
    return result;
}

template<typename T>
matrix<T> operator-(const matrix<T>& A, const matrix<T>& B)
{
#ifdef DEBUG
    
    if (A.colCount() != B.rowCount()) {
        std::cout << A.colCount() << "  " << B.rowCount() << "\n";
        throw std::invalid_argument("Matrices cannot be multiplied. Check dimensions!");
    }
#endif

    using size_type = typename matrix<T>::size_type ;
    matrix<T> result(A.rowCount(), B.colCount());
    size_type AColumns = A.colCount();
    size_type ARows =  A.rowCount(); 
    for (size_type i = 0; i < ARows; ++i) {
        for (size_type j = 0; j < B.colCount(); ++j) {
                result(i,j) = A(i,j) -  B(i,j);
        }
    }
    return result;
}

#endif /* matrix_hpp */
