//
// Created by Makram Kamaleddine on 1/21/16.
//

#ifndef QR_DECOMPOSITION_NVECTOR_HPP
#define QR_DECOMPOSITION_NVECTOR_HPP

#include <numeric>
#include "matrix.hpp"

// forward declaration
//template<typename, typename = void> struct nvector;

// partial specialization of matrix in order to add some familiar operator[] syntax
// rather than using operator() for both matrices and vectors
template<typename T>
//struct nvector<T,
//        typename std::enable_if<std::is_arithmetic<T>::value || std::is_same<T, std::complex<T>>::value>::type> {
//struct nvector<T,
//        typename std::enable_if<std::is_arithmetic<T>::value>::type> {
struct nvector {

    typedef typename matrix<T>::size_type size_type;

    nvector(size_type N)
            : vec(N, 1)
    {

    }

    nvector(size_type N, T filler)
            : vec(N, 1, filler)
    {

    }

    nvector(const nvector& other)
            : vec(other.vec)
    {

    }

    nvector(const container<T> other)
            : vec(other.size(), 1, other)
    {

    }

    nvector(size_type N, std::initializer_list<T> lst)
            : vec(N, 1, lst)
    {

    }

    nvector& operator=(const nvector& other)
    {
        vec = other.vec;
        return *this;
    }

    T& operator[](size_type i)
    {
        return vec(i, 0);
    }

    const T& operator[](size_type i) const
    {
        return vec(i, 0);
    }

    size_type size() const
    {
        return vec.rowCount();
    }

    typename container<T>::const_iterator begin() const
    {
        return vec.data().cbegin();
    }

    typename container<T>::const_iterator end() const
    {
        return vec.data().cend();
    }

    container<T> data() const
    {
        return vec.data();
    }

    void destroy()
    {
        vec.destroy();
    }  

    T norm() const
    {
        container<T> data = vec.data();
        // map
        std::for_each(data.begin(), data.end(), [](T& value) {
            value = value * value;
        });
        // fold
        return std::sqrt(std::accumulate(data.begin(), data.end(), static_cast<T>(0)));
    }

    static nvector<T> from_container(const container<T>& ctr)
    {
        nvector<T> vec(ctr.size());

        for (size_type i = 0; i < vec.size(); ++i) {
            vec[i] = ctr[i];
        }

        return vec;
    }
private:
    matrix<T> vec;
};

template<typename T>
std::ostream& operator<<(std::ostream& stream, const nvector<T>& nvec)
{
    using size_type = typename nvector<T>::size_type;

    stream << "[";
    for (size_type i = 0; i < nvec.size(); ++i)
    {
        if (i != nvec.size() - 1) {
            stream << nvec[i] << "; ";
        } else {
            stream << nvec[i];
        }
    }
    stream << "]";
    return stream;
}

//template<typename T>
//std::ostream& operator>>(std::ostream& stream, const nvector<T>& nvec)
//{
//    using size_type = typename matrix<T>::size_type;
//    stream >> nvec;
//    return stream;
//}


// overload operator* for matrix-vector multiplication
template<typename T>
nvector<T> operator*(const matrix<T>& M, const nvector<T>& v)
{
    using size_type = typename nvector<T>::size_type;
    // dimension check
    if (M.colCount() != v.size()) {
        throw std::invalid_argument("Matrix and vector dimensions are not compatible");
    }

    nvector<T> result(v.size());
    for (size_type i = 0; i < M.rowCount(); ++i)
    {
        T sum = 0;
        // do an inner product
        for (size_type j = 0; j < M.colCount(); ++j)
        {
            sum += M(i, j) * v[j];
        }
        result[i] = sum;
    }

    return result;
}

// scalar-vector multiplication
template<typename T>
nvector<T> operator*(const T& scalar, const nvector<T>& v)
{
    using size_type = typename nvector<T>::size_type;
    nvector<T> result(v.size());
    for (size_type i = 0; i < result.size(); ++i)
    {
        result[i] = v[i] * scalar;
    }

    return result;
}

// scalar-vector division
template<typename T>
nvector<T> operator/(const nvector<T>& v, const T& scalar)
{
    using size_type = typename nvector<T>::size_type;
    nvector<T> result(v.size());
    for (size_type i = 0; i < result.size(); ++i)
    {
        result[i] = v[i] / scalar;
    }

    return result;
}

template<typename T>
nvector<T> operator-(const nvector<T>& v1, const nvector<T>& v2)
{
    using size_type = typename nvector<T>::size_type;
    // dim check
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vector dimensions are not equal");
    }

    nvector<T> result(v1.size());
    for (size_type i = 0; i < v1.size(); ++i)
    {
        result[i] = v1[i] - v2[i];
    }

    return result;
}

template<typename T>
nvector<T> operator+(const nvector<T>& v1, const nvector<T>& v2)
{
    using size_type = typename nvector<T>::size_type;
    // dim check
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vector dimensions are not equal");
    }

    nvector<T> result(v1.size());
    for (size_type i = 0; i < v1.size(); ++i)
    {
        result[i] = v1[i] + v2[i];
    }

    return result;
}

#endif //QR_DECOMPOSITION_NVECTOR_HPP
