//
// Created by Makram Kamaleddine on 1/21/16.
//

#ifndef QR_DECOMPOSITION_COMMON_HPP
#define QR_DECOMPOSITION_COMMON_HPP

#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <complex>

// alias template for a container type
template<typename T>
using container = std::vector<T>;

std::ofstream pout;


template<typename T>
std::ostream& operator<<(std::ostream& stream, const container<T>& nvec)
{
    stream << "[";
//	std::vector<T>::const_iterator it;
    for (auto it = nvec.cbegin(); it != nvec.cend(); ++it)
    {
        if (it != nvec.cend() - 1) {
            stream << " "<< *it << "; ";
        } else {
            stream << " "<< *it;
        }

    }
    stream << "]";
    return stream;
}

template<typename T>
T l2_norm(container<T> const& nvec) {
    T accum = 0.;
    for (auto it = nvec.begin(); it != nvec.end(); ++it)
    {
        accum += pow(*it,2);
    }
    return sqrt(accum);
}
#endif //QR_DECOMPOSITION_COMMON_HPP
