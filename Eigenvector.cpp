#include <iostream>
#include "Eigenvector.h"

std::istream& operator>>(std::istream& is, Eigenvector& vec)
{
    for (int64_t k = 0; k < vec.dim_; k++)
    {
        int64_t value;
        is >> value;
        vec.coefficients_.push_back(value);
    }
    return is;
}

std::ostream& operator<<(std::ostream& os, const Eigenvector& vec)
{
    for (int64_t k = 0; k < vec.dim_; k++)
    {
        os << vec.coefficients_[k] << " ";
    }
    return os;
}
