#ifndef __EIGENVECTOR_H_
#define __EIGENVECTOR_H_

#include <vector>
#include <iostream>
#include <fstream>

class Eigenvector
{
public:
    Eigenvector(int64_t dim) : dim_(dim) {}
    inline int64_t at(int64_t pos) const;
    inline int64_t operator[](int64_t pos) const;
    inline const std::vector<int64_t>& coefficients(void) const;

private:
    int64_t dim_;
    std::vector<int64_t> coefficients_;
    friend std::istream& operator>>(std::istream& is, Eigenvector& vec);
    friend std::ostream& operator<<(std::ostream& os, const Eigenvector& vec);
};

int64_t Eigenvector::at(int64_t pos) const
{
    return this->coefficients_[pos];
}

int64_t Eigenvector::operator[](int64_t pos) const
{
    return this->coefficients_[pos];
}

const std::vector<int64_t>& Eigenvector::coefficients(void) const
{
    return this->coefficients_;
}

#endif // __EIGENVECTOR_H_
