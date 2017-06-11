#ifndef __EIGENVECTOR_H_
#define __EIGENVECTOR_H_

#include <vector>
#include <iostream>
#include <fstream>

class Eigenvector
{
public:
    Eigenvector(int64_t dim, int64_t index) : dim_(dim), index_(index) {}
    inline int64_t at(int64_t pos) const;
    inline int64_t operator[](int64_t pos) const;
    inline const std::vector<int64_t>& coefficients(void) const;
    inline bool operator==(const Eigenvector& vec) const;
    inline int64_t index(void) const;

private:
    int64_t dim_;
    int64_t index_;
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

bool Eigenvector::operator==(const Eigenvector& vec) const
{
    return vec.index_ == this->index_;
}

int64_t Eigenvector::index(void) const
{
    return this->index_;
}

#endif // __EIGENVECTOR_H_
