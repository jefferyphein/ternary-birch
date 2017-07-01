#ifndef __HECKE_OPERATOR_H_
#define __HECKE_OPERATOR_H_

#include <iostream>
#include "SparseMatrix.h"
#include "Genus.h"

template<typename R, typename F>
class Genus;

template<typename R, typename F>
class HeckeOperator
{
public:
    HeckeOperator() = default;
    HeckeOperator(const Genus<R,F>& genus, const Character<R,F>& chi);
    void print(std::ostream& os) const;
    void update_row(int64_t rowNumber, const std::map<int64_t, int64_t>& theRow);
    void add_row(int64_t rowNumber, const std::map<int64_t, int64_t>& theRow);
    int64_t num_rows(void) const;
    void import(std::ifstream& is);

    const SparseMatrix& matrix(void) const;
    int64_t at(int64_t row, int64_t col) const;

private:
    int64_t dim_;
    SparseMatrix m_;
};

template<typename R, typename F>
HeckeOperator<R,F>::HeckeOperator(const Genus<R,F>& genus, const Character<R,F>& chi)
{
    this->dim_ = genus.dimension(chi);
    this->m_ = SparseMatrix(this->dim_, this->dim_);
}

template<typename R, typename F>
void HeckeOperator<R,F>::print(std::ostream& os) const
{
    os << this->m_.num_rows() << " " << this->dim_ << std::endl;
    os << this->m_;
}

template<typename R, typename F>
void HeckeOperator<R,F>::update_row(int64_t rowNumber, const std::map<int64_t, int64_t>& theRow)
{
    this->m_.update_row(rowNumber, theRow);
}

template<typename R, typename F>
void HeckeOperator<R,F>::add_row(int64_t rowNumber, const std::map<int64_t, int64_t>& theRow)
{
    this->m_.add_row(rowNumber, theRow);
}

template<typename R, typename F>
int64_t HeckeOperator<R,F>::num_rows(void) const
{
    return this->m_.num_rows();
}

template<typename R, typename F>
std::ostream& operator<<(std::ostream& os, const HeckeOperator<R,F>& hecke)
{
    hecke.print(os);
    return os;
}

template<typename R, typename F>
void HeckeOperator<R,F>::import(std::ifstream& is)
{
    int64_t rows;
    is >> rows >> this->dim_;
    for (int64_t n = 0; n < rows; n++)
    {
        int64_t row, cols;
        is >> row >> cols;

        std::map<int64_t, int64_t> theRow;

        for (int64_t m = 0; m < cols; m++)
        {
            int64_t col, value;
            is >> col >> value;
            theRow[col] = value;
        }
        this->add_row(row, theRow);
    }
}

template<typename R, typename F>
const SparseMatrix& HeckeOperator<R,F>::matrix(void) const
{
    return this->m_;
}

template<typename R, typename F>
int64_t HeckeOperator<R,F>::at(int64_t row, int64_t col) const
{
    return this->m_.at(row, col);
}

#endif // __HECKE_OPERATOR_H_
