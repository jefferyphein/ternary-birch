#ifndef __HECKE_OPERATOR_H_
#define __HECKE_OPERATOR_H_

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
    void print(void) const;
    void update_row(int64_t rowNumber, const std::map<int64_t, int64_t>& theRow);
    void add_row(int64_t rowNumber, const std::map<int64_t, int64_t>& theRow);

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
void HeckeOperator<R,F>::print(void) const
{
    std::cout << this->dim_ << " " << this->dim_ << std::endl;
    this->m_.print();
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

#endif // __HECKE_OPERATOR_H_
