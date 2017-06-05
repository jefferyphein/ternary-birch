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
    void update(const mpz_class& row, const mpz_class& col, const mpz_class& value);
    void merge(const HeckeOperator<R,F>& hecke);
private:
    mpz_class dim_;
    std::shared_ptr<SparseMatrix> m_;
};

template<typename R, typename F>
HeckeOperator<R,F>::HeckeOperator(const Genus<R,F>& genus, const Character<R,F>& chi)
{
    this->dim_ = genus.dimension(chi);
    this->m_ = std::make_shared<SparseMatrix>(this->dim_, this->dim_);
}

template<typename R, typename F>
void HeckeOperator<R,F>::print(void) const
{
    std::cout << this->dim_ << " " << this->dim_ << std::endl;
    this->m_->print();
}

template<typename R, typename F>
void HeckeOperator<R,F>::update(const mpz_class& row, const mpz_class& col, const mpz_class& value)
{
    this->m_->update(row, col, value);
}

template<typename R, typename F>
void HeckeOperator<R,F>::merge(const HeckeOperator<R,F>& hecke)
{
    this->m_->merge(*hecke.m_);
}

#endif // __HECKE_OPERATOR_H_
