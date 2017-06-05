#ifndef __SPARSE_MATRIX_H_
#define __SPARSE_MATRIX_H_

#include <memory>

class SparseMatrix
{
public:
    SparseMatrix() = default;
    SparseMatrix(const mpz_class& rows, const mpz_class& cols);

    void update(const mpz_class& row, const mpz_class& col, const mpz_class& delta);
    void print(void) const;
    void merge(const SparseMatrix& mat);

private:
    mpz_class rows_;
    mpz_class cols_;
    std::shared_ptr<std::map<mpz_class, std::map<mpz_class, mpz_class>>> data_;
};

#endif // __SPARSE_MATRIX_H_
