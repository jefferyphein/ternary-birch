#include <map>
#include <iostream>
#include <iomanip>
#include <gmpxx.h>
#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(const mpz_class& rows, const mpz_class& cols)
{
    this->rows_ = rows;
    this->cols_ = cols;
    this->data_ = std::make_shared<std::map<mpz_class, std::map<mpz_class, mpz_class>>>();
}

void SparseMatrix::update(const mpz_class& row, const mpz_class& col, const mpz_class& delta)
{
    if (this->data_->count(row) > 0)
    {
        if ((*this->data_)[row].count(col) > 0)
        {
            (*this->data_)[row][col] += delta;
        }
        else
        {
            (*this->data_)[row][col] = delta;
        }
    }
    else
    {
        std::map<mpz_class, mpz_class> temp;
        (*this->data_)[row] = temp;
        (*this->data_)[row][col] = delta;
    } 
}

void SparseMatrix::print(void) const
{
    for (auto& row : *this->data_)
    {
        std::cout << row.first;
        for (auto& col : row.second)
        {
            std::cout << " " << col.first << " " << col.second;
        }
        std::cout << std::endl;
    }
}

void SparseMatrix::merge(const SparseMatrix& mat)
{
    for (auto& row : *mat.data_)
    {
        for (auto& col : row.second)
        {
            this->update(row.first, col.first, col.second);
        }
    }
}
