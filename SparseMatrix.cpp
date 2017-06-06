#include <map>
#include <iostream>
#include <iomanip>
#include <gmpxx.h>
#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(int64_t rows, int64_t cols)
{
    this->rows_ = rows;
    this->cols_ = cols;
    this->data_ = std::map<int64_t, std::map<int64_t, int64_t>>();
}

void SparseMatrix::print(void) const
{
    for (auto& row : this->data_)
    {
        std::cout << row.first;
        for (auto& col : row.second)
        {
            std::cout << " " << col.first << " " << col.second;
        }
        std::cout << std::endl;
    }
}

void SparseMatrix::update_row(int64_t rowNumber,
                              const std::map<int64_t, int64_t>& theRow)
{
    if (this->data_.count(rowNumber) == 0)
    {
        this->data_[rowNumber] = theRow;
    }
    else
    {
        std::map<int64_t, int64_t>& storedRow = this->data_[rowNumber];
        for (auto& temp : theRow)
        {
            storedRow[temp.first] += temp.second;
        }
    }
}

void SparseMatrix::add_row(int64_t rowNumber,
                           const std::map<int64_t, int64_t>& theRow)
{
    this->data_[rowNumber] = theRow;
}
