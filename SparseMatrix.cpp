#include <iostream>
#include <gmpxx.h>
#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(int64_t rows, int64_t cols)
{
    this->rows_ = rows;
    this->cols_ = cols;
    this->data_ = std::map<int64_t, std::map<int64_t, int64_t>>();
}

void SparseMatrix::print(std::ostream& os) const
{
    for (auto& row : this->data_)
    {
        os << row.first << " " << row.second.size();
        for (auto& col : row.second)
        {
            os << " " << col.first << " " << col.second;
        }
        os << std::endl;
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

int64_t SparseMatrix::num_rows(void) const
{
    return this->data_.size();
}

int64_t SparseMatrix::at(int64_t row, int64_t col) const
{
    // Check for invalid rows.
    if (row < 0 || row >= this->rows_)
    {
        throw std::range_error("Invalid row.");
    }

    // Check for invalid columns.
    if (col < 0 || col >= this->cols_)
    {
        throw std::range_error("Invalid column.");
    }

    // If the specified row doesn't exist, the value is zero.
    auto it1 = this->data_.find(row);
    if (it1 == this->data_.end()) { return 0; }

    // Get the row.
    const std::map<int64_t, int64_t>& thisRow = it1->second;

    // If the specified column doesn't exist for this row, the value is zero.
    auto it2 = thisRow.find(col);
    if (it2 == thisRow.end()) { return 0; }

    // Return the value.
    return it2->second;
}

std::ostream& operator<<(std::ostream& os, const SparseMatrix& mat)
{
    mat.print(os);
    return os;
}
