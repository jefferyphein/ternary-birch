#ifndef __SPARSE_MATRIX_H_
#define __SPARSE_MATRIX_H_

#include <memory>
#include <map>

class SparseMatrix
{
public:
    typedef std::map<int64_t, std::map<int64_t, int64_t>> dataMap;

    SparseMatrix() = default;
    SparseMatrix(int64_t rows, int64_t cols);

    void print(std::ostream& os) const;
    void update_row(int64_t rowNumber, const std::map<int64_t, int64_t>& theRow);
    void add_row(int64_t rowNumber, const std::map<int64_t, int64_t>& theRow);
    int64_t num_rows(void) const;

    inline const dataMap& data(void) const;
    int64_t at(int64_t row, int64_t col) const;

private:
    int64_t rows_;
    int64_t cols_;
    dataMap data_;
};

const SparseMatrix::dataMap& SparseMatrix::data(void) const
{
    return this->data_;
}

std::ostream& operator<<(std::ostream& os, const SparseMatrix& mat);

#endif // __SPARSE_MATRIX_H_
