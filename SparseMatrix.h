#ifndef __SPARSE_MATRIX_H_
#define __SPARSE_MATRIX_H_

#include <memory>
#include <unordered_map>

class SparseMatrix
{
public:
    SparseMatrix() = default;
    SparseMatrix(int64_t rows, int64_t cols);

    void print(void) const;
    void update_row(int64_t rowNumber, const std::map<int64_t, int64_t>& theRow);
    void add_row(int64_t rowNumber, const std::map<int64_t, int64_t>& theRow);
    int64_t num_rows(void) const;

private:
    int64_t rows_;
    int64_t cols_;
    std::map<int64_t, std::map<int64_t, int64_t>> data_;
};

#endif // __SPARSE_MATRIX_H_
