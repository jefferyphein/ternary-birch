#ifndef __EIGENVECTOR_H_
#define __EIGENVECTOR_H_

#include <algorithm>
#include "SetCover.h"

template<typename R>
class Eigenvector
{
public:
    Eigenvector() = default;
    Eigenvector(const Eigenvector<R>& other) = default;
    Eigenvector(Eigenvector<R>&& other) = default;

    Eigenvector<R>& operator=(const Eigenvector<R>& other) = default;
    Eigenvector<R>& operator=(Eigenvector<R>&& other) = default;

    Eigenvector(std::vector<Z32>&& data, W64 conductor_index)
    {
        this->data_ = std::vector<Z32>(data);
        this->conductor_index_ = conductor_index;
    }

    const std::vector<Z32>& data(void) const
    {
        return this->data_;
    }

    size_t size(void) const
    {
        return this->data_.size();
    }

    W64 conductor_index(void) const
    {
        return this->conductor_index_;
    }

    void rep_index(size_t pos)
    {
        this->rep_index_ = pos;
    }

    size_t rep_index(void) const
    {
        return this->rep_index_;
    }

    Z32 operator[](size_t pos) const
    {
        return this->data_[pos];
    }

private:
    std::vector<Z32> data_;
    W64 conductor_index_;
    size_t rep_index_;
};

template<typename R>
class EigenvectorManager
{
friend class Genus<R>;

public:
    EigenvectorManager() = default;

    void add_eigenvector(const Eigenvector<R>& vector)
    {
        if (this->finalized)
        {
            throw std::logic_error("Cannot add eigenvectors once finalized.");
        }

        if (this->dimension > 0)
        {
            if (this->dimension != vector.size())
            {
                throw std::invalid_argument("Eigenvector dimensions must match.");
            }
        }
        else
        {
            this->dimension = vector.size();
        }

        this->eigenvectors.push_back(vector);
    }

    size_t size(void) const
    {
        return this->eigenvectors.size();
    }

    void finalize(void)
    {
        if (this->finalized)
        {
            throw std::logic_error("Cannot finalize again.");
        }

        size_t num_vecs = this->eigenvectors.size();

        // If there are no eigenvectors, then there's nothing to do.
        if (num_vecs == 0)
        {
            this->finalized = true;
            return;
        }

        this->conductors.reserve(num_vecs);

        // First, we need to determine which coordinates will allow us to most
        // efficiently compute eigenvalues. If there is a coordinate which is
        // nonzero in each eigenvector, we can compute all eigenvalues using
        // only the associated genus representative. Otherwise, we need to solve
        // a set covering problem to find a subset of coordinates with nonzero
        // values, and use multiple genus representatives.

        size_t num_words = (this->dimension + 63) / 64;
        std::vector<std::vector<W64>> covers;
        covers.reserve(this->size());

        // Construct the set covers from the eigenvectors.
        for (const Eigenvector<R>& eigenvector : this->eigenvectors)
        {
            const std::vector<Z32>& data = eigenvector.data();

            std::vector<W64> cover;
            cover.reserve(this->dimension);

            size_t pos = 0;
            W64 word = 0;
            W64 mask = 1;
            for (Z32 value : data)
            {
                if (value) word |= mask;
                mask <<= 1;
                ++pos;

                if (pos % 64 == 0)
                {
                    cover.push_back(word);
                    word = 0;
                    mask = 1;
                }
            }

            if (this->dimension % 64 != 0)
            {
                cover.push_back(word);
            }

            covers.push_back(cover);
        }

        // Find a set cover and sort the positions.
        SetCover cover(this->dimension, covers, SetCover::METHOD_KINDA_GREEDY);
        this->indices = cover.positions();
        std::sort(this->indices.begin(), this->indices.end());
        size_t num_indices = this->indices.size();

        // Store the position of each eigenvector associated to each index.
        this->position_lut.resize(num_indices);

        // For each eigenvector, set the genus representative to be used for
        // computing eigenvalues.
        for (size_t n=0; n<num_vecs; n++)
        {
            Eigenvector<R>& eigenvector = this->eigenvectors[n];
            for (size_t pos=0; pos<num_indices; pos++)
            {
                size_t index = this->indices[pos];
                if (eigenvector.data()[index])
                {
                    eigenvector.rep_index(index);
                    this->position_lut[pos].push_back(n);
                    break;
                }
            }
        }

        // Next, we stride the eigenvector coodinate data so that we can compute
        // eigenvalues in a cache friendly way. If we have the following
        // eigenvectors:
        //  ( 0, 1, 1, -1)
        //  ( 1, 0, 2, -2)
        //  ( 3, 0, 1,  0)
        //
        // This will get strided as
        //  (0, 1, 3, ..., 1, 0, 0, ..., 1, 2, 1, ..., -1, -2, 0, ...)
        // where "..." indicates some number of zeros so that the front of each
        // coordinate boundary is at the beginning of a cache line.

        // We assume a 64-byte cache line.
        this->stride = ((num_vecs + 15) / 16) * 16;

        // Allocate memory for the strided eigenvectors.
        this->strided_eigenvectors.resize(this->stride * this->dimension);

        // Interleave the eigenvectors.
        size_t offset = 0;
        for (const Eigenvector<R>& eigenvector : this->eigenvectors)
        {
            this->conductors.push_back(eigenvector.conductor_index());

            const std::vector<Z32> vec = eigenvector.data();
            for (size_t pos=0; pos<this->dimension; pos++)
            {
                this->strided_eigenvectors[pos * this->stride + offset] = vec[pos];
            }
            ++offset;
        }

        // Reduce all conductors using bitwise-or to determine which primitive
        // characters are needed for computing eigenvalues.
        this->conductor_primes = std::accumulate(
            this->conductors.begin(),
            this->conductors.end(), 0, std::bit_or<W64>());

        // Set the finalized flag.
        this->finalized = true;
    }

    const Eigenvector<R>& operator[](size_t index) const
    {
        return this->eigenvectors[index];
    }

private:
    bool finalized = false;
    size_t dimension = 0;
    size_t stride = 0;
    std::vector<Eigenvector<R>> eigenvectors;
    std::vector<Z32> strided_eigenvectors;
    std::vector<W64> conductors;
    std::vector<std::vector<Z64>> position_lut;
    std::vector<Z64> indices;
    W64 conductor_primes;
};

#endif // __EIGENVECTOR_H_
