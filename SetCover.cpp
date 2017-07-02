#include <algorithm>
#include <iostream>
#include <memory>
#include "SetCover.h"

SetCover::SetCover(int64_t total, const std::vector<std::vector<uint64_t>>& vectors, int64_t method)
{
    this->total_ = total;

    int64_t numVectors = vectors.size();
    if (numVectors == 0)
    {
        throw std::runtime_error("No vectors provided.");
    }

    int64_t numLimbs = vectors[0].size();

    // First we attempt to identify a single column with all 1's, and so let's
    // "and" all the limbs together.
    std::vector<uint64_t> anded(numLimbs, ~0L);
    for (const std::vector<uint64_t>& vec : vectors)
    {
        int64_t n = 0;
        for (uint64_t value : vec)
        {
            anded[n++] &= value;
        }
    }

    // Now let's look for a nonzero limb, and identify the position where a
    // 1 bit occurs.
    int64_t n = 0;
    for (uint64_t value : anded)
    {
        if (value != 0)
        {
            int64_t pos = 0;
            while ((value & 1) == 0)
            {
                value = (value >> 1);
                ++pos;
            }
            this->positions_.push_back(64 * n + pos);
            return;
        }
        ++n;
    }

    /* We may now assume that there is no single column of bits which are all
     * set to 1. */

    switch (method)
    {
        case SetCover::METHOD_GREEDY:
            this->greedy(vectors);
            break;
        case SetCover::METHOD_KINDA_GREEDY:
            this->kinda_greedy(vectors);
            break;
        case SetCover::METHOD_BRUTE_FORCE:
            this->brute_force(vectors);
            break;
        default:
            throw std::runtime_error("Unrecognized set cover method.");
    }
}

void SetCover::greedy(const std::vector<std::vector<uint64_t>>& vectors)
{
    // A flag indicating whether we have reached the end of our loops.
    bool done = false;

    // The number of limbs to consider, as well as the number of vectors.
    int64_t numLimbs = vectors[0].size();
    int64_t numVectors = vectors.size();

    // Initialize the pivot index and weight.
    int64_t pivotIndex = -1;
    int64_t pivotWeight = 0;

    // Loop over all limbs and all 8-bit chunks of each limb.
    for (int64_t n = 0; !done && n < numLimbs; n++)
    {
        for (int64_t b = 0; !done && b < 8; b++)
        {
            // Determine column count for this byte position.
            int64_t counts[256] = {0};
            for (const std::vector<uint64_t>& vec : vectors)
            {
                ++counts[(vec[n] >> (8*b)) & 0xff];
            }

            // Determine column totals for this byte.
            int64_t totals[8] = {0};
            for (int64_t k = 0; k < 256; k++)
            {
                if (counts[k] == 0) { continue; }

                // Update column totals if counts[k] is nonzero.
                if (k & 0x01) { totals[0] += counts[k]; }
                if (k & 0x02) { totals[1] += counts[k]; }
                if (k & 0x04) { totals[2] += counts[k]; }
                if (k & 0x08) { totals[3] += counts[k]; }
                if (k & 0x10) { totals[4] += counts[k]; }
                if (k & 0x20) { totals[5] += counts[k]; }
                if (k & 0x40) { totals[6] += counts[k]; }
                if (k & 0x80) { totals[7] += counts[k]; }
            }

            // Loop over totals, updating best pivot index and weight.
            for (int64_t k = 0; !done && k < 8; k++)
            {
                // Determine the absolute position of this column and flag that
                // we're done if it equals or exceeds the maximum number of columns
                // to be considered.
                int64_t pos = 64*n + 8*b + k;
                if (pos >= this->total_)
                {
                    done = true;
                    continue;
                }

                // If we found a column with a 1 in each bit position, we're done.
                // Push the position onto the list, and return.
                if (totals[k] == numVectors)
                {
                    this->positions_.push_back(pos);
                    return;
                }

                // Update index and weight, if necessary.
                if (totals[k] > pivotWeight)
                {
                    pivotWeight = totals[k];
                    pivotIndex = pos;
                }
            }
        }
    }

    // Push the best pivot index onto the list of positions.
    this->positions_.push_back(pivotIndex);

    // Determine the subset of vectors which have zero value at the pivotIndex
    // we just chose.
    auto subset = std::make_shared<std::vector<std::vector<uint64_t>>>();
    for (const std::vector<uint64_t>& vec : vectors)
    {
        if ((vec[pivotIndex / 64] & (1L << (pivotIndex % 64))) == 0)
        {
            subset->push_back(vec);
        }
    }

    // Repeat the greedy algorithm on the subset.
    this->greedy(*subset);
}

void SetCover::kinda_greedy(const std::vector<std::vector<uint64_t>>& vectors)
{
    std::vector<int64_t> bestPositions;

    for (int64_t n = 0; n < this->total_; n++)
    {
        // Form a subset of vector whose n-th bit is zero.
        auto subset = std::make_shared<std::vector<std::vector<uint64_t>>>();
        for (const std::vector<uint64_t>& vec : vectors)
        {
            if ((vec[n / 64] & (1L << (n % 64))) == 0)
            {
                subset->push_back(vec);
            }
        }

        // Apply the greedy algorithm to the remaining subset of vectors.
        SetCover setcover(this->total_, *subset, SetCover::METHOD_GREEDY);
        subset.reset();

        // Update the best positions if we find a smaller covering set.
        size_t num = setcover.num_positions();
        if (bestPositions.empty() || num < bestPositions.size())
        {
            std::vector<int64_t> temp(setcover.positions());
            temp.push_back(n);
            bestPositions = std::move(temp);
        }
    }

    // Assign the best positions based on what we found.
    this->positions_ = std::move(bestPositions);
}

void SetCover::brute_force(const std::vector<std::vector<uint64_t>>&)
{
    throw std::runtime_error("Brute force method not yet implemented.");
}

bool SetCover::is_set_cover(const std::vector<std::vector<uint64_t>>& vectors, const std::vector<int64_t>& pivots) const
{
    // Determine whether each vector has a pivot at one of the specified pivots.
    for (const std::vector<uint64_t>& vec : vectors)
    {
        bool valid = false;
        for (int64_t pivot : pivots)
        {
            if (!valid && (vec[pivot / 64] & (1L << (pivot % 64))) != 0)
            {
                valid = true;
            }
        }

        if (!valid) { return false; }
    }

    return true;
}

const std::vector<int64_t>& SetCover::positions(void) const
{
    return this->positions_;
}

size_t SetCover::num_positions(void) const
{
    return this->positions_.size();
}
