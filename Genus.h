#ifndef __GENUS_H_
#define __GENUS_H

#include <memory>
#include <iomanip>
#include <cassert>
#include <set>
#include <functional>
#include "QuadForm.h"
#include "Prime.h"
#include "NeighborIterator.h"

template<typename R, typename F>
class Genus
{
public:
    typedef std::shared_ptr<QuadForm<R,F>> QuadFormPtr;

    Genus(const QuadForm<R,F>& q);

    QuadFormPtr quad_form(void) const;

    std::shared_ptr<Prime<R,F>> smallest_good_prime(void) const;

    void print(void) const;

    size_t size(void) const;

private:
    std::set<QuadForm<R,F>> genusSet_;

    /* Computes the full genus. */
    void compute_genus(void);

    /* A vector of shared QuadForm pointers. An unordered list of the genus. */
    std::vector<QuadFormPtr> genusReps_;

    /* A pointer to the originating form. */
    QuadFormPtr q_;

    /* The discriminant. */
    R disc_;

    /* Flag indicating whether the genus has been computed in full. */
    bool computed_ = false;
};

template<typename R, typename F>
Genus<R,F>::Genus(const QuadForm<R,F>& q)
{
    this->q_ = QuadForm<R,F>::reduce(q, false);
    this->disc_ = this->q_->discriminant();
    this->compute_genus();
}

template<typename R, typename F>
std::shared_ptr<QuadForm<R,F>> Genus<R,F>::quad_form(void) const
{
    return this->q_;
}

template<typename R, typename F>
size_t Genus<R,F>::size(void) const
{
    return this->genusSet_.size();
}

template<typename R, typename F>
void Genus<R,F>::compute_genus(void)
{
    // Do nothing if the genus has already been computed.
    if (this->computed_) { return; }

    // Insert the initial quadratic form.
    this->genusSet_.insert(*this->q_);
    this->genusReps_.push_back(this->q_);

    // Determine the smallest good prime.
    std::shared_ptr<Prime<R,F>> p = this->smallest_good_prime();

    // Loop over all genus representatives as the genus is built.
    int64_t index = 0;
    while (index < this->genusReps_.size())
    {
        // The current genus representative.
        QuadFormPtr cur = this->genusReps_[index++];

        // The neighbor iterator.
        NeighborIterator<R,F> it(cur, p);

        // The first neighbor.
        QuadFormPtr pn = it.next_neighbor();

        // Loop over all p-neighbors.
        while (pn != nullptr)
        {
            // Reduce the p-neighbor.
            QuadFormPtr qq = QuadForm<R,F>::reduce(*pn);

            // Check whether this form is already in the genus.
            if (this->genusSet_.count(*qq) == 0)
            {
                this->genusSet_.insert(*qq);
                this->genusReps_.push_back(qq);
            }

            // Get the next p-neighbor.
            pn = it.next_neighbor();
        }
    }
}

template<typename R, typename F>
void Genus<R,F>::print(void) const
{
    for (auto& q : this->genusSet_)
    {
        std::cout << q << std::endl;
    }
}

#endif // __GENUS_H_
