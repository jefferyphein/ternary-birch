#ifndef __GENUS_H_
#define __GENUS_H

#include <memory>
#include <cassert>
#include "QuadForm.h"
#include "Prime.h"
#include "NeighborIterator.h"

static const bool RED = 0;
static const bool BLACK = 1;

/* Forward declaration. */
template<typename R, typename F>
class Genus;

template<typename R, typename F>
class GNODE
{
friend class Genus<R,F>;

public:
    GNODE(const Genus<R,F>& g, std::shared_ptr<QuadForm<R,F>> q);
    
private:
    std::shared_ptr<QuadForm<R,F>> q_;
    int64_t n_;
    std::weak_ptr<GNODE> parent_;
    std::shared_ptr<GNODE> left_;
    std::shared_ptr<GNODE> right_;
    std::shared_ptr<GNODE> next_;
    bool color_;
};

template<typename R, typename F>
class Genus
{
public:
    Genus(const QuadForm<R,F>& q);

    std::shared_ptr<QuadForm<R,F>> quad_form(void) const;

    std::shared_ptr<Prime<R,F>> smallest_good_prime(void) const;

    size_t size(void) const;

private:
    void compute_genus(void);

    std::shared_ptr<QuadForm<R,F>> q_;
    R disc_;
    std::vector<std::shared_ptr<QuadForm<R,F>>> genusReps_;
    bool computed_ = false;
    std::shared_ptr<GNODE<R,F>> root_;
    std::shared_ptr<GNODE<R,F>> head_;
    std::shared_ptr<GNODE<R,F>> tail_;
};

template<typename R, typename F>
GNODE<R,F>::GNODE(const Genus<R,F>& g, std::shared_ptr<QuadForm<R,F>> q)
{
    this->q_ = q;
    this->n_ = g.size()+1;
    this->color_= RED;
}

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
    return this->genusReps_.size();
}

template<typename R, typename F>
void Genus<R,F>::compute_genus(void)
{
    // Do nothing if the genus has already been computed.
    if (this->computed_) { return; }

    // Clear genus representatives.
    this->genusReps_.resize(0);

    // Create the initial node.
    std::shared_ptr<GNODE<R,F>> x = std::make_shared<GNODE<R,F>>(*this, this->q_);
    x->color_ = BLACK;
    this->root_ = x;
    this->head_ = x;
    this->tail_ = x;

    // Determine the smallest prime not dividing the discriminant.
    std::shared_ptr<Prime<R,F>> p = this->smallest_good_prime();

    // Starting node.
    std::shared_ptr<GNODE<R,F>> ptr = x;

    // Loop over all genus representatives.
    do
    {
        // The current genus representative quadratic form.
        std::shared_ptr<QuadForm<R,F>> cur = ptr->q_;

        // Instantiate a neighbor iterator.
        NeighborIterator<R,F> it(cur, p);

        // Get the first p-neighbor.
        std::shared_ptr<QuadForm<R,F>> pn = it.next_neighbor();

        // Loop over all p-neighbors.
        while (pn != nullptr)
        {
#ifdef DEBUG
            assert( pn->isometry()->is_isometry(*this->q_, *pn) );
#endif

            // Get the reduced form of this neighbor.
            std::shared_ptr<QuadForm<R,F>> qq = QuadForm<R,F>::reduce(*pn);

            // Multiply the p-neighbor isometry on the right by the reduction
            // isometry. This now represents an isometry from the original form
            // to the genus representative isometric to this p-neighbor.
            pn->isometry()->multiply_on_right_by(qq->isometry());

#ifdef DEBUG
            assert( pn->isometry()->is_isometry(*this->q_, *qq) );
#endif

            // Get the next p-neighbor.
            pn = it.next_neighbor();
        }

        ptr = ptr->next_;
    }
    while (ptr != nullptr);

    this->computed_ = true;
}

#endif // __GENUS_H_
