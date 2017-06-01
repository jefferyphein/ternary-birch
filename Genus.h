#ifndef __GENUS_H_
#define __GENUS_H

#include <memory>
#include <iomanip>
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
    typedef std::shared_ptr<QuadForm<R,F>> QuadFormPtr;
    typedef std::shared_ptr<GNODE<R,F>> GNODEPtr;

    GNODE(const Genus<R,F>& g, QuadFormPtr q);
    
private:
    QuadFormPtr q_;
    int64_t n_;
    std::weak_ptr<GNODE> parent_;
    GNODEPtr left_;
    GNODEPtr right_;
    GNODEPtr next_;
    bool color_;
};

template<typename R, typename F>
class Genus
{
public:
    typedef std::shared_ptr<QuadForm<R,F>> QuadFormPtr;
    typedef std::shared_ptr<GNODE<R,F>> GNODEPtr;

    Genus(const QuadForm<R,F>& q);

    QuadFormPtr quad_form(void) const;

    std::shared_ptr<Prime<R,F>> smallest_good_prime(void) const;

    void print(void) const;

    size_t size(void) const;

private:
    void compute_genus(void);
    GNODEPtr find_representative(GNODEPtr x,
                                 QuadFormPtr q);
    GNODEPtr add_representative(QuadFormPtr q);
    void add_representative_fixup(GNODEPtr z);
    void left_rotate(GNODEPtr x);
    void right_rotate(GNODEPtr x);

    void print(GNODEPtr x, int64_t level=0) const;

    QuadFormPtr q_;
    R disc_;
    std::vector<QuadFormPtr> genusReps_;
    bool computed_ = false;
    GNODEPtr root_;
    GNODEPtr head_;
    GNODEPtr tail_;
};

template<typename R, typename F>
GNODE<R,F>::GNODE(const Genus<R,F>& g, QuadFormPtr q)
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

    // Add the defining quadratic form to the list of representatives.
    this->genusReps_.push_back(this->q_);

    // Create the initial node.
    GNODEPtr x = std::make_shared<GNODE<R,F>>(*this, this->q_);
    x->color_ = BLACK;
    this->root_ = x;
    this->head_ = x;
    this->tail_ = x;

    // Determine the smallest prime not dividing the discriminant.
    std::shared_ptr<Prime<R,F>> p = this->smallest_good_prime();

    // Starting node.
    GNODEPtr ptr = x;

    // Loop over all genus representatives.
    do
    {
        // The current genus representative quadratic form.
        QuadFormPtr cur = ptr->q_;

        // Instantiate a neighbor iterator.
        NeighborIterator<R,F> it(cur, p);

        // Get the first p-neighbor.
        QuadFormPtr pn = it.next_neighbor();

        // Loop over all p-neighbors.
        while (pn != nullptr)
        {
#ifdef DEBUG
            assert( pn->isometry()->is_isometry(*this->q_, *pn) );
#endif

            // Get the reduced form of this neighbor.
            QuadFormPtr qq = QuadForm<R,F>::reduce(*pn);

            // Look for this reduced form in the genus.
            auto x = this->find_representative(this->root_, qq);

            if (x == nullptr)
            {
                // Compose to obtain a global isometry between this genus
                // representative and the original form.
                qq->isometry()->multiply_on_left_by(pn->isometry());

#ifdef DEBUG
                assert( qq->isometry()->is_isometry(*this->q_, *qq) );
#endif

                // Add the reduced form to the genus.
                this->add_representative(qq);
            }

            // Get the next p-neighbor.
            pn = it.next_neighbor();
        }

        ptr = ptr->next_;
    }
    while (ptr != nullptr);

    this->computed_ = true;
}

template<typename R, typename F>
std::shared_ptr<GNODE<R,F>> Genus<R,F>::find_representative(GNODEPtr x, QuadFormPtr q)
{
    if (x == nullptr) { return nullptr; }

    int64_t c = q->compare(x->q_);

    if (c == 0) { return x; }
    else if (c == 1) { return this->find_representative(x->left_, q); }
    else { return this->find_representative(x->right_, q); }
}

template<typename R, typename F>
std::shared_ptr<GNODE<R,F>> Genus<R,F>::add_representative(QuadFormPtr q)
{
    GNODEPtr x = this->root_;
    GNODEPtr y = nullptr;
    GNODEPtr z = std::make_shared<GNODE<R,F>>(*this, q);

    // Traverse the BST to figure out where to insert the new node.
    while (x != nullptr)
    {
        y = x;
        if (z->q_->compare(x->q_) == 1) { x = x->left_; }
        else { x = x->right_; }
    }

    // Set the parent node of our new node as the leaf we visited last.
    z->parent_ = y;

    // Should the new node be inserted to the left or to the right of the leaf?
    if (y == nullptr) { this->root_ = z; }
    else if (z->q_->compare(y->q_) == 1) { y->left_ = z; }
    else { y->right_ = z; }

    // Fix the red-black tree.
    this->add_representative_fixup(z);

    // Add pointer to new representative to our vector.
    this->genusReps_.push_back(q);

    // Update the linked list of genus representatives.
    this->tail_->next_ = z;
    this->tail_ = z;

    return z;
}

template<typename R, typename F>
void Genus<R,F>::add_representative_fixup(GNODEPtr z)
{
    // Make a shared pointer out of the weak pointer.
    GNODEPtr zp = std::shared_ptr<GNODE<R,F>>(z->parent_);

    while (zp != nullptr && zp->color_ == RED)
    {
        // Make a shared pointer out of a weak pointer.
        GNODEPtr zpp(zp->parent_);

        if (zpp != nullptr && zp == zpp->left_)
        {
            GNODEPtr y = zpp->right_;
            if (y != nullptr && y->color_ == RED)
            {
                zp->color_ = BLACK;
                y->color_ = BLACK;
                zpp->color_ = RED;
                z = zpp;
            }
            else if (z == zp->right_)
            {
                z = zp;
                this->left_rotate(z);
            }
            else
            {
                zp->color_ = BLACK;
                zpp->color_ = RED;
                this->right_rotate(zpp);
            }
        }
        else
        {
            GNODEPtr y = zpp != nullptr ? zpp->left_ : nullptr;
            if (y != nullptr && y->color_ == RED)
            {
                zp->color_ = BLACK;
                y->color_ = BLACK;
                zpp->color_ = RED;
                z = zpp;
            }
            else if (z == zp->left_)
            {
                z = zp;
                this->right_rotate(z);
            }
            else
            {
                zp->color_ = BLACK;
                zpp->color_ = RED;
                this->left_rotate(zpp);
            }
        }

        if (z->parent_.use_count() == 0) { zp = nullptr; }
        else { zp = std::shared_ptr<GNODE<R,F>>(z->parent_); }
    }

    this->root_->color_ = BLACK;
}

template<typename R, typename F>
void Genus<R,F>::left_rotate(GNODEPtr x)
{
    GNODEPtr y = x->right_;
    x->right_ = y->left_;
    if (y->left_ != nullptr) { y->left_->parent_ = x; }
    y->parent_ = x->parent_;

    if (x->parent_.use_count() == 0) { this->root_ = y; }
    else
    {
        GNODEPtr xp(x->parent_);
        if (x == xp->left_) { xp->left_ = y; }
        else { xp->right_ = y; }
    }

    y->left_ = x;
    x->parent_ = y;
}

template<typename R, typename F>
void Genus<R,F>::right_rotate(GNODEPtr x)
{
    GNODEPtr y = x->left_;
    x->left_ = y->right_;
    if (y->right_ != nullptr) { y->right_->parent_ = x; }
    y->parent_ = x->parent_;

    if (x->parent_.use_count() == 0) { this->root_ = y; }
    else
    {
        GNODEPtr xp(x->parent_);
        if (x == xp->left_) { xp->left_ = y; }
        else { xp->right_ = y; }
    }

    y->right_ = x;
    x->parent_ = y;
}

template<typename R, typename F>
void Genus<R,F>::print(GNODEPtr x, int64_t level) const
{
    if (x == nullptr) { return; }

    this->print(x->left_, level+1);
    std::cout << std::setw(4) << level << ":   " << x->q_ << std::endl;
    this->print(x->right_, level+1);
}

template<typename R, typename F>
void Genus<R,F>::print(void) const
{
    this->print(this->root_);
}

#endif // __GENUS_H_
