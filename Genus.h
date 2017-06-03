#ifndef __GENUS_H_
#define __GENUS_H_

#include <memory>
#include <iomanip>
#include <cassert>
#include <unordered_set>
#include <set>
#include "QuadForm.h"
#include "NeighborIterator.h"
#include "Character.h"

template<typename R, typename F>
class GenusRep
{
public:
    typedef std::shared_ptr<QuadForm<R,F>> QuadFormPtr;

    GenusRep(QuadFormPtr qPtr);

    int64_t pos(void) const;
    void pos(int64_t x);
    QuadFormPtr quad_form(void) const;

    bool operator==(const GenusRep<R,F>& rep) const;
    bool operator<(const GenusRep<R,F>& rep) const;
private:
    QuadFormPtr q_;
    int64_t pos_ = -1;
};

template<typename R, typename F>
class Genus
{
public:
    typedef std::shared_ptr<QuadForm<R,F>> QuadFormPtr;

    Genus(const QuadForm<R,F>& q);

    QuadFormPtr quad_form(void) const;

    R smallest_good_prime(void) const;

    void add_character(const Character<R,F>& chi);

    void print(void) const;

    size_t size(void) const;

private:
    /* An unordered set which stores distinct genus representatives. */
    std::unordered_set<GenusRep<R,F>> genusReps_;

    /* Computes the full genus. */
    void compute_genus(void);

    /* A vector of shared QuadForm pointers. An unordered list of the genus. */
    std::vector<QuadFormPtr> genusVec_;

    /* A pointer to the originating form. */
    QuadFormPtr q_;

    /* The discriminant. */
    R disc_;

    /* Flag indicating whether the genus has been computed in full. */
    bool computed_ = false;

    /* A set of characters to compute. */
    std::set<Character<R,F>> reprSet_;
};

template<typename R, typename F>
GenusRep<R,F>::GenusRep(QuadFormPtr q)
{
    this->q_ = q;
}

template<typename R, typename F>
int64_t GenusRep<R,F>::pos(void) const
{
    return this->pos_;
}

template<typename R, typename F>
void GenusRep<R,F>::pos(int64_t x)
{
    this->pos_ = x;
}

template<typename R, typename F>
std::shared_ptr<QuadForm<R,F>> GenusRep<R,F>::quad_form(void) const
{
    return this->q_;
}

template<typename R, typename F>
bool GenusRep<R,F>::operator==(const GenusRep<R,F>& rep) const
{
    return *(rep.quad_form()) == *this->q_;
}

template<typename R, typename F>
bool GenusRep<R,F>::operator<(const GenusRep<R,F>& rep) const
{
    return *this->q_ < *(rep.quad_form());
}

namespace std
{
    template<typename R, typename F>
    struct hash<GenusRep<R,F>>
    {
        int64_t operator()(const GenusRep<R,F>& rep) const
        {
            return rep.quad_form()->hash_value();
        }
    };
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
    return this->genusSet_.size();
}

template<typename R, typename F>
void Genus<R,F>::compute_genus(void)
{
    // Do nothing if the genus has already been computed.
    if (this->computed_) { return; }

    // The initial genus representative.
    GenusRep<R,F> gr(this->q_);

    // Initialize the genus.
    this->genusReps_.insert(gr);
    this->genusVec_.push_back(this->q_);

    // Determine the smallest good prime.
    R p = this->smallest_good_prime();

    // Loop over all genus representatives as the genus is built.
    int64_t index = 0;
    while (index < this->genusVec_.size())
    {
        // The current genus representative.
        QuadFormPtr cur = this->genusVec_[index++];

        // The neighbor iterator.
        NeighborIterator<R,F> it(cur, p);

        // The first neighbor.
        QuadFormPtr pn = it.next_neighbor();

        // Loop over all p-neighbors.
        while (pn != nullptr)
        {
            // Reduce the p-neighbor.
            QuadFormPtr qq = QuadForm<R,F>::reduce(*pn);

            // Create a genus representative object.
            GenusRep<R,F> rep(qq);

            // Check whether this form is already in the genus.
            if (this->genusReps_.count(rep) == 0)
            {
                this->genusReps_.insert(rep);
                this->genusVec_.push_back(qq);
            }

            // Get the next p-neighbor.
            pn = it.next_neighbor();
        }
    }
}

template<typename R, typename F>
void Genus<R,F>::add_character(const Character<R,F>& repr)
{
    if (this->reprSet_.count(repr) == 0)
    {
        this->reprSet_.insert(repr);
    }
}

template<typename R, typename F>
void Genus<R,F>::print(void) const
{
    //for (auto& q : this->genusVec_)
    //{
    //    auto vec = q->automorphisms();
    //    for (auto& aut : vec)
    //    {
    //        assert( aut->is_automorphism(*q) );
    //    }
    //}
    std::cout << this->genusVec_.size() << std::endl;
}

#endif // __GENUS_H_
