#ifndef __NEIGHBOR_ITERATOR_H_
#define __NEIGHBOR_ITERATOR_H_

#include <memory>
#include "QuadForm.h"

template<typename R, typename F>
class NeighborIterator
{
public:
    NeighborIterator(std::shared_ptr<QuadForm<R,F>> q, const R& p);

    inline mpz_class num_neighbors(void) const;

    std::shared_ptr<QuadForm<R,F>> next_neighbor(void);

private:
    std::shared_ptr<QuadForm<R,F>> build_neighbor(std::vector<R>& vec);

    std::shared_ptr<QuadForm<R,F>> q_;
    R p_;
    mpz_class numNeighbors_;

    mpz_class pos_;
    std::vector<R> isotropicVector_;
    std::shared_ptr<ChangeOfBasis<R,F>> s_;
    R aniso_;
};

template<typename R, typename F>
mpz_class NeighborIterator<R,F>::num_neighbors(void) const
{
    return this->numNeighbors_;
}

#endif // __NEIGHBOR_ITERATOR_H_
