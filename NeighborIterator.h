#ifndef __NEIGHBOR_ITERATOR_H_
#define __NEIGHBOR_ITERATOR_H_

#include <memory>
#include "QuadForm.h"
#include "Prime.h"

template<typename R, typename F>
class NeighborIterator
{
public:
    NeighborIterator(std::shared_ptr<QuadForm<R,F>> q, std::shared_ptr<Prime<R,F>> pR);

    inline mpz_class num_neighbors(void) const;

    std::shared_ptr<QuadForm<R,F>> next_neighbor(void);

private:
    std::shared_ptr<QuadForm<R,F>> build_neighbor(std::vector<R>& vec);

    std::shared_ptr<QuadForm<R,F>> q_;
    std::shared_ptr<Prime<R,F>> p_;
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
