#ifndef __NEIGHBOR_ITERATOR_H_
#define __NEIGHBOR_ITERATOR_H_

#include <memory>
#include "QuadForm.h"

template<typename R, typename F>
class NeighborIterator
{
public:
    NeighborIterator() = default;
    NeighborIterator(std::shared_ptr<QuadForm<R,F>> q, const R& p);

    inline int64_t num_neighbors(void) const;
    static int64_t num_neighbors(const R& p);
    inline const mpz_class& prime(void) const;

    inline std::shared_ptr<QuadForm<R,F>> quad_form(void) const;

    std::shared_ptr<QuadForm<R,F>> get_neighbor(int64_t pos) const;

    std::shared_ptr<QuadForm<R,F>> next_neighbor(void);

private:
    const std::shared_ptr<QuadForm<R,F>> build_neighbor(std::vector<R>& vec) const;

    std::shared_ptr<QuadForm<R,F>> q_;
    R p_;
    int64_t numNeighbors_;

    int64_t pos_;
    std::vector<R> isotropicVector_;
    std::shared_ptr<ChangeOfBasis<R,F>> s_;
    R aniso_;
};

template<typename R, typename F>
int64_t NeighborIterator<R,F>::num_neighbors(void) const
{
    return this->numNeighbors_;
}

template<typename R, typename F>
const mpz_class& NeighborIterator<R,F>::prime(void) const
{
    return this->p_;
}

template<typename R, typename F>
std::shared_ptr<QuadForm<R,F>> NeighborIterator<R,F>::quad_form(void) const
{
    return this->q_;
}

#endif // __NEIGHBOR_ITERATOR_H_
