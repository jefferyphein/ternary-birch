#ifndef __PRIME_H_
#define __PRIME_H_

#include <vector>
#include "Ideal.h"

template <typename R, typename F>
class Prime : public Ideal<R,F>
{
public:
    Prime(const R& gen) : Ideal<R,F>(gen) {}
    Prime(const std::vector<R>& gens) : Ideal<R,F>(gens) {}
};

#endif // __PRIME_H_
