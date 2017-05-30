#ifndef __IDEAL_H_
#define __IDEAL_H_

#include <vector>
#include <gmpxx.h>

template<typename R, typename F>
class Ideal
{
public:
    Ideal(const R& gen);
    Ideal(const std::vector<R>& gens);

    bool is_principal(void) const;

    R principal_generator(void) const;

    mpz_class norm(void) const;

private:
    std::vector<R> gens_;
};

template<typename R, typename F>
Ideal<R,F>::Ideal(const R& gen)
{
    this->gens_.push_back(gen);
}

template<typename R, typename F>
Ideal<R,F>::Ideal(const std::vector<R>& gens)
{
    for (const R& gen : gens)
    {
        this->gens_.push_back(gen);
    }
}

template<typename R, typename F>
bool Ideal<R,F>::is_principal(void) const
{
    return this->gens_.size() <= 1;
}

template<typename R, typename F>
R Ideal<R,F>::principal_generator(void) const
{
    size_t len = this->gens_.size();

    if (len > 1)
    {
        throw std::runtime_error("Ideal is not principal.");
    }

    if (len == 1)
    {
        return this->gens_[0];
    }

    return R(1);
}

#endif // __IDEAL_H_
