#ifndef __CHARACTER_H_
#define __CHARACTER_H_

#include "Isometry.h"
#include "QuadForm.h"

template<typename R, typename F>
class Character
{
public:
    Character() : ps_(std::vector<R>()), cond_(R(1)) {};
    Character(const std::vector<R>& ps);
    Character(const R& cond);

    mpz_class rho(const Isometry<R,F>& s, const QuadForm<R,F>& q) const;
    const R& conductor(void) const;
    const std::vector<R>& primes(void) const;
    bool operator<(const Character<R,F>& repr) const;
    
private:
    std::vector<R> ps_;
    R cond_;
};

template<typename R, typename F>
Character<R,F>::Character(const std::vector<R>& ps)
{
    this->ps_ = ps;

    this->cond_ = R(1);
    for (const R& p : ps)
    {
        this->cond_ *= p;
    }
}

template<typename R, typename F>
Character<R,F>::Character(const R& cond)
{
    this->cond_ = cond;

    // TODO: ADD A FACTORIZATION METHOD HERE.
    this->ps_ = std::vector<R>(1, cond);
}

template<typename R, typename F>
const std::vector<R>& Character<R,F>::primes(void) const
{
    return this->ps_;
}

template<typename R, typename F>
const R& Character<R,F>::conductor(void) const
{
    return this->cond_;
}

template<typename R, typename F>
bool Character<R,F>::operator<(const Character<R,F>& chi) const
{
    return this->cond_ < chi.conductor();
}

#endif // __CHARACTER_H_
