#ifndef __REPRESENTATION_H
#define __REPRESENTATION_H

#include "Isometry.h"
#include "QuadForm.h"

template<typename R, typename F>
class Representation
{
public:
    Representation() = default;

    virtual R rho(const Isometry<R,F>& s, const QuadForm<R,F>& q) const;
    virtual R conductor(void) const;
};

template<typename R, typename F>
class Character : public Representation<R,F>
{
public:
    Character(const std::vector<R>& ps);
    Character(const R& d);

    R rho(const Isometry<R,F>& s, const QuadForm<R,F>& q) const;
private:
    std::vector<R> ps_;
    R cond_;
};

template<typename R, typename F>
Character<R, F>::Character(const std::vector<R>& ps)
{
    this->ps_ = ps;

    this->cond_ = R(1);
    for (const R& p : ps)
    {
        this->cond_ *= p;
    }
}

#endif // __REPRESENTATION_H
