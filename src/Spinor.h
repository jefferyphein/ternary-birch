#ifndef __SPINOR_H_
#define __SPINOR_H_

#include "birch.h"

template<typename R>
class Spinor
{
    template<typename T>
    friend class Spinor;

public:
    Spinor(const std::vector<R>& primes)
    {
        this->primes_ = primes;
        this->twist = (1LL << this->primes_.size()) - 1;
    }

    Z64 norm(const QuadForm<R>& q, const Isometry<R>& s, const R& scalar) const
    {
        R tr = s.a11 + s.a22 + s.a33;
        if (tr != -scalar)
        {
            return this->compute_vals(tr + scalar);
        }

        R delta = 2 * q.a() * (s.a22 + s.a33) - (q.g() * s.a31 + q.h() * s.a21);
        if (delta != 0)
        {
            return this->compute_vals(delta) ^ this->twist;
        }

        R abh = 4 * q.a() * q.b() - q.h() * q.h();
        R abhm33 = abh * s.a33 + s.a32 * (q.g() * q.h() - 2 * q.a() * q.f()) + s.a31 * (q.f() * q.h() - 2 * q.b() * q.g());
        if (abhm33 != abh * scalar)
        {
            return this->compute_vals(abhm33 - abh*scalar) ^ this->compute_vals(2 * q.a()) ^ this->twist;
        }

        return this->compute_vals(abh * scalar);
    }

    const std::vector<R> primes(void) const
    {
        return this->primes_;
    }

private:
    std::vector<R> primes_;
    Z64 twist;

    Z64 compute_vals(R x) const
    {
        Z64 val = 0;
        Z64 mask = 1;
        for (const R& p : this->primes_)
        {
            while (x % p == 0)
            {
                x /= p;
                val ^= mask;
            }
            mask <<= 1;
        }
        return val;
    }
};

#endif // __SPINOR_H_
