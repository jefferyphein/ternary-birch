#include "Genus.h"
#include "Math.h"
#include <memory>

typedef Genus<mpz_class, mpq_class> GenusZZ;

template<>
mpz_class GenusZZ::smallest_good_prime(void) const
{
    mpz_class d = this->disc_;

    if (d % 2 != 0)
    {
        return mpz_class(2);
    }

    mpz_class p = 3;

    while (Math::gcd(d, p) != 1) { p += 2; }

    return p;
}
