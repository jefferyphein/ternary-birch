#include "Genus.h"
#include "Prime.h"
#include "Math.h"
#include <memory>

typedef Genus<mpz_class, mpq_class> GenusZZ;
typedef Prime<mpz_class, mpq_class> PrimeZZ;

template<>
std::shared_ptr<PrimeZZ> GenusZZ::smallest_good_prime(void) const
{
    mpz_class d = this->disc_;

    if (d % 2 != 0)
    {
        return std::make_shared<PrimeZZ>(mpz_class(2));
    }

    mpz_class p = 3;

    while (Math::gcd(d, p) != 1) { p += 2; }

    return std::make_shared<PrimeZZ>(p);
}
