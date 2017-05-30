#include "Ideal.h"

typedef Ideal<mpz_class, mpq_class> IdealZZ;

template<>
mpz_class IdealZZ::norm(void) const
{
    return abs(this->principal_generator());
}
