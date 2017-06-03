#include <memory>
#include <map>
#include "Genus.h"
#include "Math.h"

typedef Genus<mpz_class, mpq_class> GenusZZ;
typedef GenusRep<mpz_class, mpq_class> GenusRepZZ;
typedef Character<mpz_class, mpq_class> CharacterZZ;
typedef QuadForm<mpz_class, mpq_class> QuadFormZZ;

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

template<>
void GenusZZ::set_dimensions(GenusRepZZ& rep)
{
    auto auts = rep.quad_form()->automorphisms();

    // Initialize the dimensions for each character to be 1. This will be set
    // to zero if we ever see a rho value of -1.
    std::map<CharacterZZ, mpz_class> dimensionMap_;
    for (auto& chi : this->charSet_)
    {
        dimensionMap_[chi] = 1;
    }

    // For each automorphism, compute the conjugate so that it becomes an
    // automorphism on the original quadratic form.
    for (auto& aut : auts)
    {
        auto conj = aut->conjugate_by(*rep.quad_form()->isometry());

#ifdef DEBUG
        assert( conj->is_isometry(*this->q_, *this->q_) );
#endif

        // Loop over all primitive characters and compute their
        // representation values with respect to this automorphism.
        std::map<mpz_class, mpz_class> primeValues;
        for (auto& chi : this->primeCharSet_)
        {
#ifdef DEBUG
            assert( rep.quad_form()->isometry()->is_isometry(*this->q_, *rep.quad_form()) );
#endif

            primeValues[chi.conductor()] = chi.rho(*conj, *this->q_);
        }

        // Loop over all characters and compute their value based on the
        // values of the primitive characters we just computed.
        for (auto& chi : this->charSet_)
        {
            if (dimensionMap_[chi] == 0) { continue; }

            mpz_class value = 1;
            auto ps = chi.primes();
            for (auto& p : ps)
            {
                value *= primeValues[p];
            }

            // This character represents -1, hence the subspace fixed by
            // automorphisms is zero.
            if (value == -1)
            {
                dimensionMap_[chi] = 0;
            }
        }
    }

    // Set the dimensions for each character.
    for (auto& chi : this->charSet_)
    {
        rep.set_dimension(chi, dimensionMap_[chi]);
    }
}

