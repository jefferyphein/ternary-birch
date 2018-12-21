#include "birch.h"
#include "Fp.h"

template<>
template<>
W16 W16_Fp::mod(const Z& a) const
{
    Z r;
    return mpz_mod_ui(r.get_mpz_t(), a.get_mpz_t(), p);
}

template<>
template<>
W32 W32_Fp::mod(const Z& a) const
{
    Z r;
    return mpz_mod_ui(r.get_mpz_t(), a.get_mpz_t(), p);
}

template<>
template<>
W64 W64_Fp::mod(const Z& a) const
{
    Z r;
    return mpz_mod_ui(r.get_mpz_t(), a.get_mpz_t(), p);
}
