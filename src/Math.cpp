#include "Math.h"

template class Math<Z>;
template class Math<Z64>;

const std::vector<int> hilbert_lut_odd = { 1, 1, 1, 1,
                                           1, 1,-1,-1,
                                           1,-1, 1,-1,
                                           1,-1,-1, 1 };

const std::vector<int> hilbert_lut_p2 = { 1, 1, 1, 1, 1, 1, 1, 1,
                                          1,-1, 1,-1,-1, 1,-1, 1,
                                          1, 1, 1, 1,-1,-1,-1,-1,
                                          1,-1, 1,-1, 1,-1, 1,-1,
                                          1,-1,-1, 1, 1,-1,-1, 1,
                                          1, 1,-1,-1,-1,-1,-1, 1,
                                          1,-1,-1, 1,-1, 1, 1,-1,
                                          1, 1,-1,-1, 1, 1,-1,-1 };

template<>
int Z_Math::hilbert_symbol(Z a, Z b, const Z& p)
{
    int a_val = 0;
    while (a % p == 0)
    {   
        ++a_val;
        a /= p;
    }   

    int b_val = 0;
    while (b % p == 0)
    {   
        ++b_val;
        b /= p;
    }   

    if (p == 2)
    {   
        int aa = (mpz_class(a%8).get_si() >> 1) & 0x3;
        int bb = (mpz_class(b%8).get_si() >> 1) & 0x3;
        int index = ((a_val&0x1)<<5) | (aa << 3) | ((b_val&0x1)<<2) | bb; 
        return hilbert_lut_p2[index];
    }   

    int a_notsqr = mpz_legendre(a.get_mpz_t(), p.get_mpz_t()) == -1; 
    int b_notsqr = mpz_legendre(b.get_mpz_t(), p.get_mpz_t()) == -1; 

    int index = ((a_val&0x1)<<3) | (a_notsqr<<2) | ((b_val&0x1)<<1) | b_notsqr;
    if (((index & 0xa) == 0xa) && ((p%4) == 0x3))
    {   
        return -hilbert_lut_odd[index];
    }   
    else
    {   
        return hilbert_lut_odd[index];
    }   
}

template<>
int Z64_Math::hilbert_symbol(Z64 a, Z64 b, const Z64& p)
{
    return Z_Math::hilbert_symbol(a, b, p);
}
