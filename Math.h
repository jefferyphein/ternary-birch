#ifndef __MATH_H_
#define __MATH_H_

#include <cmath>
#include <gmpxx.h>

class Math
{
public:
    static mpz_class gcd(const mpz_class& a, const mpz_class& b)
    {
        mpz_class g;
        mpz_gcd(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
        return g;
    }

    static mpz_class pow(const mpz_class& a, const mpz_class& e, const mpz_class& p)
    {
        if (e == 0) { return 1; }

        if (e < 0) { return Math::modinv(Math::pow(a, -e, p), p); }

        // Handle the even power case.
        if ((e & 1) == 0)
        {
            mpz_class v = Math::pow(a, e/2, p);
            return mpz_class((v * v) % p);
        }

        // Handle the odd power case.
        mpz_class v = Math::pow(a, (e-1)/2, p);
        return mpz_class((a * ((v * v) % p)) % p);
    }

    static mpz_class modinv(const mpz_class& a, const mpz_class& p)
    {
        mpz_class g;
        mpz_class s;
        mpz_class t;
        mpz_gcdext(g.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t(), a.get_mpz_t(), p.get_mpz_t());
        return s;
    }

    static void fix_vector(std::vector<mpz_class>& vec, const mpz_class& p)
    {
        for (mpz_class& x : vec)
        {
            if (x > (p-1)/2) { x -= p; }
            else if (-x > (p-1)/2) { x += p; }
        }
    }

    static int64_t kronecker(const mpz_class& a, const mpz_class& b)
    {
        return static_cast<int64_t>(mpz_kronecker(a.get_mpz_t(), b.get_mpz_t()));
    }

    static mpz_class square_root(mpz_class a, const mpz_class& p)
    {
        // Make a positive.
        a = a % p; if (a < 0) { a += p; }

        // Check the obvious case.
        if (a == 1) { return 1; }
        if (a == 0) { return 0; }

        // Set up the Tonelli-Shanks algorithm.
        mpz_class Q = p-1;
        mpz_class S = 0;
        while (Q % 2 == 0) { Q /= 2; S++; }

        // If S == 1, i.e. p = 3 (mod 4), then we have a solution!
        if (S == 1) { return Math::pow(a, (p+1)/4, p); }

        // Need to pick a quadratic nonresidue, we do so randomly.
        mpz_class z = -1;
        do { z = 1 + (rand() % (p-1)); }
        while (Math::kronecker(z, p) == 1);

        // The Tonelli-Shanks algrithm follows.

        mpz_class c = Math::pow(z, Q, p);
        mpz_class R = Math::pow(a, (Q+1)/2, p);
        mpz_class t = Math::pow(a, Q, p);
        mpz_class M = S;

        while (true)
        {
            if (t == 1) { return R; }
            mpz_class i = 0;
            mpz_class t1 = t;
            while (t1 != 1) { t1 = (t1 * t1) % p; i++; }

            mpz_class _pow = 1;
            mpz_class j; for (j = 0; j < M-i-1; j++) { _pow *= 2; }
            mpz_class b = Math::pow(c, _pow, p) % p;
            R = (R * b) % p;
            t = (t * b * b) % p;
            c = (b * b) % p;
            M = i;
        }

        // Make the compiler happy.
        return 0;
    }
};

#endif // __MATH_H_
