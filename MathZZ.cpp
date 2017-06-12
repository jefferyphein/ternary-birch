#include <cmath>
#include <gmpxx.h>
#include <iostream>
#include "Math.h"

typedef Math<mpz_class, mpq_class> MathZZ;

template<>
mpz_class MathZZ::gcd(const mpz_class& a, const mpz_class& b)
{
    mpz_class g;
    mpz_gcd(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return g;
}

template<>
mpz_class MathZZ::modinv(const mpz_class& a, const mpz_class& p)
{
    mpz_class g;
    mpz_class s;
    mpz_class t;
    mpz_gcdext(g.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t(), a.get_mpz_t(), p.get_mpz_t());
    return s;
}

template<>
mpz_class MathZZ::pow(const mpz_class& a, const mpz_class& e, const mpz_class& p)
{
    if (e == 0) { return 1; }

    if (e < 0) { return MathZZ::modinv(MathZZ::pow(a, -e, p), p); }

    // Handle the even power case.
    if ((e & 1) == 0)
    {
        mpz_class v = MathZZ::pow(a, e/2, p);
        return mpz_class((v * v) % p);
    }

    // Handle the odd power case.
    mpz_class v = MathZZ::pow(a, (e-1)/2, p);
    return mpz_class((a * ((v * v) % p)) % p);
}

template<>
void MathZZ::fix_vector(std::vector<mpz_class>& vec, const mpz_class& p)
{
    for (mpz_class& x : vec)
    {
        if (p == 2 && x != 0) { x = 1; }
        else if (x > (p-1)/2) { x -= p; }
        else if (-x > (p-1)/2) { x += p; }
    }
}

template<>
void MathZZ::normalize_vector(std::vector<mpz_class>& vec, const mpz_class& p)
{
    if (vec.size() != 3)
    {
        throw std::runtime_error("Vector must be three-dimensional.");
    }

    if (vec[0] != 0)
    {
        mpz_class inv = MathZZ::modinv(vec[0], p);
        vec[0] = 1;
        vec[1] = (vec[1] * inv) % p;
        vec[2] = (vec[2] * inv) % p;
    }
    else if (vec[1] != 0)
    {
        mpz_class inv = MathZZ::modinv(vec[1], p);
        vec[1] = 1;
        vec[2] = (vec[2] * inv) % p;
    }
    else
    {
        vec[2] = 1;
    }
}

template<>
int64_t MathZZ::kronecker(const mpz_class& a, const mpz_class& b)
{
    return static_cast<int64_t>(mpz_kronecker(a.get_mpz_t(), b.get_mpz_t()));
}

template<>
mpz_class MathZZ::square_root(mpz_class a, const mpz_class& p)
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
    if (S == 1) { return MathZZ::pow(a, (p+1)/4, p); }

    // Need to pick a quadratic nonresidue, we do so randomly.
    mpz_class z = -1;
    do { z = 1 + (rand() % (p-1)); }
    while (MathZZ::kronecker(z, p) == 1);

    // The Tonelli-Shanks algrithm follows.

    mpz_class c = MathZZ::pow(z, Q, p);
    mpz_class R = MathZZ::pow(a, (Q+1)/2, p);
    mpz_class t = MathZZ::pow(a, Q, p);
    mpz_class M = S;

    while (true)
    {
        if (t == 1) { return R; }
        mpz_class i = 0;
        mpz_class t1 = t;
        while (t1 != 1) { t1 = (t1 * t1) % p; i++; }

        mpz_class _pow = 1;
        mpz_class j; for (j = 0; j < M-i-1; j++) { _pow *= 2; }
        mpz_class b = MathZZ::pow(c, _pow, p) % p;
        R = (R * b) % p;
        t = (t * b * b) % p;
        c = (b * b) % p;
        M = i;
    }

    // Make the compiler happy.
    return 0;
}

template<>
int64_t MathZZ::valuation(mpq_class x, mpz_class p)
{
    mpz_class num = x.get_num();
    if (num == 0) { return 0; }

    mpz_class den = x.get_den();
    if (num == 0) { return 0; }

    int64_t count = 0;
    while (num % p == 0) { ++count; num /= p; }
    while (den % p == 0) { --count; den /= p; }
    return count;
}

template<>
std::vector<mpz_class> MathZZ::prime_divisors_naive(mpz_class n)
{
    std::vector<mpz_class> ps;

    while (n % 2 == 0)
    {
        ps.push_back(2);
        do
        {
            n /= 2;
        }
        while (n % 2 == 0);
    }

    mpz_class upper = sqrt(n);
    for (mpz_class p = 3; p <= upper; p += 2)
    {
        if (n % p == 0)
        {
            ps.push_back(p);
            do
            {
                n /= p;
            }
            while (n % p == 0);
        }
        upper = sqrt(n);
    }

    if (n != 1)
    {
        ps.push_back(n);
    }

    return ps;
}

template<>
std::vector<mpz_class> MathZZ::prime_divisors(const mpz_class& n)
{
    return MathZZ::prime_divisors_naive(n);
}

template<>
std::vector<std::vector<mpz_class>> MathZZ::squarefree_divisors(const mpz_class& n)
{
    std::vector<mpz_class> ps = MathZZ::prime_divisors(n);
    std::vector<std::vector<mpz_class>> divs(1, std::vector<mpz_class>());

    for (auto& p : ps)
    {
        size_t m = divs.size();
        for (size_t k = 0; k < m; k++)
        {
            auto vec = divs[k];
            vec.push_back(p);
            divs.push_back(vec);
        }
    }

    return divs;
}

template<>
std::vector<mpz_class> MathZZ::primes_up_to(const mpz_class& upTo, const mpz_class& coprimeTo)
{
    // Do not allow primes beyond the size of the largest int64_t value.
    if (upTo > mpz_class(std::numeric_limits<int64_t>::max()))
    {
        throw std::runtime_error("Supplied bound is too large.");
    }

    int64_t N = mpz_get_si(upTo.get_mpz_t());

    // Create the prime list.
    std::vector<mpz_class> primes;
    if (N >= 2 && coprimeTo % 2 != 0) { primes.push_back(2); }
    if (N >= 3 && coprimeTo % 3 != 0) { primes.push_back(3); }
    if (N < 4) { return primes; }

    // Initialize vector to start the sieve.
    std::vector<bool> isPrime(N/3, true);

    // We're going to iterate over values which are not divisible two 2 or 3.
    // This will step as follows: 5, 7, 11, 13, 17, 19, 23, 25, 29, 31, 35, etc.
    // In other words, starting at 5, the sequence steps by an increment of +2
    // then an increment of +4, repeatedly.
    for (int64_t k = 5, t = 2; k < N; k += t, t = 6-t)
    {
        if (isPrime[k/3])
        {
            if (MathZZ::gcd(k, coprimeTo) == abs(1))
            {
                primes.push_back(k);
            }

            for (int64_t j = k*k, v = t; j < N; j += v*k, v = 6-v)
            {
                isPrime[j/3] = false;
            }
        }
    }

    return primes;
}

template<>
bool MathZZ::is_squarefree(const mpz_class& x)
{
    std::vector<mpz_class> primes = std::move(MathZZ::prime_divisors_naive(x));
    mpz_class prod = 1;
    for (const mpz_class& p : primes)
    {
        prod *= p;
    }
    return prod == x;
}

template<>
bool MathZZ::is_positive(const mpz_class& x)
{
    return abs(x) == x;
}
