#include "birch.h"
#include "Isometry.h"
#include "QuadForm.h"
#include "Math.h"

template<>
Z Z_QuadForm::discriminant(void) const
{
    mpz_t disc;
    mpz_t temp;
    mpz_init(disc);
    mpz_init(temp);

    mpz_mul(disc, this->b_.get_mpz_t(), this->c_.get_mpz_t());      // bc
    mpz_mul_2exp(disc, disc, 2);                                    // 4bc
    mpz_submul(disc, this->f_.get_mpz_t(), this->f_.get_mpz_t());   // 4bc-ff
    mpz_mul(disc, disc, this->a_.get_mpz_t());                      // 4abc-aff
    mpz_mul(temp, this->b_.get_mpz_t(), this->g_.get_mpz_t());      // bg
    mpz_submul(disc, temp, this->g_.get_mpz_t());                   // 4abc-aff-bgg
    mpz_mul(temp, this->f_.get_mpz_t(), this->g_.get_mpz_t());      // fg
    mpz_submul(temp, this->c_.get_mpz_t(), this->h_.get_mpz_t());   // fg-ch
    mpz_addmul(disc, temp, this->h_.get_mpz_t());                   // discriminant

    Z ret(disc);

    mpz_clear(temp);
    mpz_clear(disc);

    return ret;
}

template<>
bool Z_QuadForm::operator==(const Z_QuadForm& q) const
{
    return mpz_cmp(a_.get_mpz_t(), q.a_.get_mpz_t()) == 0 &&
           mpz_cmp(b_.get_mpz_t(), q.b_.get_mpz_t()) == 0 &&
           mpz_cmp(c_.get_mpz_t(), q.c_.get_mpz_t()) == 0 &&
           mpz_cmp(f_.get_mpz_t(), q.f_.get_mpz_t()) == 0 &&
           mpz_cmp(g_.get_mpz_t(), q.g_.get_mpz_t()) == 0 &&
           mpz_cmp(h_.get_mpz_t(), q.h_.get_mpz_t()) == 0;
}

template<>
Z Z_QuadForm::evaluate(const Z& x, const Z& y, const Z& z) const
{
    mpz_t value;
    mpz_t temp;
    mpz_init(value);
    mpz_init(temp);

    mpz_mul(value, this->a_.get_mpz_t(), x.get_mpz_t());    // ax
    mpz_addmul(value, this->g_.get_mpz_t(), z.get_mpz_t()); // ax+gz
    mpz_addmul(value, this->h_.get_mpz_t(), y.get_mpz_t()); // ax+gz+hy
    mpz_mul(value, value, x.get_mpz_t());                   // axx+gxz+hxy
    mpz_mul(temp, this->b_.get_mpz_t(), y.get_mpz_t());     // by
    mpz_addmul(temp, this->f_.get_mpz_t(), z.get_mpz_t());  // by+fz
    mpz_addmul(value, temp, y.get_mpz_t());                 // axx+byy+fyz+gxz+hxy
    mpz_mul(temp, this->c_.get_mpz_t(), z.get_mpz_t());     // cz
    mpz_addmul(value, temp, z.get_mpz_t());                 // evaluate

    Z ret(value);

    mpz_clear(temp);
    mpz_clear(value);

    return ret;
}

template<>
Z_QuadForm Z_QuadForm::reduce(const Z_QuadForm& q, Z_Isometry& s)
{
    mpz_t a, b, c, f, g, h;
    mpz_init_set(a, q.a_.get_mpz_t());
    mpz_init_set(b, q.b_.get_mpz_t());
    mpz_init_set(c, q.c_.get_mpz_t());
    mpz_init_set(f, q.f_.get_mpz_t());
    mpz_init_set(g, q.g_.get_mpz_t());
    mpz_init_set(h, q.h_.get_mpz_t());

    // TODO: Move these into some kind of "scratch" space so that we can
    // avoid reallocating these every time we call this function.
    mpz_t t, num, den, temp, temp2, temp3;
    mpz_inits(t, num, den, temp, temp2, temp3, NULL);

    int flag = 1;
    while (flag)
    {
        mpz_add(t, a, b);
        mpz_add(t, t, f);
        mpz_add(t, t, g);
        mpz_add(t, t, h);
        if (mpz_cmp_ui(t, 0) < 0)
        {
            s.A101011001();
            mpz_add(c, c, t);
            mpz_add(f, f, h);
            mpz_add(f, f, b);
            mpz_add(f, f, b);
            mpz_add(g, g, h);
            mpz_add(g, g, a);
            mpz_add(g, g, a);
        }

        mpz_add(den, a, a);
        if (mpz_cmp(a, h) >= 0)
        {
            mpz_sub(num, a, h);
            mpz_fdiv_q(t, num, den);
        }
        else
        {
            mpz_add(num, a, h);
            mpz_sub_ui(num, num, 1);
            mpz_fdiv_q(t, num, den);
            mpz_neg(t, t);
        }
        if (mpz_cmp_ui(t, 0) != 0)
        {
            s.A1t0010001(Z(t));
            mpz_mul(temp, a, t);
            mpz_add(h, h, temp);
            mpz_addmul(b, h, t);
            mpz_addmul(f, g, t);
            mpz_add(h, h, temp);
        }

        mpz_add(den, b, b);
        if (mpz_cmp(b, f) >= 0)
        {
            mpz_sub(num, b, f);
            mpz_fdiv_q(t, num, den);
        }
        else
        {
            mpz_add(num, b, f);
            mpz_sub_ui(num, num, 1);
            mpz_fdiv_q(t, num, den);
            mpz_neg(t, t);
        }
        if (mpz_cmp_ui(t, 0) != 0)
        {
            s.A10001t001(Z(t));
            mpz_mul(temp, b, t);
            mpz_add(f, f, temp);
            mpz_addmul(c, f, t);
            mpz_addmul(g, h, t);
            mpz_add(f, f, temp);
        }

        mpz_add(den, a, a);
        if (mpz_cmp(a, g) >= 0)
        {
            mpz_sub(num, a, g);
            mpz_fdiv_q(t, num, den);
        }
        else
        {
            mpz_add(num, a, g);
            mpz_sub_ui(num, num, 1);
            mpz_fdiv_q(t, num, den);
            mpz_neg(t, t);
        }
        if (mpz_cmp_ui(t, 0) != 0)
        {
            s.A10t010001(Z(t));
            mpz_mul(temp, a, t);
            mpz_add(g, g, temp);
            mpz_addmul(c, g, t);
            mpz_addmul(f, h, t);
            mpz_add(g, g, temp);
        }

        if (mpz_cmp(a, b) > 0 || (mpz_cmp(a, b) == 0 && mpz_cmpabs(f, g) > 0))
        {
            s.A0n0n0000n();
            mpz_swap(a, b);
            mpz_swap(f, g);
        }

        if (mpz_cmp(b, c) > 0 || (mpz_cmp(b, c) == 0 && mpz_cmpabs(g, h) > 0))
        {
            s.An0000n0n0();
            mpz_swap(b, c);
            mpz_swap(g, h);
        }

        if (mpz_cmp(a, b) > 0 || (mpz_cmp(a, b) == 0 && mpz_cmpabs(f, g) > 0))
        {
            s.A0n0n0000n();
            mpz_swap(a, b);
            mpz_swap(f, g);
        }

        int fgh = (mpz_cmp_ui(f, 0) != 0) &&
                  (mpz_cmp_ui(g, 0) != 0) &&
                  (mpz_cmp_ui(h, 0) != 0);
        if (fgh)
        {
            if (mpz_cmp_ui(f, 0) < 0) fgh = !fgh;
            if (mpz_cmp_ui(g, 0) < 0) fgh = !fgh;
            if (mpz_cmp_ui(h, 0) < 0) fgh = !fgh;
        }

        if (fgh)
        {
            if (mpz_cmp_ui(f, 0) < 0)
            {
                s.An00010001();
                mpz_neg(f, f);
            }

            if (mpz_cmp_ui(g, 0) < 0)
            {
                s.A1000n0001();
                mpz_neg(g, g);
            }

            if (mpz_cmp_ui(h, 0) < 0)
            {
                s.A10001000n();
                mpz_neg(h, h);
            }
        }
        else
        {
            int s1 = mpz_cmp_ui(f, 0) > 0;
            int s2 = mpz_cmp_ui(g, 0) > 0;
            int s3 = mpz_cmp_ui(h, 0) > 0;

            if ((s1+s2+s3) % 2 == 1)
            {
                if (mpz_cmp_ui(f, 0) == 0) s1 = 1;
                else
                {
                    if (mpz_cmp_ui(g, 0) == 0) s2 = 1;
                    else if (mpz_cmp_ui(h, 0) == 0) s3 = 1;
                }
            }

            if (s1 == 1)
            {
                s.An00010001();
                mpz_neg(f, f);
            }

            if (s2 == 1)
            {
                s.A1000n0001();
                mpz_neg(g, g);
            }

            if (s3 == 1)
            {
                s.A10001000n();
                mpz_neg(h, h);
            }
        }

        mpz_add(temp, a, b);
        mpz_add(temp, temp, f);
        mpz_add(temp, temp, g);
        mpz_add(temp, temp, h);

        flag = !(mpz_cmpabs(f, b) <= 0 && mpz_cmpabs(g, a) <= 0 &&
                 mpz_cmpabs(h, a) <= 0 && mpz_cmp_ui(temp, 0) >= 0);
    }

    mpz_add(temp3, a, a);
    mpz_add(temp3, temp3, h);
    mpz_add(temp2, temp3, g);
    mpz_add(temp2, temp2, g);

    if (mpz_cmp_ui(temp, 0) == 0 && mpz_cmp_ui(temp2, 0) > 0)
    {
        s.An010n1001();
        mpz_add(f, f, h);
        mpz_add(f, f, b);
        mpz_add(f, f, b);
        mpz_neg(f, f);
        mpz_add(g, g, temp3);
        mpz_neg(g, g);
    }

    if (mpz_cmp_ui(g, 0) != 0)
    {
        mpz_neg(temp, h);
        if (mpz_cmp(a, temp) == 0)
        {
            s.Ann00n0001();
            mpz_add(f, f, g);
            mpz_neg(f, f);
            mpz_neg(g, g);
            mpz_neg(h, h);
        }
    }

    if (mpz_cmp_ui(h, 0) != 0)
    {
        mpz_neg(temp, g);
        if (mpz_cmp(a, temp) == 0)
        {
            s.An0n01000n();
            mpz_add(f, f, h);
            mpz_neg(f, f);
            mpz_neg(h, h);
            mpz_add(g, g, a);
            mpz_add(g, g, a);
        }
    }

    if (mpz_cmp_ui(h, 0) != 0)
    {
        mpz_neg(temp, f);
        if (mpz_cmp(b, temp) == 0)
        {
            s.A1000nn00n();
            mpz_add(g, g, h);
            mpz_neg(g, g);
            mpz_neg(h, h);
            mpz_add(f, f, b);
            mpz_add(f, f, b);
        }
    }

    if (mpz_cmp(a, h) == 0)
    {
        mpz_add(temp, f, f);
        if (mpz_cmp(g, temp) > 0)
        {
            s.Ann001000n();
            mpz_sub(f, g, f);
        }
    }

    if (mpz_cmp(a, g) == 0)
    {
        mpz_add(temp, f, f);
        if (mpz_cmp(h, temp) > 0)
        {
            s.An0n0n0001();
            mpz_sub(f, h, f);
        }
    }

    if (mpz_cmp(b, f) == 0)
    {
        mpz_add(temp, g, g);
        if (mpz_cmp(h, temp) > 0)
        {
            s.An000nn001();
            mpz_sub(g, h, g);
        }
    }

    if (mpz_cmp(a, b) == 0 && mpz_cmpabs(f, g) > 0)
    {
        s.A0n0n0000n();
        mpz_swap(a, b);
        mpz_swap(f, g);
    }

    if (mpz_cmp(b, c) == 0 && mpz_cmpabs(g, h) > 0)
    {
        s.An0000n0n0();
        mpz_swap(g, h);
    }

    if (mpz_cmp(a, b) == 0 && mpz_cmpabs(f, g) > 0)
    {
        s.A0n0n0000n();
        mpz_swap(f, g);
    }

    Z_QuadForm ret = Z_QuadForm(Z(a), Z(b), Z(c), Z(f), Z(g), Z(h));

    mpz_clears(a, b, c, f, g, h, t, num, den, temp, temp2, temp3, NULL);

    return ret;
}

static W64 sign_vector(const Z& x, const Z& det, const std::vector<Z_PrimeSymbol>& primes)
{
    W64 vec = (x == -1);
    for (const Z_PrimeSymbol& symb : primes)
    {
        int value = Z_Math::hilbert_symbol(x, -det, symb.p);
        vec = (vec << 1) | (value == -1);
    }
    return vec;
}

// A naive GF(2) solver that looks for solutions to a specific linear equation
// from a specified starting point.
static W64 GF2_solve_naive(const std::vector<W64>& vecs, W64 start, W64 target)
{
    W64 upper = 1LL << vecs.size();
    size_t num_vecs = vecs.size();
    for (W64 i=start+1; i<upper; i++)
    {   
        W64 x = 0;
        W64 mask = upper >> 1;
        for (size_t j=0; j<num_vecs; j++)
        {
            if (i & mask) x ^= vecs[j];
            mask >>= 1;
        }

        if (x == target) return i;
    }   

    return 0;
}

// A simple brute force search for p^2-isotropic vectors. This can probably be
// rewritten to avoid an exhaustive search, but since we expect the primes to
// be small, this should work for now.
static Z_Vector3 Z_isotropic_mod_pp(const Z_QuadForm& q, const Z& p)
{
    Z pp = p*p;
    Z_Vector3 vec = {0,0,1};

    // Try (0 0 1) first.
    if (q.evaluate(vec) % pp == 0)
        return vec;

    // Try (0 1 x) next.
    vec.y = 1;
    for (vec.z = 0; vec.z < pp; vec.z++)
        if (q.evaluate(vec) % pp == 0)
            return vec;

    // Lastly, try (1 x y).
    vec.x = 1;
    for (vec.y = 0; vec.y < pp; vec.y++)
        for (vec.z = 0; vec.z < pp; vec.z++)
            if (q.evaluate(vec) % pp == 0)
                return vec;

    // Otherwise, return the zero vector to indicate no solution found.
    vec = {0};
    return vec;
}

template<>
Z_QuadForm Z_QuadForm::get_quad_form(const std::vector<Z_PrimeSymbol>& input)
{
    Z det = 1;
    Z disc = 1;
    bool two_specified = false;

    int num_ramified = 0;
    for (const Z_PrimeSymbol& symb : input)
    {
        if (symb.power % 2 == 0)
        {
            throw std::invalid_argument("Prime powers must be odd.");
        }

        if (symb.p == 2)
        {
            two_specified = true;
        }

        if (symb.ramified)
        {
            ++num_ramified;
        }

        det *= symb.p;
        for (int k=0; k<symb.power; k++)
        {
            disc *= symb.p;
        }
    }

    if (num_ramified % 2 == 0)
    {
        throw std::invalid_argument("Must specify an odd number of ramified primes.");
    }

    // Make a copy of the prime symbols to include 2 if it wasn't already
    // included. We also build the target vector encoding the desired
    // ramification properties.
    std::vector<Z_PrimeSymbol> primes;
    W64 target = 1; // Ramified at the infinite place.

    if (!two_specified)
    {
        Z_PrimeSymbol s;
        s.p = 2;
        s.power = 0;
        s.ramified = 0;
        primes.push_back(s);
        target <<= 1;
    }
    for (const Z_PrimeSymbol& symb : input)
    {
        primes.push_back(symb);
        target <<= 1;
        if (symb.ramified)
        {
            target |= 1;
        }
    }

    // The list of primes used in our factor baes.
    std::vector<Z> fullbase;

    // Create an std::vector consisting of GF(2)-vectors encoding Hilbert
    // symbol values.
    std::vector<W64> signs;

    // Add the relation for the infinite prime.
    fullbase.push_back(-1);
    signs.push_back(sign_vector(-1, det, primes));

    for (const Z_PrimeSymbol& symb : primes)
    {
        signs.push_back(sign_vector(symb.p, det, primes));
        fullbase.push_back(symb.p);
    }

    Z p = 2;
    W64 solution;
    bool done = false;
    bool added_to_end = false;

    while (!done)
    {
        solution = 0;
        do
        {
            solution = GF2_solve_naive(signs, solution, target);

            if (solution)
            {
                W64 mask = 1LL << fullbase.size();
                Z b = -1;
                Z det2 = det;

                // Construct the quadratic space associated to the solution
                // we've found.
                for (const Z& q : fullbase)
                {
                    mask >>= 1;
                    if (solution & mask)
                    {
                        b *= q;
                        if (q > 0)
                        {
                            if (det2 % q == 0) det2 /= q;
                            else
                            {
                                det2 *= q;
                            }
                        }
                    }
                }

                // Verify that the Hilbert symbols are correct at inf, 2, and
                // the primes dividing the discriminant.
                mask = 1LL << primes.size();
                for (const Z_PrimeSymbol& symb : primes)
                {
                    mask >>= 1;
                    int sign = (target & mask) ? -1 : 1;
                    if (Z_Math::hilbert_symbol(-b, -disc, symb.p) != sign)
                    {
                        throw std::runtime_error("Hilbert symbols do not check out. How did this happen?");
                    }
                }

                // We assume that the quadratic space is good outside of the
                // initial specifications. If we detect that one of our
                // auxilliary primes is ramified, we immediately flag this
                // space as not good and try another solution.
                int good = true;
                for (size_t n=primes.size()+1; n<fullbase.size(); n++)
                {
                    int sign = Z_Math::hilbert_symbol(-b, -disc, fullbase[n]);
                    if (sign == -1) good = false;
                }

                // If good, we specify as such and break out of the while loop
                // searching for solutions.
                if (good)
                {
                    done = true;
                    break;
                }
            }
        }
        while (solution != 0);

        // If we flagged done within our while-loop above, i.e. we found a
        // quadratic space with the correct ramification behavior, then we
        // break out of the factor base building loop.
        if (done) break;

        // Get the next prime not dividing the discriminant...
        mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
        while (disc % p == 0)
        {
            mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
        }

        // If we've added a prime to the end of our factor base list, we
        // remove it if it didn't lead to a solution. This helps avoid the
        // factor base growing too large, although it can lead to the size of
        // the primes in the factor base becoming larger than we'd like. (This
        // does not appear to be a problem in practice.)
        if (added_to_end)
        {
            signs.pop_back();
            fullbase.pop_back();
        }

        // ...and push it's sign vector onto the list.
        signs.push_back(sign_vector(p, det, primes));
        fullbase.push_back(p);
        added_to_end = true;
    }

    // Construct the coefficients of the diagonalized quadratic space with
    // the desired signs.
    Z b = -1;
    W64 mask = 1LL << fullbase.size();

    // Set up the base primes. These are primes which do not divide the
    // discriminant, but were included (as squares) in the discriminant of
    // the resulting quadratic space. We also include 2 in this list, as
    // additional powers of 2 can sometimes occur, even if 2 divides the
    // discriminant.
    std::vector<Z> base;
    base.push_back(2);

    for (const Z& p : fullbase)
    {
        mask >>= 1;
        if (solution & mask)
        {
            b *= p;
            if (p > 0)
            {
                if (det % p == 0) det /= p;
                else
                {
                    // Add p to the base primes.
                    det *= p;
                    base.push_back(p);
                }
            }
        }
    }

    // Verify that the Hilbert symbols are correct.
    mask = 1LL << primes.size();
    for (const Z_PrimeSymbol& symb : primes)
    {
        mask >>= 1;
        int sign = (target & mask) ? -1 : 1;
        if (Z_Math::hilbert_symbol(-b, -disc, symb.p) != sign)
        {
            throw std::runtime_error("Hilbert symbols do not check out. How did this happen?");
        }
    }

    // Start by building a diagonal form.
    Z a = 1;
    Z c = det;
    Z f = 0;
    Z g = 0;
    Z h = 0;
    Z_QuadForm q(a, b, c, f, g, h);
    Z N = q.discriminant();

    // Remove the prime squares from the discriminant for those primes not
    // dividing the intended discriminant.
    for (const Z& p : base)
    {
        Z pp = p*p;

        while (N % pp == 0)
        {
            // Find an isotropic vector mod p^2, make it a basis vector,
            // and then divide p^2 out of the discriminant.
            Z_Vector3 vec = Z_isotropic_mod_pp(q, p);

            #ifdef DEBUG
            assert( q.evaluate(vec) % pp == 0 );
            #endif

            if (vec.x == 0 && vec.y == 0 && vec.z == 0) break;
            else if (vec.x == 0 && vec.y == 0)
            {
                #ifdef DEBUG
                assert( vec.z == 1 );
                assert( g % p == 0 );
                assert( f % p == 0 );
                assert( c % pp == 0 );
                #endif

                c /= pp;
                f /= p;
                g /= p;
            }
            else if (vec.x == 0)
            {
                b += (c*vec.z*vec.z + f*vec.z);
                f += (2*c*vec.z);
                h += (g*vec.z);

                #ifdef DEBUG
                assert( vec.y == 1 );
                assert( b % pp == 0 );
                assert( f % p == 0 );
                assert( h % p == 0 );
                #endif

                b /= pp;
                f /= p;
                h /= p;
            }
            else
            {
                a += (b*vec.y*vec.y + c*vec.z*vec.z + f*vec.y*vec.z + g*vec.z + h*vec.y);
                g += (2*c*vec.z + f*vec.y);
                h += (2*b*vec.y + f*vec.z);

                #ifdef DEBUG
                assert( vec.x == 1 );
                assert( a % pp == 0 );
                assert( g % p == 0 );
                assert( h % p == 0 );
                #endif

                a /= pp;
                g /= p;
                h /= p;
            }

            q = Z_QuadForm(a, b, c, f, g, h);
            N = q.discriminant();
        }
    }

    // Reduce the resulting quadratic form, scale up to the correct
    // discriminant, then reduce again. The resulting form will have the
    // correct local behavior as well as the correct discriminant.
    Z_Isometry s;
    q = Z_QuadForm::reduce(q, s);
    for (const Z_PrimeSymbol& symb : primes)
    {
        for (int j=3; j<=symb.power; j+=2)
        {
            q.f_ *= symb.p;
            q.g_ *= symb.p;
            q.c_ *= symb.p * symb.p;
        }
    }
    q = Z_QuadForm::reduce(q, s);

    // Do one final verification that the symbols are correct.
    Z x = 4 * q.a_ * q.b_ - q.h_ * q.h_;
    mask = 1LL << primes.size();
    for (const Z_PrimeSymbol& symb : primes)
    {
        mask >>= 1;
        int sign = (target & mask) ? -1 : 1;
        if (Z_Math::hilbert_symbol(-x, -disc, symb.p) != sign)
        {
            throw std::runtime_error("Hilbert symbols do not check out. How did this happen?");
        }
    }

    return q;
}

template<>
W64 Z_QuadForm::hash_value(void) const
{
    W64 fnv = FNV_OFFSET;
    fnv = (fnv ^ mpz_get_si(this->a_.get_mpz_t())) * FNV_PRIME;
    fnv = (fnv ^ mpz_get_si(this->b_.get_mpz_t())) * FNV_PRIME;
    fnv = (fnv ^ mpz_get_si(this->c_.get_mpz_t())) * FNV_PRIME;
    fnv = (fnv ^ mpz_get_si(this->f_.get_mpz_t())) * FNV_PRIME;
    fnv = (fnv ^ mpz_get_si(this->g_.get_mpz_t())) * FNV_PRIME;
    fnv = (fnv ^ mpz_get_si(this->h_.get_mpz_t())) * FNV_PRIME;
    return fnv;
}

template<>
W64 Z64_QuadForm::hash_value(void) const
{
    W64 fnv = FNV_OFFSET;
    fnv = (fnv ^ this->a_) * FNV_PRIME;
    fnv = (fnv ^ this->b_) * FNV_PRIME;
    fnv = (fnv ^ this->c_) * FNV_PRIME;
    fnv = (fnv ^ this->f_) * FNV_PRIME;
    fnv = (fnv ^ this->g_) * FNV_PRIME;
    fnv = (fnv ^ this->h_) * FNV_PRIME;
    return fnv;
}
