#include <gmpxx.h>
#include <cassert>
#include "QuadForm.h"
#include "Math.h"
#include "AutomorphismZZ.h"

// A constant used to determine whether int64_t or mpz_class data types should
// be used for reducing quadratic forms. If all coefficients of the supplied
// quadratic form are smaller than this value, we will use int64_t, however if
// any one coefficient is larger, we will use mpz_class.
//
// TODO: Perform some analysis and determine a tight bound on this value, so
// that we can guarantee not to exceed std::numeric_limits<int64_t>::{min,max}
// during a reduction call. The current value is merely an experimental value.
static constexpr int64_t UPPERBOUND = 1L << 40;

typedef QuadForm<mpz_class, mpq_class> QuadFormZZ;
typedef Isometry<mpz_class, mpq_class> IsometryQQ;
typedef std::shared_ptr<IsometryQQ> IsometryQQPtr;
typedef Math<mpz_class, mpq_class> MathZZ;

// This function is templatized to allow for native int64_t arithmetic. This
// is significantly faster than mpz_class arithmetic, but may lead to overflow
// issues when the coefficients are large enough. It is therefore important to
// only set T=int64_t when the coefficients are sufficiently small; in all
// other cases, set T=mpz_class.
//
// NOTE: The coefficients a, b, c, f, g, h are never multiplied together
// during this routine, they are only ever added or subtracted. The only
// exception being the computation of the variable `t` during the first while
// loop.
template<typename R, typename F, typename T>
std::shared_ptr<QuadFormZZ> reduceT(const QuadFormZZ& q,
                                    T a, T b, T c, T f, T g, T h,
                                    bool saveIsometry)
{
    // Isometry coefficients.
    T a11 = 1;
    T a12 = 0;
    T a13 = 0;
    T a21 = 0;
    T a22 = 1;
    T a23 = 0;
    T a31 = 0;
    T a32 = 0;
    T a33 = 1;

    // Flag controlling the initial reduction loop.
    bool flag = true;

    // Temporary variable.
    T temp;

    while (flag)
    {
        T t = a + b + f + g + h;
        if (t < 0)
        {
            if (saveIsometry)
            {
                // transform on the right by...
                // |  1  0  1 |
                // |  0  1  1 |
                // |  0  0  1 |
                a13 += (a11 + a12);
                a23 += (a21 + a22);
                a33 += (a31 + a32);
            }

            // apply the isometry.
            c += t;
            f += (h + b + b);
            g += (h + a + a);
        }

        t = (a-h) / (a+a);
        if (a < h && (a-h) % (a+a) != 0) { t--; }
        if (t != 0)
        {
            if (saveIsometry)
            {
                // transform on the right by...
                // |  1  t  0 |
                // |  0  1  0 |
                // |  0  0  1 |
                a12 += t * a11;
                a22 += t * a21;
                a32 += t * a31;
            }

            // apply the isometry.
            f += g * t;
            b += t * (a * t + h);
            h += t * (a + a);
        }

        ///////////////////////

        t = (b-f) / (b+b);
        if (b < f && (b-f) % (b+b) != 0) { t--; }
        if (t != 0)
        {
            if (saveIsometry)
            {
                // transform on the right by...
                // |  1  0  0 |
                // |  0  1  t |
                // |  0  0  1 |
                a13 += t * a12;
                a23 += t * a22;
                a33 += t * a32;
            }

            // apply the isometry.
            g += h * t;
            c += t * (b * t + f);
            f += t * (b + b);
        }

        ///////////////////////

        t = (a-g) / (a+a);
        if (a < g && (a-g) % (a+a) != 0) { t--; }
        if (t != 0)
        {
            if (saveIsometry)
            {
                // transform on the right by...
                // |  1  0  t |
                // |  0  1  0 |
                // |  0  0  1 |
                a13 += t * a11;
                a23 += t * a21;
                a33 += t * a31;
            }

            // apply the isometry.
            f += h * t;
            c += t * (a * t + g);
            g += t * (a + a);
        }

        if (a > b || (a == b && abs(f) > abs(g)))
        {
            if (saveIsometry)
            {
                // transform on the right by...
                // |  0 -1  0 |
                // | -1  0  0 |
                // |  0  0 -1 |
                temp = a11; a11 = -a12; a12 = -temp; a13 = -a13;
                temp = a21; a21 = -a22; a22 = -temp; a23 = -a23;
                temp = a31; a31 = -a32; a32 = -temp; a33 = -a33;
            }

            // apply the isometry.
            std::swap(a, b);
            std::swap(f, g);
        }

        if (b > c || (b == c && abs(g) > abs(h)))
        {
            if (saveIsometry)
            {
                // transform on the right by...
                // | -1  0  0 |
                // |  0  0 -1 |
                // |  0 -1  0 |
                temp = a12; a12 = -a13; a13 = -temp; a11 = -a11;
                temp = a22; a22 = -a23; a23 = -temp; a21 = -a21;
                temp = a32; a32 = -a33; a33 = -temp; a31 = -a31;
            }

            // apply the isometry.
            std::swap(b, c);
            std::swap(g, h);
        }

        if (a > b || (a == b && abs(f) > abs(g)))
        {
            if (saveIsometry)
            {
                // transform on the right by...
                // |  0 -1  0 |
                // | -1  0  0 |
                // |  0  0 -1 |
                temp = a11; a11 = -a12; a12 = -temp; a13 = -a13;
                temp = a21; a21 = -a22; a22 = -temp; a23 = -a23;
                temp = a31; a31 = -a32; a32 = -temp; a33 = -a33;
            }

            // apply the isometry.
            std::swap(a, b);
            std::swap(f, g);
        }

        // Determine the sign of f*g*h.
        bool fgh = (f != 0 && g != 0 && h != 0); 
        if (fgh)
        {
            if (f < 0) fgh = !fgh;
            if (g < 0) fgh = !fgh;
            if (h < 0) fgh = !fgh;
        }

        if (fgh)
        {
            // note that for us to enter this block, we must have either a pair of
            //  f, g, h negative, or none of them. in what follows, we do not
            //  actually apply the transformation listed to the coefficients, but
            //  rather update values according to what the net result will be once
            //  we exit this block.
            if (f < 0)
            {
                if (saveIsometry)
                {
                    // transform on the right by...
                    // | -1  0  0 |
                    // |  0  1  0 |
                    // |  0  0  1 |
                    a11 = -a11;
                    a21 = -a21;
                    a31 = -a31;
                }

                // negate f.
                f = -f;
            }
            if (g < 0)
            {
                if (saveIsometry)
                {
                    // transform on the right by...
                    // |  1  0  0 |
                    // |  0 -1  0 |
                    // |  0  0  1 |
                    a12 = -a12;
                    a22 = -a22;
                    a32 = -a32;
                }

                // negate g.
                g = -g;
            }
            if (h < 0)
            {
                if (saveIsometry)
                {
                    // transform on the right by...
                    // |  1  0  0 |
                    // |  0  1  0 |
                    // |  0  0 -1 |
                    a13 = -a13;
                    a23 = -a23;
                    a33 = -a33;
                }

                // negate h.
                h = -h;
            }
        }
        else
        {
            // the same reasoning from the previous block applies here. we must
            //  have either a nonzero coefficient or an odd number of them
            //  negative. we build the isometry, but only update the coefficients
            //  as needed to make them correct once we exit this block.

            short s1 = f > 0;
            short s2 = g > 0;
            short s3 = h > 0;

            if ((s1 + s2 + s3) % 2 == 1)
            {
                if (f == 0) { s1 = 1; }
                else
                {
                    if (g == 0) { s2 = 1; }
                    else
                    {
                        if (h == 0) { s3 = 1; }
                    }
                }
            }

            if (s1 == 1)
            {
                if (saveIsometry)
                {
                    // transform on the right by...
                    // | -1  0  0 |
                    // |  0  1  0 |
                    // |  0  0  1 |
                    a11 = -a11;
                    a21 = -a21;
                    a31 = -a31;
                }

                // negate f.
                f = -f;
            }

            if (s2 == 1)
            {
                if (saveIsometry)
                {
                    // transform on the right by...
                    // |  1  0  0 |
                    // |  0 -1  0 |
                    // |  0  0  1 |
                    a12 = -a12;
                    a22 = -a22;
                    a32 = -a32;
                }

                // negate g.
                g = -g;
            }

            if (s3 == 1)
            {
                if (saveIsometry)
                {
                    // transform on the right by...
                    // |  1  0  0 |
                    // |  0  1  0 |
                    // |  0  0 -1 |
                    a13 = -a13;
                    a23 = -a23;
                    a33 = -a33;
                }

                // negate h.
                h = -h;
            }
        }

        flag = !(abs(f) <= b && abs(g) <= a && abs(h) <= a && a+b+f+g+h >= 0);
    }

    if (a+b+f+g+h == 0 && a+a+g+g+h > 0)
    {
        if (saveIsometry)
        {
            // transform on the right by...
            // | -1  0  1 |
            // |  0 -1  1 |
            // |  0  0  1 |
            a13 += (a11 + a12); a11 = -a11; a12 = -a12;
            a23 += (a21 + a22); a21 = -a21; a22 = -a22;
            a33 += (a31 + a32); a31 = -a31; a32 = -a32;
        }

        // apply the isometry.
        c += (a + b + f + g + h);
        f += (h + b + b); f = -f;
        g += (h + a + a); g = -g;
    }

    if (a == -h && g != 0)
    {
        if (saveIsometry)
        {
            // transform on the right by...
            // | -1 -1  0 |
            // |  0 -1  0 |
            // |  0  0  1 |
            a12 += a11; a12 = -a12; a11 = -a11;
            a22 += a21; a22 = -a22; a21 = -a21;
            a32 += a31; a32 = -a32; a31 = -a31;
        }

        // apply the isometry.
        f += g; f = -f;
        g = -g;
        h = -h;
    }

    if (a == -g && h != 0)
    {
        if (saveIsometry)
        {
            // transform on the right by...
            // | -1  0 -1 |
            // |  0  1  0 |
            // |  0  0 -1 |
            a13 += a11; a13 = -a13; a11 = -a11;
            a23 += a21; a23 = -a23; a21 = -a21;
            a33 += a31; a33 = -a33; a31 = -a31;
        }

        // apply the isometry.
        f += h; f = -f;
        h = -h;
        g += (a + a);
    }

    if (b == -f && h != 0)
    {
        if (saveIsometry)
        {
            // transform on the right by...
            // |  1  0  0 |
            // |  0 -1 -1 |
            // |  0  0 -1 |
            a13 += a12; a13 = -a13; a12 = -a12;
            a23 += a22; a23 = -a23; a22 = -a22;
            a33 += a32; a33 = -a33; a32 = -a32;
        }

        // apply the isometry.
        g += h; g = -g;
        h = -h;
        f += (b + b);
    }

    if (a == h && g > f + f)
    {
        if (saveIsometry)
        {
            // transform on the right by...
            // | -1 -1  0 |
            // |  0  1  0 |
            // |  0  0 -1 |
            a12 -= a11; a11 = -a11; a13 = -a13;
            a22 -= a21; a21 = -a21; a23 = -a23;
            a32 -= a31; a31 = -a31; a33 = -a33;
        }

        // apply the isometry.
        f = g - f;
    }

    if (a == g && h > f + f)
    {
        if (saveIsometry)
        {
            // transform on the right by...
            // | -1  0 -1 |
            // |  0 -1  0 |
            // |  0  0  1 |
            a13 -= a11; a11 = -a11; a12 = -a12;
            a23 -= a21; a21 = -a21; a22 = -a22;
            a33 -= a31; a31 = -a31; a32 = -a32;
        }

        // apply the isometry.
        f = h - f;
    }

    if (b == f && h > g + g)
    {
        if (saveIsometry)
        {
            // transform on the right by...
            // | -1  0  0 |
            // |  0 -1 -1 |
            // |  0  0  1 |
            a13 -= a12; a12 = -a12; a11 = -a11;
            a23 -= a22; a22 = -a22; a21 = -a21;
            a33 -= a32; a32 = -a32; a31 = -a31;
        }

        // apply the isometry.
        g = h - g;
    }

    if (a == b && abs(f) > abs(g))
    {
        if (saveIsometry)
        {
            // transform on the right by...
            // |  0 -1  0 |
            // | -1  0  0 |
            // |  0  0 -1 |
            temp = -a11; a11 = -a12; a12 = temp; a13 = -a13;
            temp = -a21; a21 = -a22; a22 = temp; a23 = -a23;
            temp = -a31; a31 = -a32; a32 = temp; a33 = -a33;
        }

        // apply the isometry.
        std::swap(a, b);
        std::swap(g, f);
    }

    if (b == c && abs(g) > abs(h))
    {
        if (saveIsometry)
        {
            // transform on the right by...
            // | -1  0  0 |
            // |  0  0 -1 |
            // |  0 -1  0 |
            temp = a12; a12 = -a13; a13 = -temp; a11 = -a11;
            temp = a22; a22 = -a23; a23 = -temp; a21 = -a21;
            temp = a32; a32 = -a33; a33 = -temp; a31 = -a31;
        }

        // apply the isometry.
        std::swap(g, h);
    }

    if (a == b && abs(f) > abs(g))
    {
        if (saveIsometry)
        {
            // transform on the right by...
            // |  0 -1  0 |
            // | -1  0  0 |
            // |  0  0 -1 |
            temp = a11; a11 = -a12; a12 = -temp; a13 = -a13;
            temp = a21; a21 = -a22; a22 = -temp; a23 = -a23;
            temp = a31; a31 = -a32; a32 = -temp; a33 = -a33;
        }

        // apply the isometry.
        std::swap(g, f);
    }

    // The reduced quadratic form.
    std::shared_ptr<QuadFormZZ> qq = std::make_shared<QuadFormZZ>(
        q.discriminant(),
        mpz_class(a),
        mpz_class(b),
        mpz_class(c),
        mpz_class(f),
        mpz_class(g),
        mpz_class(h));

    if (saveIsometry)
    {
        // The isometry accompanying this reduction.
        std::shared_ptr<IsometryQQ> s = std::make_shared<IsometryQQ>(
            mpq_class(a11),
            mpq_class(a12),
            mpq_class(a13),
            mpq_class(a21),
            mpq_class(a22),
            mpq_class(a23),
            mpq_class(a31),
            mpq_class(a32),
            mpq_class(a33));

        // Assign the isometry.
        qq->isometry(s);

#ifdef DEBUG
        assert( qq->isometry()->is_isometry(q, *qq) );
#endif
    }
    else
    {
        qq->isometry(std::make_shared<IsometryQQ>(true));
    }

    return qq;
}

template<>
std::shared_ptr<QuadFormZZ> QuadFormZZ::reduce(const QuadFormZZ& q,
                                               bool saveIsometry)
{
    // The coefficients.
    mpz_class a = q.a_;
    mpz_class b = q.b_;
    mpz_class c = q.c_;
    mpz_class f = q.f_;
    mpz_class g = q.g_;
    mpz_class h = q.h_;

    // The pointer that we will eventually return.
    std::shared_ptr<QuadFormZZ> qq;

    // If the size of all coefficients is reasonably small, we can perform all
    // computations with native int64_t data types instead of mpz_class
    // objects.
    if (abs(a) < UPPERBOUND && abs(b) < UPPERBOUND && abs(c) < UPPERBOUND &&
        abs(f) < UPPERBOUND && abs(g) < UPPERBOUND && abs(h) < UPPERBOUND)
    {
        // The coefficients as int64_t data types.
        int64_t aa = mpz_get_si(a.get_mpz_t());
        int64_t bb = mpz_get_si(b.get_mpz_t());
        int64_t cc = mpz_get_si(c.get_mpz_t());
        int64_t ff = mpz_get_si(f.get_mpz_t());
        int64_t gg = mpz_get_si(g.get_mpz_t());
        int64_t hh = mpz_get_si(h.get_mpz_t());

        // Reduce with T=int64_t.
        qq = reduceT<mpz_class, mpq_class, int64_t>(q, aa, bb, cc, ff, gg, hh, saveIsometry);
    }
    else
    {
        // Reduce with T=mpz_class.
        qq = reduceT<mpz_class, mpq_class, mpz_class>(q, a, b, c, f, g, h, saveIsometry);
    }

    // Set the reduction flag and return the shared pointer.
    qq->reduced_ = true;
    return qq;
}

template<>
std::vector<mpz_class> QuadFormZZ::isotropic_vector(const mpz_class& p) const
{
    std::vector<mpz_class> vec(3, 0);

    if (MathZZ::gcd(this->disc_, p) != 1)
    {
        throw std::runtime_error("Prime ideal is not coprime to discriminant.");
    }

    // We first reduce the coefficients modulo p and check for trivial
    // isotropic vectors.
    mpz_class a = this->a_ % p;
    if (a == 0)
    {
        vec[0] = 1;
        return vec;
    }
    if (a < 0) { a += p; }

    mpz_class b = this->b_ % p;
    if (b == 0)
    {
        vec[1] = 1;
        return vec;
    }
    if (b < 0) { b += p; }

    mpz_class c = this->c_ % p;
    if (c == 0)
    {
        vec[2] = 1;
        return vec;
    }
    if (c < 0) { c += p; }

    mpz_class f = this->f_ % p; if (f < 0) { f += p; }
    mpz_class g = this->g_ % p; if (g < 0) { g += p; }
    mpz_class h = this->h_ % p; if (h < 0) { h += p; }

    if (p == 2)
    {
        if (f == 0) { vec[1] = 1; vec[2] = 1; return vec; }
        if (g == 0) { vec[0] = 1; vec[2] = 1; return vec; }
        if (h == 0) { vec[0] = 1; vec[1] = 1; return vec; }

        // There is no return statement here, since this would require all
        // coefficients to be 1 mod 2. In this case, the discriminant must be
        // divisible by 2, invalidating the fact that p must be coprime to the
        // discriminant.
    }

    // A coefficient arising from diagonalizing the quadratic form.
    mpz_class v = (h*h-4*a*b) % p;

    if (v == 0)
    {
        vec[0] = (-h * MathZZ::modinv(2*a, p)) % p;
        vec[1] = 1;
        MathZZ::fix_vector(vec, p);
        return vec;
    }

    int64_t ll = MathZZ::kronecker(v, p);


    if (ll == 1)
    {
        mpz_class sqr = MathZZ::square_root(v, p);

        mpz_class temp = MathZZ::modinv(sqr, p);
        vec[0] = (1 - h*temp) % p;
        vec[1] = (2*a*temp) % p;

        if (vec[0] != 0)
        {
            mpz_class temp = MathZZ::modinv(vec[0], p);
            vec[0] = 1;
            vec[1] = (vec[1] * temp) % p;
        }
        else
        {
            vec[1] = 1;
        }

        MathZZ::fix_vector(vec, p);
        return vec;
    }

    mpz_class A = (2*a) % p;
    mpz_class B = ((-v) * MathZZ::modinv(2*a, p)) % p;
    mpz_class C = ((2 * this->disc_) * MathZZ::modinv(-v, p)) % p;

    mpz_class num, x, y;
    do
    {
        // Randomly choose values modulo p in the hopes of constructing a
        // local square represented by the quadratic form.
        x = rand() % p;
        y = rand() % (x == 0 ? p-1 : p);
        if (x == 0) { y++; }

        // The value we wish to test.
        num = -(((((A * x) % p) * x) % p + (((B * y) % p) * y) % p) *
                MathZZ::modinv(C, p)) % p;
    }
    while (num == 0 || MathZZ::kronecker(num, p) == -1);

    // Compute its square root modulo p.
    mpz_class sqr = MathZZ::square_root(num, p);

    // Avoid recomputing this value.
    mpz_class temp = MathZZ::modinv(-v, p);

    // Build the isotropic vector.
    vec[0] = (((x - (h * y * MathZZ::modinv(2*a, p))) % p) +
              sqr * ((f * h - 2 * b * g) % p) * temp) % p;
    vec[1] = (y + sqr * (g * h - 2 * a *f) * temp) % p;
    vec[2] = sqr;

    // Normalize the vector.
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

    MathZZ::fix_vector(vec, p);
    return vec;
}

template<>
int64_t QuadFormZZ::hash_value() const
{
    int64_t h = 0;
    h = ((h << 12) ^ (mpz_get_si(this->a_.get_mpz_t()) * 37));
    h = ((h << 12) ^ (mpz_get_si(this->b_.get_mpz_t()) * 101));
    h = ((h << 12) ^ (mpz_get_si(this->c_.get_mpz_t()) * 199));
    h = ((h << 12) ^ (mpz_get_si(this->f_.get_mpz_t()) * 683));
    h = ((h << 12) ^ (mpz_get_si(this->g_.get_mpz_t()) * 827));
    h = ((h << 12) ^ (mpz_get_si(this->h_.get_mpz_t()) * 941));
    return h;
}

// This code is adapted form Sage's codebase.
bool border(const QuadFormZZ& q, int64_t n)
{
    switch (n)
    {
        case 1:
            return (q.a() == q.h()) && (q.g() == q.f()*2);
        case 2:
            return (q.a() == q.g()) && (q.h() == q.f()*2);
        case 3:
            return (q.b() == q.f()) && (q.h() == q.g()*2);
        case 4:
            return (q.a() == -q.h());
        case 5:
            return (q.a() == -q.g());
        case 6:
            return (q.b() == -q.f());
        case 7:
            return (q.a() + q.b() + q.f() + q.g() + q.h() == 0) &&
                   (q.a()*2 + q.g()*2 + q.h() == 0);
        case 8:
            return (q.a() == q.b()) && (q.f() == q.g());
        case 9:
            return (q.b() == q.c()) && (q.g() == q.h());
        case 10:
            return (q.f() == q.g()) && (q.f() == 0);
        case 11:
            return (q.f() == q.h()) && (q.f() == 0);
        case 12:
            return (q.g() == q.h()) && (q.g() == 0);
        case 13:
            return (q.f() == q.g()) && (q.g() == q.h()) &&
                   (q.h() == q.a());
        case 14:
            return (q.a() == q.g()) && (q.a() == q.h());
        case 15:
            return (q.a() == q.b()) &&
                   (q.a() + q.b() + q.f() + q.g() + q.h() == 0);
        case 16:
            return (q.a() == q.b()) && (q.b() == q.c()) &&
                   (q.a() + q.b() + q.f() + q.g() + q.h() == 0);
        default:
            return false;
    }
}

template<>
const std::vector<IsometryQQPtr>& QuadFormZZ::automorphisms(void)
{
    if (!this->reduced_)
    {
        throw std::runtime_error("A quadratic form must be reduced before"
            " its automorphisms can be computed.");
    }

    if (this->auts_.size() > 0)
    {
        return this->auts_;
    }

    if (border(*this, 1))
    {
        if (border(*this, 2))
        {
            if (border(*this, 14))
            {
                if (border(*this, 9))
                {
                    this->numAuts_ = 16;
                    this->auts_.push_back(AutomorphismZZ::A100010001);
                    this->auts_.push_back(AutomorphismZZ::ANNN001010);
                    this->auts_.push_back(AutomorphismZZ::ANN001000N);
                    this->auts_.push_back(AutomorphismZZ::AN0N0N0001);
                    this->auts_.push_back(AutomorphismZZ::AN0000N0N0);
                    this->auts_.push_back(AutomorphismZZ::A10100N010);
                    this->auts_.push_back(AutomorphismZZ::A1100010N0);
                    this->auts_.push_back(AutomorphismZZ::A1110N000N);
                    return this->auts_;
                }
                else
                {
                    this->numAuts_ = 8;
                    this->auts_.push_back(AutomorphismZZ::A100010001);
                    this->auts_.push_back(AutomorphismZZ::ANN001000N);
                    this->auts_.push_back(AutomorphismZZ::AN0N0N0001);
                    this->auts_.push_back(AutomorphismZZ::A1110N000N);
                    return this->auts_;
                }
            }
        }
        else
        {
            this->numAuts_ = 4;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::ANN001000N);
            return this->auts_;
        }
    }

    if (border(*this, 2))
    {
        this->numAuts_ = 4;
        this->auts_.push_back(AutomorphismZZ::A100010001);
        this->auts_.push_back(AutomorphismZZ::AN0N0N0001);
        return this->auts_;
    }

    if (border(*this, 3))
    {
        this->numAuts_ = 4;
        this->auts_.push_back(AutomorphismZZ::A100010001);
        this->auts_.push_back(AutomorphismZZ::AN000NN001);
        return this->auts_;
    }

    if (border(*this, 4))
    {
        if (border(*this, 10))
        {
            if (border(*this, 8))
            {
                this->numAuts_ = 24;
                this->auts_.push_back(AutomorphismZZ::A100010001);
                this->auts_.push_back(AutomorphismZZ::AN00N1000N);
                this->auts_.push_back(AutomorphismZZ::AN000N0001);
                this->auts_.push_back(AutomorphismZZ::AN10N00001);
                this->auts_.push_back(AutomorphismZZ::AN1001000N);
                this->auts_.push_back(AutomorphismZZ::A0N0N0000N);
                this->auts_.push_back(AutomorphismZZ::A0N01N0001);
                this->auts_.push_back(AutomorphismZZ::A010N10001);
                this->auts_.push_back(AutomorphismZZ::A01010000N);
                this->auts_.push_back(AutomorphismZZ::A1N00N000N);
                this->auts_.push_back(AutomorphismZZ::A1N0100001);
                this->auts_.push_back(AutomorphismZZ::A1001N000N);
                return this->auts_;
            }
            else
            {
                this->numAuts_ = 8;
                this->auts_.push_back(AutomorphismZZ::A100010001);
                this->auts_.push_back(AutomorphismZZ::AN000N0001);
                this->auts_.push_back(AutomorphismZZ::AN1001000N);
                this->auts_.push_back(AutomorphismZZ::A1N00N000N);
                return this->auts_;
            }
        }
        else
        {
            this->numAuts_ = 4;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::A1N00N000N);
            return this->auts_;
        }
    }

    if (border(*this, 5))
    {
        if (border(*this, 6))
        {
            if (border(*this, 7))
            {
                if (border(*this, 8))
                {
                    if (border(*this, 15))
                    {
                        this->numAuts_ = 16;
                        this->auts_.push_back(AutomorphismZZ::A100010001);
                        this->auts_.push_back(AutomorphismZZ::AN0001N00N);
                        this->auts_.push_back(AutomorphismZZ::AN010N1001);
                        this->auts_.push_back(AutomorphismZZ::A0N0N0000N);
                        this->auts_.push_back(AutomorphismZZ::A0N1100001);
                        this->auts_.push_back(AutomorphismZZ::A01N10N00N);
                        this->auts_.push_back(AutomorphismZZ::A010N01001);
                        this->auts_.push_back(AutomorphismZZ::A10N0N000N);
                        return this->auts_;
                    }
                }
                else
                {
                    this->numAuts_ = 8;
                    this->auts_.push_back(AutomorphismZZ::A100010001);
                    this->auts_.push_back(AutomorphismZZ::AN0001N00N);
                    this->auts_.push_back(AutomorphismZZ::AN010N1001);
                    this->auts_.push_back(AutomorphismZZ::A10N0N000N);
                    return this->auts_;
                }
            }
        }
        else if (border(*this, 11))
        {
            this->numAuts_ = 8;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::AN0001000N);
            this->auts_.push_back(AutomorphismZZ::AN010N0001);
            this->auts_.push_back(AutomorphismZZ::A10N0N000N);
            return this->auts_;
        }
        else
        {
            this->numAuts_ = 4;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::A10N0N000N);
            return this->auts_;
        }
    }

    if (border(*this, 6))
    {
        if (border(*this, 12))
        {
            if (border(*this, 9))
            {
                this->numAuts_ = 24;
                this->auts_.push_back(AutomorphismZZ::A100010001);
                this->auts_.push_back(AutomorphismZZ::AN000N00N1);
                this->auts_.push_back(AutomorphismZZ::AN000N1001);
                this->auts_.push_back(AutomorphismZZ::AN0000N0N0);
                this->auts_.push_back(AutomorphismZZ::AN00001010);
                this->auts_.push_back(AutomorphismZZ::AN0001N00N);
                this->auts_.push_back(AutomorphismZZ::AN0001001N);
                this->auts_.push_back(AutomorphismZZ::A1000N000N);
                this->auts_.push_back(AutomorphismZZ::A1000N10N0);
                this->auts_.push_back(AutomorphismZZ::A10000N01N);
                this->auts_.push_back(AutomorphismZZ::A1000010N1);
                this->auts_.push_back(AutomorphismZZ::A10001N010);
                return this->auts_;
            }
            else
            {
                this->numAuts_ = 8;
                this->auts_.push_back(AutomorphismZZ::A100010001);
                this->auts_.push_back(AutomorphismZZ::AN000N1001);
                this->auts_.push_back(AutomorphismZZ::AN0001N00N);
                this->auts_.push_back(AutomorphismZZ::A1000N000N);
                return this->auts_;
            }
        }
        else
        {
            this->numAuts_ = 4;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::AN0001N00N);
            return this->auts_;
        }
    }

    if (border(*this, 7))
    {
        if (border(*this, 8) && border(*this, 15))
        {
            if (border(*this, 16))
            {
                if (border(*this, 9))
                {
                    this->numAuts_ = 48;
                    this->auts_.push_back(AutomorphismZZ::A100010001);
                    this->auts_.push_back(AutomorphismZZ::AN00N01N10);
                    this->auts_.push_back(AutomorphismZZ::AN0000N0N0);
                    this->auts_.push_back(AutomorphismZZ::AN01N10N00);
                    this->auts_.push_back(AutomorphismZZ::AN010N1001);
                    this->auts_.push_back(AutomorphismZZ::AN10N00N01);
                    this->auts_.push_back(AutomorphismZZ::AN1001001N);
                    this->auts_.push_back(AutomorphismZZ::A0N0N0000N);
                    this->auts_.push_back(AutomorphismZZ::A0N01N00N1);
                    this->auts_.push_back(AutomorphismZZ::A0N10N01N0);
                    this->auts_.push_back(AutomorphismZZ::A0N1001N01);
                    this->auts_.push_back(AutomorphismZZ::A00N0N0N00);
                    this->auts_.push_back(AutomorphismZZ::A00N01N10N);
                    this->auts_.push_back(AutomorphismZZ::A001N010N1);
                    this->auts_.push_back(AutomorphismZZ::A001100010);
                    this->auts_.push_back(AutomorphismZZ::A01NN10010);
                    this->auts_.push_back(AutomorphismZZ::A01N10N00N);
                    this->auts_.push_back(AutomorphismZZ::A010001100);
                    this->auts_.push_back(AutomorphismZZ::A01001NN10);
                    this->auts_.push_back(AutomorphismZZ::A1N00N10N0);
                    this->auts_.push_back(AutomorphismZZ::A1N010N100);
                    this->auts_.push_back(AutomorphismZZ::A10N00N01N);
                    this->auts_.push_back(AutomorphismZZ::A10N1001N0);
                    this->auts_.push_back(AutomorphismZZ::A1001N010N);
                    return this->auts_;
                }
                else
                {
                    this->numAuts_ = 16;
                    this->auts_.push_back(AutomorphismZZ::A100010001);
                    this->auts_.push_back(AutomorphismZZ::AN00N01N10);
                    this->auts_.push_back(AutomorphismZZ::AN010N1001);
                    this->auts_.push_back(AutomorphismZZ::A0N0N0000N);
                    this->auts_.push_back(AutomorphismZZ::A0N10N01N0);
                    this->auts_.push_back(AutomorphismZZ::A01N10N00N);
                    this->auts_.push_back(AutomorphismZZ::A01001NN10);
                    this->auts_.push_back(AutomorphismZZ::A10N1001N0);
                    return this->auts_;
                }
            }
            else
            {
                this->numAuts_ = 8;
                this->auts_.push_back(AutomorphismZZ::A100010001);
                this->auts_.push_back(AutomorphismZZ::AN010N1001);
                this->auts_.push_back(AutomorphismZZ::A0N0N0000N);
                this->auts_.push_back(AutomorphismZZ::A01N10N00N);
                return this->auts_;
            }
        }
        else if (border(*this, 9))
        {
            this->numAuts_ = 12;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::AN0000N0N0);
            this->auts_.push_back(AutomorphismZZ::AN010N1001);
            this->auts_.push_back(AutomorphismZZ::AN1001001N);
            this->auts_.push_back(AutomorphismZZ::A1N00N10N0);
            this->auts_.push_back(AutomorphismZZ::A10N00N01N);
            return this->auts_;
        }
        else
        {
            this->numAuts_ = 4;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::AN010N1001);
            return this->auts_;
        }
    }

    if (border(*this, 8))
    {
        if (border(*this, 9))
        {
            if (border(*this, 10) && border(*this, 11) && border(*this, 12))
            {
                this->numAuts_ = 48;
                this->auts_.push_back(AutomorphismZZ::A100010001);
                this->auts_.push_back(AutomorphismZZ::AN000N0001);
                this->auts_.push_back(AutomorphismZZ::AN0000N0N0);
                this->auts_.push_back(AutomorphismZZ::AN00001010);
                this->auts_.push_back(AutomorphismZZ::AN0001000N);
                this->auts_.push_back(AutomorphismZZ::A0N0N0000N);
                this->auts_.push_back(AutomorphismZZ::A0N000N100);
                this->auts_.push_back(AutomorphismZZ::A0N0001N00);
                this->auts_.push_back(AutomorphismZZ::A0N0100001);
                this->auts_.push_back(AutomorphismZZ::A00NN00010);
                this->auts_.push_back(AutomorphismZZ::A00N0N0N00);
                this->auts_.push_back(AutomorphismZZ::A00N010100);
                this->auts_.push_back(AutomorphismZZ::A00N1000N0);
                this->auts_.push_back(AutomorphismZZ::A001N000N0);
                this->auts_.push_back(AutomorphismZZ::A0010N0100);
                this->auts_.push_back(AutomorphismZZ::A001010N00);
                this->auts_.push_back(AutomorphismZZ::A001100010);
                this->auts_.push_back(AutomorphismZZ::A010N00001);
                this->auts_.push_back(AutomorphismZZ::A01000NN00);
                this->auts_.push_back(AutomorphismZZ::A010001100);
                this->auts_.push_back(AutomorphismZZ::A01010000N);
                this->auts_.push_back(AutomorphismZZ::A1000N000N);
                this->auts_.push_back(AutomorphismZZ::A10000N010);
                this->auts_.push_back(AutomorphismZZ::A1000010N0);
                return this->auts_;
            }
            else if (border(*this, 13) && border(*this, 14))
            {
                this->numAuts_ = 48;
                this->auts_.push_back(AutomorphismZZ::A100010001);
                this->auts_.push_back(AutomorphismZZ::ANNN001010);
                this->auts_.push_back(AutomorphismZZ::ANNN010100);
                this->auts_.push_back(AutomorphismZZ::ANNN100001);
                this->auts_.push_back(AutomorphismZZ::AN000N0111);
                this->auts_.push_back(AutomorphismZZ::AN0000N0N0);
                this->auts_.push_back(AutomorphismZZ::AN0011100N);
                this->auts_.push_back(AutomorphismZZ::A0N0N0000N);
                this->auts_.push_back(AutomorphismZZ::A0N000N111);
                this->auts_.push_back(AutomorphismZZ::A0N0111N00);
                this->auts_.push_back(AutomorphismZZ::A00NN00111);
                this->auts_.push_back(AutomorphismZZ::A00N0N0N00);
                this->auts_.push_back(AutomorphismZZ::A00N1110N0);
                this->auts_.push_back(AutomorphismZZ::A001NNN100);
                this->auts_.push_back(AutomorphismZZ::A001010NNN);
                this->auts_.push_back(AutomorphismZZ::A001100010);
                this->auts_.push_back(AutomorphismZZ::A010NNN001);
                this->auts_.push_back(AutomorphismZZ::A010001100);
                this->auts_.push_back(AutomorphismZZ::A010100NNN);
                this->auts_.push_back(AutomorphismZZ::A100NNN010);
                this->auts_.push_back(AutomorphismZZ::A100001NNN);
                this->auts_.push_back(AutomorphismZZ::A111N000N0);
                this->auts_.push_back(AutomorphismZZ::A1110N000N);
                this->auts_.push_back(AutomorphismZZ::A11100NN00);
                return this->auts_;
            }
            else
            {
                this->numAuts_ = 12;
                this->auts_.push_back(AutomorphismZZ::A100010001);
                this->auts_.push_back(AutomorphismZZ::AN0000N0N0);
                this->auts_.push_back(AutomorphismZZ::A0N0N0000N);
                this->auts_.push_back(AutomorphismZZ::A00N0N0N00);
                this->auts_.push_back(AutomorphismZZ::A001100010);
                this->auts_.push_back(AutomorphismZZ::A010001100);
                return this->auts_;
            }
        }
        else if (border(*this, 10))
        {
            if (border(*this, 11) && border(*this, 12))
            {
                this->numAuts_ = 16;
                this->auts_.push_back(AutomorphismZZ::A100010001);
                this->auts_.push_back(AutomorphismZZ::AN000N0001);
                this->auts_.push_back(AutomorphismZZ::AN0001000N);
                this->auts_.push_back(AutomorphismZZ::A0N0N0000N);
                this->auts_.push_back(AutomorphismZZ::A0N0100001);
                this->auts_.push_back(AutomorphismZZ::A010N00001);
                this->auts_.push_back(AutomorphismZZ::A01010000N);
                this->auts_.push_back(AutomorphismZZ::A1000N000N);
                return this->auts_;
            }
            else
            {
                this->numAuts_ = 8;
                this->auts_.push_back(AutomorphismZZ::A100010001);
                this->auts_.push_back(AutomorphismZZ::AN000N0001);
                this->auts_.push_back(AutomorphismZZ::A0N0N0000N);
                this->auts_.push_back(AutomorphismZZ::A01010000N);
                return this->auts_;
            }
        }
        else if (border(*this, 14))
        {
            this->numAuts_ = 12;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::ANNN100001);
            this->auts_.push_back(AutomorphismZZ::AN0011100N);
            this->auts_.push_back(AutomorphismZZ::A0N0N0000N);
            this->auts_.push_back(AutomorphismZZ::A010NNN001);
            this->auts_.push_back(AutomorphismZZ::A1110N000N);
            return this->auts_;
        }
        else
        {
            this->numAuts_ = 4;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::A0N0N0000N);
            return this->auts_;
        }
    }

    if (border(*this, 9))
    {
        if (border(*this, 12))
        {
            if (border(*this, 10) && border(*this, 11))
            {
                this->numAuts_ = 16;
                this->auts_.push_back(AutomorphismZZ::A100010001);
                this->auts_.push_back(AutomorphismZZ::AN000N0001);
                this->auts_.push_back(AutomorphismZZ::AN0000N0N0);
                this->auts_.push_back(AutomorphismZZ::AN00001010);
                this->auts_.push_back(AutomorphismZZ::AN0001000N);
                this->auts_.push_back(AutomorphismZZ::A1000N000N);
                this->auts_.push_back(AutomorphismZZ::A10000N010);
                this->auts_.push_back(AutomorphismZZ::A1000010N0);
                return this->auts_;
            }
            else
            {
                this->numAuts_ = 8;
                this->auts_.push_back(AutomorphismZZ::A100010001);
                this->auts_.push_back(AutomorphismZZ::AN0000N0N0);
                this->auts_.push_back(AutomorphismZZ::AN00001010);
                this->auts_.push_back(AutomorphismZZ::A1000N000N);
                return this->auts_;
            }
        }
        else if (border(*this, 14))
        {
            if (border(*this, 13))
            {
                this->numAuts_ = 8;
                this->auts_.push_back(AutomorphismZZ::A100010001);
                this->auts_.push_back(AutomorphismZZ::ANNN001010);
                this->auts_.push_back(AutomorphismZZ::AN0000N0N0);
                this->auts_.push_back(AutomorphismZZ::A1110N000N);
                return this->auts_;
            }
            else
            {
                this->numAuts_ = 8;
                this->auts_.push_back(AutomorphismZZ::A100010001);
                this->auts_.push_back(AutomorphismZZ::ANNN001010);
                this->auts_.push_back(AutomorphismZZ::AN0000N0N0);
                this->auts_.push_back(AutomorphismZZ::A1110N000N);
                return this->auts_;
            }
        }
        else if (border(*this, 15))
        {
            this->numAuts_ = 16;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::AN00N01N10);
            this->auts_.push_back(AutomorphismZZ::AN0000N0N0);
            this->auts_.push_back(AutomorphismZZ::A0N10N01N0);
            this->auts_.push_back(AutomorphismZZ::A0N1001N01);
            this->auts_.push_back(AutomorphismZZ::A01NN10010);
            this->auts_.push_back(AutomorphismZZ::A01N10N00N);
            this->auts_.push_back(AutomorphismZZ::A1001N010N);
            return this->auts_;
        }
        else
        {
            this->numAuts_ = 4;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::AN0000N0N0);
            return this->auts_;
        }
    }

    if (border(*this, 10))
    {
        if (border(*this, 11) && border(*this, 12))
        {
            this->numAuts_ = 8;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::AN000N0001);
            this->auts_.push_back(AutomorphismZZ::AN0001000N);
            this->auts_.push_back(AutomorphismZZ::A1000N000N);
            return this->auts_;
        }
        else
        {
            this->numAuts_ = 4;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::AN000N0001);
            return this->auts_;
        }
    }

    if (border(*this, 11))
    {
        this->numAuts_ = 4;
        this->auts_.push_back(AutomorphismZZ::A100010001);
        this->auts_.push_back(AutomorphismZZ::AN0001000N);
        return this->auts_;
    }

    if (border(*this, 12))
    {
        this->numAuts_ = 4;
        this->auts_.push_back(AutomorphismZZ::A100010001);
        this->auts_.push_back(AutomorphismZZ::A1000N000N);
        return this->auts_;
    }

    if (border(*this, 13) && border(*this, 14))
    {
        this->numAuts_ = 4;
        this->auts_.push_back(AutomorphismZZ::A100010001);
        this->auts_.push_back(AutomorphismZZ::A1110N000N);
        return this->auts_;
    }

    if (border(*this, 14))
    {
        this->numAuts_ = 4;
        this->auts_.push_back(AutomorphismZZ::A100010001);
        this->auts_.push_back(AutomorphismZZ::A1110N000N);
        return this->auts_;
    }

    if (border(*this, 15))
    {
        if (border(*this, 16))
        {
            this->numAuts_ = 8;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::AN00N01N10);
            this->auts_.push_back(AutomorphismZZ::A0N10N01N0);
            this->auts_.push_back(AutomorphismZZ::A01N10N00N);
            return this->auts_;
        }
        else
        {
            this->numAuts_ = 4;
            this->auts_.push_back(AutomorphismZZ::A100010001);
            this->auts_.push_back(AutomorphismZZ::A01N10N00N);
            return this->auts_;
        }
    }

    this->numAuts_ = 2;
    this->auts_.push_back(AutomorphismZZ::A100010001);
    return this->auts_;
}
