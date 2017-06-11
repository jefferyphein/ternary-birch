#include <gmpxx.h>
#include <cassert>
#include "QuadForm.h"
#include "Math.h"
#include "AutomorphismZZ.h"

typedef QuadForm<mpz_class, mpq_class> QuadFormZZ;
typedef Isometry<mpz_class, mpq_class> IsometryQQ;
typedef std::shared_ptr<IsometryQQ> IsometryQQPtr;
typedef Math<mpz_class, mpq_class> MathZZ;

template<>
std::shared_ptr<QuadFormZZ> QuadFormZZ::reduce(const QuadFormZZ& q,
                                               bool saveIsometry)
{
    mpz_class a = q.a();
    mpz_class b = q.b();
    mpz_class c = q.c();
    mpz_class f = q.f();
    mpz_class g = q.g();
    mpz_class h = q.h();

    auto s = std::make_shared<IsometryQQ>(true);

    // flag controlling the initial reduction loop.
    bool flag = true;
    
    // temporary variable.
    mpq_class temp;
    
    while (flag)
    {
        mpz_class t = a + b + f + g + h;
        if (t < 0)
        {
            if (saveIsometry)
            {
                // transform on the right by...
                // |  1  0  1 |
                // |  0  1  1 |
                // |  0  0  1 |
                s->a13 += (s->a11 + s->a12);
                s->a23 += (s->a21 + s->a22);
                s->a33 += (s->a31 + s->a32);
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
                s->a12 += t * s->a11;
                s->a22 += t * s->a21;
                s->a32 += t * s->a31;
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
                s->a13 += t * s->a12;
                s->a23 += t * s->a22;
                s->a33 += t * s->a32;
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
                s->a13 += t * s->a11;
                s->a23 += t * s->a21;
                s->a33 += t * s->a31;
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
                temp = s->a11; s->a11 = -s->a12; s->a12 = -temp; s->a13 = -s->a13;
                temp = s->a21; s->a21 = -s->a22; s->a22 = -temp; s->a23 = -s->a23;
                temp = s->a31; s->a31 = -s->a32; s->a32 = -temp; s->a33 = -s->a33;
            }
            
            // apply the isometry.
            a.swap(b);
            f.swap(g);
        }
        
        if (b > c || (b == c && abs(g) > abs(h)))
        {
            if (saveIsometry)
            {
                // transform on the right by...
                // | -1  0  0 |
                // |  0  0 -1 |
                // |  0 -1  0 |
                temp = s->a12; s->a12 = -s->a13; s->a13 = -temp; s->a11 = -s->a11;
                temp = s->a22; s->a22 = -s->a23; s->a23 = -temp; s->a21 = -s->a21;
                temp = s->a32; s->a32 = -s->a33; s->a33 = -temp; s->a31 = -s->a31;
            }
            
            // apply the isometry.
            b.swap(c);
            g.swap(h);
        }
        
        if (a > b || (a == b && abs(f) > abs(g)))
        {
            if (saveIsometry)
            {
                // transform on the right by...
                // |  0 -1  0 |
                // | -1  0  0 |
                // |  0  0 -1 |
                temp = s->a11; s->a11 = -s->a12; s->a12 = -temp; s->a13 = -s->a13;
                temp = s->a21; s->a21 = -s->a22; s->a22 = -temp; s->a23 = -s->a23;
                temp = s->a31; s->a31 = -s->a32; s->a32 = -temp; s->a33 = -s->a33;
            }
            
            // apply the isometry.
            a.swap(b);
            f.swap(g);
        }
        
        if (f * g * h > 0)
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
                    s->a11 = -s->a11;
                    s->a21 = -s->a21;
                    s->a31 = -s->a31;
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
                    s->a12 = -s->a12;
                    s->a22 = -s->a22;
                    s->a32 = -s->a32;
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
                    s->a13 = -s->a13;
                    s->a23 = -s->a23;
                    s->a33 = -s->a33;
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
                    s->a11 = -s->a11;
                    s->a21 = -s->a21;
                    s->a31 = -s->a31;
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
                    s->a12 = -s->a12;
                    s->a22 = -s->a22;
                    s->a32 = -s->a32;
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
                    s->a13 = -s->a13;
                    s->a23 = -s->a23;
                    s->a33 = -s->a33;
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
            s->a13 += (s->a11 + s->a12); s->a11 = -s->a11; s->a12 = -s->a12;
            s->a23 += (s->a21 + s->a22); s->a21 = -s->a21; s->a22 = -s->a22;
            s->a33 += (s->a31 + s->a32); s->a31 = -s->a31; s->a32 = -s->a32;
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
            s->a12 += s->a11; s->a12 = -s->a12; s->a11 = -s->a11;
            s->a22 += s->a21; s->a22 = -s->a22; s->a21 = -s->a21;
            s->a32 += s->a31; s->a32 = -s->a32; s->a31 = -s->a31;
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
            s->a13 += s->a11; s->a13 = -s->a13; s->a11 = -s->a11;
            s->a23 += s->a21; s->a23 = -s->a23; s->a21 = -s->a21;
            s->a33 += s->a31; s->a33 = -s->a33; s->a31 = -s->a31;
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
            s->a13 += s->a12; s->a13 = -s->a13; s->a12 = -s->a12;
            s->a23 += s->a22; s->a23 = -s->a23; s->a22 = -s->a22;
            s->a33 += s->a32; s->a33 = -s->a33; s->a32 = -s->a32;
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
            s->a12 -= s->a11; s->a11 = -s->a11; s->a13 = -s->a13;
            s->a22 -= s->a21; s->a21 = -s->a21; s->a23 = -s->a23;
            s->a32 -= s->a31; s->a31 = -s->a31; s->a33 = -s->a33;
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
            s->a13 -= s->a11; s->a11 = -s->a11; s->a12 = -s->a12;
            s->a23 -= s->a21; s->a21 = -s->a21; s->a22 = -s->a22;
            s->a33 -= s->a31; s->a31 = -s->a31; s->a32 = -s->a32;
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
            s->a13 -= s->a12; s->a12 = -s->a12; s->a11 = -s->a11;
            s->a23 -= s->a22; s->a22 = -s->a22; s->a21 = -s->a21;
            s->a33 -= s->a32; s->a32 = -s->a32; s->a31 = -s->a31;
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
            temp = -s->a11; s->a11 = -s->a12; s->a12 = temp; s->a13 = -s->a13;
            temp = -s->a21; s->a21 = -s->a22; s->a22 = temp; s->a23 = -s->a23;
            temp = -s->a31; s->a31 = -s->a32; s->a32 = temp; s->a33 = -s->a33;
        }
        
        // apply the isometry.
        a.swap(b);
        g.swap(f);
    }
    
    if (b == c && abs(g) > abs(h))
    {
        if (saveIsometry)
        {
            // transform on the right by...
            // | -1  0  0 |
            // |  0  0 -1 |
            // |  0 -1  0 |
            temp = s->a12; s->a12 = -s->a13; s->a13 = -temp; s->a11 = -s->a11;
            temp = s->a22; s->a22 = -s->a23; s->a23 = -temp; s->a21 = -s->a21;
            temp = s->a32; s->a32 = -s->a33; s->a33 = -temp; s->a31 = -s->a31;
        }
        
        // apply the isometry.
        g.swap(h);
    }
    
    if (a == b && abs(f) > abs(g))
    {
        if (saveIsometry)
        {
            // transform on the right by...
            // |  0 -1  0 |
            // | -1  0  0 |
            // |  0  0 -1 |
            temp = s->a11; s->a11 = -s->a12; s->a12 = -temp; s->a13 = -s->a13;
            temp = s->a21; s->a21 = -s->a22; s->a22 = -temp; s->a23 = -s->a23;
            temp = s->a31; s->a31 = -s->a32; s->a32 = -temp; s->a33 = -s->a33;
        }
        
        // apply the isometry.
        g.swap(f);
printf("t\n");
    }

#ifdef DEBUG
    if (saveIsometry)
    {
        assert( s->is_isometry(q, a, b, c, f, g, h) );
    }
#endif

    auto qq = std::make_shared<QuadFormZZ>(q.disc_, a, b, c, f, g, h);

    qq->s_ = saveIsometry ? s : std::make_shared<IsometryQQ>(true);
    qq->reduced_ = true;

#ifdef DEBUG
    if (saveIsometry)
    {
        assert( qq->isometry()->is_isometry(q, *qq) );
    }
#endif

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
