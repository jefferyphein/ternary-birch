#include <gmpxx.h>
#include <cassert>
#include "QuadForm.h"
#include "Prime.h"
#include "Math.h"

typedef QuadForm<mpz_class, mpq_class> QuadFormZZ;
typedef Isometry<mpz_class, mpq_class> IsometryQQ;
typedef Prime<mpz_class, mpq_class> PrimeZZ;

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

    auto s = std::shared_ptr<IsometryQQ>(
        new IsometryQQ(true)
    );

    // flag controlling the initial reduction loop.
    bool flag = true;
    
    // temporary variable.
    mpq_class temp;
    
    while (flag)
    {
        mpz_class t = a + b + f + g + h;
        if (t < 0)
        {
            // transform on the right by...
            // |  1  0  1 |
            // |  0  1  1 |
            // |  0  0  1 |
            s->a13 += (s->a11 + s->a12);
            s->a23 += (s->a21 + s->a22);
            s->a33 += (s->a31 + s->a32);
            
            // apply the isometry.
            c += t;
            f += (h + b + b);
            g += (h + a + a);
        }
        
        t = (a-h) / (a+a);
        if (a < h && (a-h) % (a+a) != 0) { t--; }
        if (t != 0)
        {
            // transform on the right by...
            // |  1  t  0 |
            // |  0  1  0 |
            // |  0  0  1 |
            s->a12 += t * s->a11;
            s->a22 += t * s->a21;
            s->a32 += t * s->a31;

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
            // transform on the right by...
            // |  1  0  0 |
            // |  0  1  t |
            // |  0  0  1 |
            s->a13 += t * s->a12;
            s->a23 += t * s->a22;
            s->a33 += t * s->a32;

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
            // transform on the right by...
            // |  1  0  t |
            // |  0  1  0 |
            // |  0  0  1 |
            s->a13 += t * s->a11;
            s->a23 += t * s->a21;
            s->a33 += t * s->a31;

            // apply the isometry.
            f += h * t;
            c += t * (a * t + g);
            g += t * (a + a);
        }
        
        if (a > b || (a == b && abs(f) > abs(g)))
        {
            // transform on the right by...
            // |  0 -1  0 |
            // | -1  0  0 |
            // |  0  0 -1 |
            temp = s->a11; s->a11 = -s->a12; s->a12 = -temp; s->a13 = -s->a13;
            temp = s->a21; s->a21 = -s->a22; s->a22 = -temp; s->a23 = -s->a23;
            temp = s->a31; s->a31 = -s->a32; s->a32 = -temp; s->a33 = -s->a33;
            
            // apply the isometry.
            a.swap(b);
            f.swap(g);
        }
        
        if (b > c || (b == c && abs(g) > abs(h)))
        {
            // transform on the right by...
            // | -1  0  0 |
            // |  0  0 -1 |
            // |  0 -1  0 |
            temp = s->a12; s->a12 = -s->a13; s->a13 = -temp; s->a11 = -s->a11;
            temp = s->a22; s->a22 = -s->a23; s->a23 = -temp; s->a21 = -s->a21;
            temp = s->a32; s->a32 = -s->a33; s->a33 = -temp; s->a31 = -s->a31;
            
            // apply the isometry.
            b.swap(c);
            g.swap(h);
        }
        
        if (a > b || (a == b && abs(f) > abs(g)))
        {
            // transform on the right by...
            // |  0 -1  0 |
            // | -1  0  0 |
            // |  0  0 -1 |
            temp = s->a11; s->a11 = -s->a12; s->a12 = -temp; s->a13 = -s->a13;
            temp = s->a21; s->a21 = -s->a22; s->a22 = -temp; s->a23 = -s->a23;
            temp = s->a31; s->a31 = -s->a32; s->a32 = -temp; s->a33 = -s->a33;
            
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
                // transform on the right by...
                // | -1  0  0 |
                // |  0  1  0 |
                // |  0  0  1 |
                s->a11 = -s->a11;
                s->a21 = -s->a21;
                s->a31 = -s->a31;
                
                // negate f.
                f = -f;
            }
            if (g < 0)
            {
                // transform on the right by...
                // |  1  0  0 |
                // |  0 -1  0 |
                // |  0  0  1 |
                s->a12 = -s->a12;
                s->a22 = -s->a22;
                s->a32 = -s->a32;
                
                // negate g.
                g = -g;
            }
            if (h < 0)
            {
                // transform on the right by...
                // |  1  0  0 |
                // |  0  1  0 |
                // |  0  0 -1 |
                s->a13 = -s->a13;
                s->a23 = -s->a23;
                s->a33 = -s->a33;
                
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
                // transform on the right by...
                // | -1  0  0 |
                // |  0  1  0 |
                // |  0  0  1 |
                s->a11 = -s->a11;
                s->a21 = -s->a21;
                s->a31 = -s->a31;
                
                // negate f.
                f = -f;
            }
            
            if (s2 == 1)
            {
                // transform on the right by...
                // |  1  0  0 |
                // |  0 -1  0 |
                // |  0  0  1 |
                s->a12 = -s->a12;
                s->a22 = -s->a22;
                s->a32 = -s->a32;
                
                // negate g.
                g = -g;
            }
            
            if (s3 == 1)
            {
                // transform on the right by...
                // |  1  0  0 |
                // |  0  1  0 |
                // |  0  0 -1 |
                s->a13 = -s->a13;
                s->a23 = -s->a23;
                s->a33 = -s->a33;
                
                // negate h.
                h = -h;
            }
        }
        
        flag = !(abs(f) <= b && abs(g) <= a && abs(h) <= a && a+b+f+g+h >= 0);
    }
    
    
    if (a+b+f+g+h == 0 && a+a+g+g+h > 0)
    {
        // transform on the right by...
        // | -1  0  1 |
        // |  0 -1  1 |
        // |  0  0  1 |
        s->a13 += (s->a11 + s->a12); s->a11 = -s->a11; s->a12 = -s->a12;
        s->a23 += (s->a21 + s->a22); s->a21 = -s->a21; s->a22 = -s->a22;
        s->a33 += (s->a31 + s->a32); s->a31 = -s->a31; s->a32 = -s->a32;
        
        // apply the isometry.
        c += (a + b + f + g + h);
        f += (h + b + b); f = -f;
        g += (h + a + a); g = -g;
    }
    
    if (a == -h && g != 0)
    {
        // transform on the right by...
        // | -1 -1  0 |
        // |  0 -1  0 |
        // |  0  0  1 |
        s->a12 += s->a11; s->a12 = -s->a12; s->a11 = -s->a11;
        s->a22 += s->a21; s->a22 = -s->a22; s->a21 = -s->a21;
        s->a32 += s->a31; s->a32 = -s->a32; s->a31 = -s->a31;
        
        // apply the isometry.
        f += g; f = -f;
        g = -g;
        h = -h;
    }
    
    if (a == -g && h != 0)
    {
        // transform on the right by...
        // | -1  0 -1 |
        // |  0  1  0 |
        // |  0  0 -1 |
        s->a13 += s->a11; s->a13 = -s->a13; s->a11 = -s->a11;
        s->a23 += s->a21; s->a23 = -s->a23; s->a21 = -s->a21;
        s->a33 += s->a31; s->a33 = -s->a33; s->a31 = -s->a31;
        
        // apply the isometry.
        f += h; f = -f;
        h = -h;
        g += (a + a);
    }
    
    if (b == -f && h != 0)
    {
        // transform on the right by...
        // |  1  0  0 |
        // |  0 -1 -1 |
        // |  0  0 -1 |
        s->a13 += s->a12; s->a13 = -s->a13; s->a12 = -s->a12;
        s->a23 += s->a22; s->a23 = -s->a23; s->a22 = -s->a22;
        s->a33 += s->a32; s->a33 = -s->a33; s->a32 = -s->a32;
        
        // apply the isometry.
        g += h; g = -g;
        h = -h;
        f += (b + b);
    }
    
    if (a == h && g > f + f)
    {
        // transform on the right by...
        // | -1 -1  0 |
        // |  0  1  0 |
        // |  0  0 -1 |
        s->a12 -= s->a11; s->a11 = -s->a11; s->a13 = -s->a13;
        s->a22 -= s->a21; s->a21 = -s->a21; s->a23 = -s->a23;
        s->a32 -= s->a31; s->a31 = -s->a31; s->a33 = -s->a33;
        
        // apply the isometry.
        f = g - f;
    }
    
    if (a == g && h > f + f)
    {
        // transform on the right by...
        // | -1  0 -1 |
        // |  0 -1  0 |
        // |  0  0  1 |
        s->a13 -= s->a11; s->a11 = -s->a11; s->a12 = -s->a12;
        s->a23 -= s->a21; s->a21 = -s->a21; s->a22 = -s->a22;
        s->a33 -= s->a31; s->a31 = -s->a31; s->a32 = -s->a32;
        
        // apply the isometry.
        f = h - f;
    }
    
    if (b == f && h > g + g)
    {
        // transform on the right by...
        // | -1  0  0 |
        // |  0 -1 -1 |
        // |  0  0  1 |
        s->a13 -= s->a12; s->a12 = -s->a12; s->a11 = -s->a11;
        s->a23 -= s->a22; s->a22 = -s->a22; s->a21 = -s->a21;
        s->a33 -= s->a32; s->a32 = -s->a32; s->a31 = -s->a31;
        
        // apply the isometry.
        g = h - g;
    }
    
    if (a == b && abs(f) > abs(g))
    {
        // transform on the right by...
        // |  0 -1  0 |
        // | -1  0  0 |
        // |  0  0 -1 |
        temp = -s->a11; s->a11 = -s->a12; s->a12 = temp; s->a13 = -s->a13;
        temp = -s->a21; s->a21 = -s->a22; s->a22 = temp; s->a23 = -s->a23;
        temp = -s->a31; s->a31 = -s->a32; s->a32 = temp; s->a33 = -s->a33;
        
        // apply the isometry.
        a.swap(b);
        g.swap(f);
    }
    
    if (b == c && abs(g) > abs(h))
    {
        // transform on the right by...
        // | -1  0  0 |
        // |  0  0 -1 |
        // |  0 -1  0 |
        temp = s->a12; s->a12 = -s->a13; s->a13 = -temp; s->a11 = -s->a11;
        temp = s->a22; s->a22 = -s->a23; s->a23 = -temp; s->a21 = -s->a21;
        temp = s->a32; s->a32 = -s->a33; s->a33 = -temp; s->a31 = -s->a31;
        
        // apply the isometry.
        g.swap(h);
    }
    
    if (a == b && abs(f) > abs(g))
    {
        // transform on the right by...
        // |  0 -1  0 |
        // | -1  0  0 |
        // |  0  0 -1 |
        temp = s->a11; s->a11 = -s->a12; s->a12 = -temp; s->a13 = -s->a13;
        temp = s->a21; s->a21 = -s->a22; s->a22 = -temp; s->a23 = -s->a23;
        temp = s->a31; s->a31 = -s->a32; s->a32 = -temp; s->a33 = -s->a33;
        
        // apply the isometry.
        g.swap(f);
printf("t\n");
    }

    assert( s->is_isometry(q, a, b, c, f, g, h) );

    auto qq = std::shared_ptr<QuadFormZZ>(
        new QuadForm(q.disc_, a, b, c, f, g, h)
    );

    qq->s_ = saveIsometry ? s : std::make_shared<IsometryQQ>(true);
    qq->reduced_ = true;

    return qq;
}

template<>
std::vector<mpz_class> QuadFormZZ::isotropic_vector(std::shared_ptr<PrimeZZ> pR) const
{
    // The prime integer.
    mpz_class p = abs(pR->principal_generator());

    std::vector<mpz_class> vec(3, 0);

    if (Math::gcd(this->disc_, p) != 1)
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
        vec[0] = (-h * Math::modinv(2*a, p)) % p;
        vec[1] = 1;
        Math::fix_vector(vec, p);
        return vec;
    }

    int64_t ll = Math::kronecker(v, p);


    if (ll == 1)
    {
        mpz_class sqr = Math::square_root(v, p);

        mpz_class temp = Math::modinv(sqr, p);
        vec[0] = (1 - h*temp) % p;
        vec[1] = (2*a*temp) % p;

        if (vec[0] != 0)
        {
            mpz_class temp = Math::modinv(vec[0], p);
            vec[0] = 1;
            vec[1] = (vec[1] * temp) % p;
        }
        else
        {
            vec[1] = 1;
        }

        Math::fix_vector(vec, p);
        return vec;
    }

    mpz_class A = (2*a) % p;
    mpz_class B = ((-v) * Math::modinv(2*a, p)) % p;
    mpz_class C = ((2 * this->disc_) * Math::modinv(-v, p)) % p;

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
                Math::modinv(C, p)) % p;
    }
    while (num == 0 || Math::kronecker(num, p) == -1);

    // Compute its square root modulo p.
    mpz_class sqr = Math::square_root(num, p);

    // Avoid recomputing this value.
    mpz_class temp = Math::modinv(-v, p);

    // Build the isotropic vector.
    vec[0] = (((x - (h * y * Math::modinv(2*a, p))) % p) +
              sqr * ((f * h - 2 * b * g) % p) * temp) % p;
    vec[1] = (y + sqr * (g * h - 2 * a *f) * temp) % p;
    vec[2] = sqr;

    // Normalize the vector.
    if (vec[0] != 0)
    {
        mpz_class inv = Math::modinv(vec[0], p);
        vec[0] = 1;
        vec[1] = (vec[1] * inv) % p;
        vec[2] = (vec[2] * inv) % p;
    }
    else if (vec[1] != 0)
    {
        mpz_class inv = Math::modinv(vec[1], p);
        vec[1] = 1;
        vec[2] = (vec[2] * inv) % p;
    }
    else
    {
        vec[2] = 1;
    }

    Math::fix_vector(vec, p);
    return vec;
}
