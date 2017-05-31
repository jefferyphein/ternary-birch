#include <memory>
#include <iostream>
#include <iomanip>
#include <cassert>
#include "NeighborIterator.h"
#include "QuadForm.h"
#include "Math.h"

typedef NeighborIterator<mpz_class, mpq_class> NeighborIteratorZZ;
typedef QuadForm<mpz_class, mpq_class> QuadFormZZ;
typedef Prime<mpz_class, mpq_class> PrimeZZ;
typedef ChangeOfBasis<mpz_class, mpq_class> ChangeOfBasisZZ;
typedef Isometry<mpz_class, mpq_class> IsometryQQ;

template<>
NeighborIteratorZZ::NeighborIterator(std::shared_ptr<QuadFormZZ> q, std::shared_ptr<PrimeZZ> pR) :
    q_(q),
    p_(pR),
    numNeighbors_(mpz_class(pR->norm() + 1)),
    pos_(mpz_class(0)),
    isotropicVector_(q->isotropic_vector(pR)),
    s_(std::make_shared<ChangeOfBasisZZ>(pR))
{
    // The prime integer.
    mpz_class p = pR->principal_generator();

    // Normalize the initial isotropic vector.
    Math::normalize_vector(this->isotropicVector_, p);

    /* Our goal at this point is to build a p-standard basis which can then
     * be utilized to help us construct all of the p-neighbors of the attached
     * quadratic form. */

    // These coefficients will help us track p-standard basis transformation.
    mpz_class a, b, c, f, g, h;

    if (this->isotropicVector_[0] == 1)
    {
        mpz_class u = this->isotropicVector_[1];
        mpz_class v = this->isotropicVector_[2];

        this->s_->a11 = 1;
        this->s_->a21 = u;
        this->s_->a31 = v;
        this->s_->a22 = 1;
        this->s_->a33 = 1;

        a = (this->q_->a() +
             this->q_->b() * u * u +
             this->q_->c() * v * v +
             this->q_->f() * u * v +
             this->q_->g() * v +
             this->q_->h() * u) % p;
        b = this->q_->b() % p;
        c = this->q_->c() % p;
        f = this->q_->f() % p;
        g = (2 * this->q_->c() * v +
                 this->q_->f() * u +
                 this->q_->g()) % p;
        h = (2 * this->q_->b() * u +
                 this->q_->f() * v +
                 this->q_->h()) % p;
    }
    else if (this->isotropicVector_[1] == 1)
    {
        mpz_class u = this->isotropicVector_[2];

        this->s_->a21 = 1;
        this->s_->a31 = u;
        this->s_->a12 = 1;
        this->s_->a33 = 1;

        a = (this->q_->b() +
             this->q_->c() * u * u +
             this->q_->f() * u) % p;
        b = this->q_->a() % p;
        c = this->q_->c() % p;
        f = this->q_->g() % p;
        g = (2 * this->q_->c() * u +
                 this->q_->f()) % p;
        h = (this->q_->g() * u +
             this->q_->h()) % p;
    }
    else
    {
        this->s_->a31 = 1;
        this->s_->a12 = 1;
        this->s_->a23 = 1;

        a = this->q_->c() % p;
        b = this->q_->a() % p;
        c = this->q_->b() % p;
        f = this->q_->h() % p;
        g = this->q_->f() % p;
        h = this->q_->g() % p;
    }

    // If (g % p) == 0, we need to swap the y and z indeterminants, this
    // forces the xz coordinate to be nonzero modulo p upon applying the
    // change of basis over QQ.
    if (g % p == 0)
    {
        this->s_->a13.swap(this->s_->a12);
        this->s_->a23.swap(this->s_->a22);
        this->s_->a33.swap(this->s_->a32);

        b.swap(c);
        g.swap(h);
    }

    assert( (g % p) != 0 );

    // Scalar to make (h % p) == 0.
    mpz_class scalar = (-h * Math::modinv(g, p)) % p;

    if (p != 2)
    {
        // Make the absolute value of this scalar as small as possible.
        if (scalar < 0 && -scalar > (p-1)/2) { scalar += p; }
        else if (scalar > 0 && scalar > (p-1)/2) { scalar -= p; }
    }

    b = (b + scalar * (f + c * scalar)) % p;
    f = (f + 2 * c * scalar) % p;
    h = (h + scalar * g) % p;

    assert( h == 0 );

    // Multiply the isometry on the right by...
    // |  1  0  0 |
    // |  0  1  0 |
    // |  0  s  1 |
    // where s is the scalar defined above.
    this->s_->a12 = (this->s_->a12 + this->s_->a13 * scalar) % p;
    this->s_->a22 = (this->s_->a22 + this->s_->a23 * scalar) % p;
    this->s_->a32 = (this->s_->a32 + this->s_->a33 * scalar) % p;

    mpz_class ginv = Math::modinv(g, p);
    this->s_->a12 = (this->s_->a12 - this->s_->a11 * f * ginv) % p;
    this->s_->a22 = (this->s_->a22 - this->s_->a21 * f * ginv) % p;
    this->s_->a32 = (this->s_->a32 - this->s_->a31 * f * ginv) % p;
    this->s_->a13 = (ginv * (this->s_->a13 - this->s_->a11 * c * ginv)) % p;
    this->s_->a23 = (ginv * (this->s_->a23 - this->s_->a21 * c * ginv)) % p;
    this->s_->a33 = (ginv * (this->s_->a33 - this->s_->a31 * c * ginv)) % p;

    // Set the anisotropic value for this p-standard basis.
    this->aniso_ = b % p;
    if (p != 2)
    {
        if (this->aniso_ < 0 && -this->aniso_ < (p-1)/2) { this->aniso_ += p; }
        else if (this->aniso_ > 0 && this->aniso_ > (p-1)/2) { this->aniso_ -= p; }
    }
}

template<>
std::shared_ptr<QuadFormZZ> NeighborIteratorZZ::build_neighbor(std::vector<mpz_class>& vec)
{
    // The prime integer.
    mpz_class p = this->p_->principal_generator();

    // Normalize the isotropic vector modulo p, then adjust its coefficients
    // so that their absolute value is as small as possible.
    Math::normalize_vector(vec, p);
    Math::fix_vector(vec, p);

    // The coefficients of the p-neighbor.
    mpz_class a, b, c, f, g, h;

    // The isometry.
    std::shared_ptr<IsometryQQ> s = std::make_shared<IsometryQQ>(true);

    if (vec[0] == 1)
    {
        mpz_class u = vec[1];
        mpz_class v = vec[2];

        // Initialize isometry to...
        // | 1  0  0 |
        // | u  1  0 |
        // | v  0  1 |
        s->a11 = 1;
        s->a21 = u;
        s->a31 = v;

        // Update coefficients.
        a = this->q_->a() +
            this->q_->b() * u * u +
            this->q_->c() * v * v +
            this->q_->f() * u * v +
            this->q_->g() * v +
            this->q_->h() * u;
        b = this->q_->b();
        c = this->q_->c();
        f = this->q_->f();
        g = 2 * this->q_->c() * v +
            this->q_->f() * u +
            this->q_->g();
        h = 2 * this->q_->b() * u +
            this->q_->f() * v +
            this->q_->h();
    }
    else if (vec[1] == 1)
    {
        mpz_class u = vec[2];

        // Initialize isometry to be...
        // | 0  1  0 |
        // | 1  0  0 |
        // | u  0  1 |
        s->a11 = 0;
        s->a22 = 0;
        s->a12 = 1;
        s->a21 = 1;
        s->a31 = u;

        // Update coefficients.
        a = this->q_->b() +
            this->q_->c() * u * u +
            this->q_->f() * u;
        b = this->q_->a();
        c = this->q_->c();
        f = this->q_->g();
        g = 2 * this->q_->c() * u +
            this->q_->f();
        h = this->q_->g() * u +
            this->q_->h();
    }
    else
    {
        // Initialize isometry to...
        // |  0  1  0 |
        // |  0  0  1 |
        // |  1  0  0 |
        s->a11 = 0;
        s->a22 = 0;
        s->a33 = 0;
        s->a31 = 1;
        s->a12 = 1;
        s->a23 = 1;

        // Update coefficients.
        a = this->q_->c();
        b = this->q_->a();
        c = this->q_->b();
        f = this->q_->h();
        g = this->q_->f();
        h = this->q_->g();
    }

    if (g % p == 0)
    {
        // Multiply isometry on the right by...
        // |  1  0  0 |
        // |  0  0  1 |
        // |  0  1  0 |
        s->a13.swap(s->a12);
        s->a23.swap(s->a22);
        s->a33.swap(s->a32);

        b.swap(c);
        g.swap(h);
    }

    // Scalar to make (h % p) == 0.
    mpz_class scalar = (-h * Math::modinv(g, p)) % p;

    if (p != 2)
    {
        if (scalar < 0 && -scalar > (p-1)/2) { scalar += p; }
        else if (scalar > 0 && scalar > (p-1)/2) { scalar -= p; }
    }

    // Multiply on the right by...
    // |  1  0  0 |
    // |  0  1  0 |
    // |  0  s  1 |
    // where s is the scalar defined above.
    s->a12 += (s->a13 * scalar);
    s->a22 += (s->a23 * scalar);
    s->a32 += (s->a33 * scalar);

    // Update coefficients.
    b += scalar * (f + c * scalar);
    f += 2 * c * scalar;
    h += scalar * g;

    assert( h % p == 0 );

    // Avoid recomputing p*p repeatedly.
    mpz_class pp = p*p;

    // Scalar to make (a % pp) == 0.
    scalar = (-a * Math::modinv(g, pp)) % pp;

    if (p != 2)
    {
        if (scalar < 0 && -scalar > (pp-1)/2) { scalar += pp; }
        else if (scalar > 0 && scalar > (pp-1)/2) { scalar -= pp; }
    }

    // Multiply on the right by...
    // |  1  0  0 |
    // |  0  1  0 |
    // |  s  0  1 |
    s->a11 += (s->a13 * scalar);
    s->a21 += (s->a23 * scalar);
    s->a31 += (s->a33 * scalar);

    // Update coefficients.
    a += (scalar * (g + c * scalar));
    g += (2 * c * scalar);
    h += (scalar * f);

    assert( a % pp == 0 );

    // Multiply on the right by...
    // | 1/p  0   0 |
    // |  0   1   0 |
    // |  0   0   p |
    s->a11 /= p;
    s->a21 /= p;
    s->a31 /= p;
    s->a13 *= p;
    s->a23 *= p;
    s->a33 *= p;

    // Update coefficients.
    a /= pp;
    c *= pp;
    f *= p;
    h /= p;

    assert( s->is_isometry(*this->q_, a, b, c, f, g, h) );

    //s->multiply_on_left_by(this->q_->isometry());

    auto qq = std::make_shared<QuadFormZZ>(this->q_->discriminant(),
                                           a, b, c, f, g, h);
    qq->isometry(s);

    return qq;
}

template<>
std::shared_ptr<QuadFormZZ> NeighborIteratorZZ::next_neighbor(void)
{
    // Return nullptr if no more neighbors to compute.
    if (this->pos_ >= this->numNeighbors_) { return nullptr; }

    mpz_class p = this->p_->principal_generator();

    // Make a copy of the initial isotropic vector.
    std::vector<mpz_class> vec(this->isotropicVector_);

    // Create a new isotropic vector based on our current position.
    if (this->pos_ < p)
    {
        vec[0] = (vec[0] + this->pos_ *
            (this->s_->a12 - this->s_->a13 * this->aniso_ * this->pos_)) % p;
        vec[1] = (vec[1] + this->pos_ *
            (this->s_->a22 - this->s_->a23 * this->aniso_ * this->pos_)) % p;
        vec[2] = (vec[2] + this->pos_ *
            (this->s_->a32 - this->s_->a33 * this->aniso_ * this->pos_)) % p;
    }
    else
    {
        vec[0] = this->s_->a13;
        vec[1] = this->s_->a23;
        vec[2] = this->s_->a33;
    }

    assert( this->q_->evaluate(vec) % p == 0 );

    auto qq = this->build_neighbor(vec);

//    std::cout << "Q_{" << p << "}(" << vec[0] << ","
//              << vec[1] << "," 
//              << vec[2] << ") = " << (this->q_->evaluate(vec) % p) << std::endl;

    ++this->pos_;

    return qq;
}
