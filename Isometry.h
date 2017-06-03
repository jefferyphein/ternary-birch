#ifndef __ISOMETRY_H_
#define __ISOMETRY_H_

#include <iostream>
#include "QuadForm.h"
#include "Prime.h"

/* Forward declaration for friend class. */
template<typename R, typename F>
class QuadForm;

template<typename R, typename F>
class NeighborIterator;

template<typename R, typename F>
class ChangeOfBasis
{
friend class NeighborIterator<R,F>;
public:
    ChangeOfBasis(std::shared_ptr<Prime<R,F>> pR) : pR_(pR) {}

    void print(std::ostream& os) const;

private:
    std::shared_ptr<Prime<R,F>> pR_;
    R a11, a12, a13, a21, a22, a23, a31, a32, a33;
};

template<typename R, typename F>
class Isometry
{
/* Allow these classes to manipulate the private members of this class. */
friend class QuadForm<R,F>;
friend class NeighborIterator<R,F>;

public:
    Isometry(bool id = false);
    Isometry(F a11, F a12, F a13,
             F a21, F a22, F a23,
             F a31, F a32, F a33);

    void print(std::ostream& os) const;

    /* Verify:
     *  Transpose(s) * Gram(q) * s == Gram(QuadForm(a,b,c,f,g,h))
     */
    bool is_isometry(const QuadForm<R,F>& q,
                     const R& a, const R& b, const R& c,
                     const R& f, const R& g, const R& h) const;
    bool is_isometry(const QuadForm<R,F>& from, const QuadForm<R,F>& to) const;

    bool is_automorphism(const QuadForm<R,F>& q) const;

    void multiply_on_right_by(std::shared_ptr<Isometry<R,F>> s);
    void multiply_on_left_by(std::shared_ptr<Isometry<R,F>> s);

    F spinor_norm(const QuadForm<R,F>& q) const;
    F trace(void) const;

    std::shared_ptr<Isometry<R,F>> inverse(void) const;

private:
    F a11, a12, a13, a21, a22, a23, a31, a32, a33;
};

template<typename R, typename F>
Isometry<R,F>::Isometry(bool id)
{
    if (id)
    {
        this->a11 = F(1);
        this->a22 = F(1);
        this->a33 = F(1);
    }
}

template<typename R, typename F>
Isometry<R,F>::Isometry(F a11, F a12, F a13,
                        F a21, F a22, F a23,
                        F a31, F a32, F a33) :
    a11(a11), a12(a12), a13(a13),
    a21(a21), a22(a22), a23(a23),
    a31(a31), a32(a32), a33(a33) {}

template<typename R, typename F>
std::ostream& operator<<(std::ostream& os, const Isometry<R,F>& s)
{
    s.print(os);
    return os;
}

template<typename R, typename F>
std::ostream& operator<<(std::ostream& os, std::shared_ptr<Isometry<R,F>> sPtr)
{
    sPtr->print(os);
    return os;
}

template<typename R, typename F>
std::ostream& operator<<(std::ostream& os, const ChangeOfBasis<R,F>& s)
{
    s.print(os);
    return os;
}

template<typename R, typename F>
std::ostream& operator<<(std::ostream& os, std::shared_ptr<ChangeOfBasis<R,F>> sPtr)
{
    sPtr->print(os);
    return os;
}

template<typename R, typename F>
void Isometry<R,F>::print(std::ostream& os) const
{
    os << "[ ";
    os << "[ " << this->a11 << ", " << this->a12 << ", " << this->a13 << " ], ";
    os << "[ " << this->a21 << ", " << this->a22 << ", " << this->a23 << " ], ";
    os << "[ " << this->a31 << ", " << this->a32 << ", " << this->a33 << " ]";
    os << " ]";
}

template<typename R, typename F>
void ChangeOfBasis<R,F>::print(std::ostream& os) const
{
    os << "[ ";
    os << "[ " << this->a11 << ", " << this->a12 << ", " << this->a13 << " ], ";
    os << "[ " << this->a21 << ", " << this->a22 << ", " << this->a23 << " ], ";
    os << "[ " << this->a31 << ", " << this->a32 << ", " << this->a33 << " ]";
    os << " ]";
}

template<typename R, typename F>
bool Isometry<R,F>::is_isometry(const QuadForm<R,F>& q,
                                const R& a, const R& b, const R& c,
                                const R& f, const R& g, const R& h) const
{
    if ((q.a()*this->a11*this->a11 +
         q.b()*this->a21*this->a21 +
         q.c()*this->a31*this->a31 +
         q.f()*this->a21*this->a31 +
         q.g()*this->a11*this->a31 +
         q.h()*this->a11*this->a21) != a) { return false; }

    if ((q.a()*this->a12*this->a12 +
         q.b()*this->a22*this->a22 +
         q.c()*this->a32*this->a32 +
         q.f()*this->a22*this->a32 +
         q.g()*this->a12*this->a32 +
         q.h()*this->a12*this->a22) != b) { return false; }

    if ((q.a()*this->a13*this->a13 +
         q.b()*this->a23*this->a23 +
         q.c()*this->a33*this->a33 +
         q.f()*this->a23*this->a33 +
         q.g()*this->a13*this->a33 +
         q.h()*this->a13*this->a23) != c) { return false; }

    if ((F(2)*q.a()*this->a12*this->a13 +
         F(2)*q.b()*this->a22*this->a23 +
         F(2)*q.c()*this->a32*this->a33 +
         q.f()*this->a22*this->a33 +
         q.f()*this->a23*this->a32 +
         q.g()*this->a12*this->a33 +
         q.g()*this->a13*this->a32 +
         q.h()*this->a12*this->a23 +
         q.h()*this->a13*this->a22) != f) { return false; }

    if ((F(2)*q.a()*this->a11*this->a13 +
         F(2)*q.b()*this->a21*this->a23 +
         F(2)*q.c()*this->a31*this->a33 +
         q.f()*this->a21*this->a33 +
         q.f()*this->a23*this->a31 +
         q.g()*this->a11*this->a33 +
         q.g()*this->a13*this->a31 +
         q.h()*this->a11*this->a23 +
         q.h()*this->a13*this->a21) != g) { return false; }

    if ((F(2)*q.a()*this->a11*this->a12 +
         F(2)*q.b()*this->a21*this->a22 +
         F(2)*q.c()*this->a31*this->a32 +
         q.f()*this->a21*this->a32 +
         q.f()*this->a22*this->a31 +
         q.g()*this->a11*this->a32 +
         q.g()*this->a12*this->a31 +
         q.h()*this->a11*this->a22 +
         q.h()*this->a12*this->a21) != h) { return false; }

    return true;
}

template<typename R, typename F>
bool Isometry<R,F>::is_automorphism(const QuadForm<R,F>& q) const
{
    return this->is_isometry(q, q);
}

template<typename R, typename F>
bool Isometry<R,F>::is_isometry(const QuadForm<R,F>& from, const QuadForm<R,F>& to) const
{
    return this->is_isometry(from, to.a(), to.b(), to.c(), to.f(), to.g(), to.h());
}

template<typename R, typename F>
void Isometry<R,F>::multiply_on_right_by(std::shared_ptr<Isometry<R,F>> s)
{
    F a11 = this->a11;
    F a12 = this->a12;
    F a13 = this->a13;
    F a21 = this->a21;
    F a22 = this->a22;
    F a23 = this->a23;
    F a31 = this->a31;
    F a32 = this->a32;
    F a33 = this->a33;
   
    this->a11 = a11*s->a11 + a12*s->a21 + a13*s->a31;
    this->a12 = a11*s->a12 + a12*s->a22 + a13*s->a32;
    this->a13 = a11*s->a13 + a12*s->a23 + a13*s->a33;
    this->a21 = a21*s->a11 + a22*s->a21 + a23*s->a31;
    this->a22 = a21*s->a12 + a22*s->a22 + a23*s->a32;
    this->a23 = a21*s->a13 + a22*s->a23 + a23*s->a33;
    this->a31 = a31*s->a11 + a32*s->a21 + a33*s->a31;
    this->a32 = a31*s->a12 + a32*s->a22 + a33*s->a32;
    this->a33 = a31*s->a13 + a32*s->a23 + a33*s->a33;
}

template<typename R, typename F>
void Isometry<R,F>::multiply_on_left_by(std::shared_ptr<Isometry<R,F>> s)
{
    F a11 = this->a11;
    F a12 = this->a12;
    F a13 = this->a13;
    F a21 = this->a21;
    F a22 = this->a22;
    F a23 = this->a23;
    F a31 = this->a31;
    F a32 = this->a32;
    F a33 = this->a33;
    
    this->a11 = a11*s->a11 + a21*s->a12 + a31*s->a13;
    this->a12 = a12*s->a11 + a22*s->a12 + a32*s->a13;
    this->a13 = a13*s->a11 + a23*s->a12 + a33*s->a13;
    this->a21 = a11*s->a21 + a21*s->a22 + a31*s->a23;
    this->a22 = a12*s->a21 + a22*s->a22 + a32*s->a23;
    this->a23 = a13*s->a21 + a23*s->a22 + a33*s->a23;
    this->a31 = a11*s->a31 + a21*s->a32 + a31*s->a33;
    this->a32 = a12*s->a31 + a22*s->a32 + a32*s->a33;
    this->a33 = a13*s->a31 + a23*s->a32 + a33*s->a33;
}

template<typename R, typename F>
std::shared_ptr<Isometry<R,F>> Isometry<R,F>::inverse(void) const
{
    F a11 = this->a11;
    F a12 = this->a12;
    F a13 = this->a13;
    F a21 = this->a21;
    F a22 = this->a22;
    F a23 = this->a23;
    F a31 = this->a31;
    F a32 = this->a32;
    F a33 = this->a33;

    Isometry<R,F> inv = std::make_shared<Isometry<R,F>>(false);
    inv->a11 = a22 * a33 - a23 * a32;
    inv->a12 = a13 * a32 - a12 * a33;
    inv->a13 = a12 * a23 - a13 * a22;
    inv->a21 = a23 * a31 - a21 * a33;
    inv->a22 = a11 * a33 - a13 * a31;
    inv->a23 = a13 * a21 - a11 * a23;
    inv->a31 = a21 * a32 - a22 * a31;
    inv->a32 = a12 * a31 - a11 * a32;
    inv->a33 = a11 * a22 - a12 * a21;

    return inv;
}

template<typename R, typename F>
F Isometry<R,F>::trace(void) const
{
    return this->a11 + this->a22 + this->a33;
}

#endif // __ISOMETRY_H_
