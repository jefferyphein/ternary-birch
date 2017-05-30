#ifndef __ISOMETRY_H_
#define __ISOMETRY_H_

#include <iostream>
#include "QuadForm.h"

/* Forward declaration for friend class. */
template<typename R, typename F>
class QuadForm;

template<typename R, typename F>
class NeighborIterator;

template<typename R, typename F>
class Isometry
{
/* Allow these classes to manipulate the private members of this class. */
friend class QuadForm<R,F>;
friend class NeighborIterator<R,F>;

public:
    Isometry(bool id = false);

    void print(std::ostream& os) const;

    /* Verify:
     *  Transpose(s) * Gram(q) * s == Gram(QuadForm(a,b,c,f,g,h))
     */
    bool is_isometry(const QuadForm<R,F>& q,
                     const R& a, const R& b, const R& c,
                     const R& f, const R& g, const R& h) const;
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
void Isometry<R,F>::print(std::ostream& os) const
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

#endif // __ISOMETRY_H_
