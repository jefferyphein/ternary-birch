#ifndef __QUAD_FORM_H_
#define __QUAD_FORM_H_

#include <memory>
#include <iomanip>
#include <gmpxx.h>
#include "Isometry.h"
#include "Prime.h"

template<typename R, typename F>
class QuadForm
{
public:
    QuadForm(const QuadForm<R,F>&) = default;

    QuadForm(const R& a, const R& b, const R& c,
             const R& f, const R& g, const R& h);
    QuadForm(const R& disc,
             const R& a, const R& b, const R& c,
             const R& f, const R& g, const R& h);

    const R& a(void) const;
    const R& b(void) const;
    const R& c(void) const;
    const R& f(void) const;
    const R& g(void) const;
    const R& h(void) const;
    const R& discriminant() const;

    std::shared_ptr<Isometry<R,F>> isometry() const;
    void isometry(std::shared_ptr<Isometry<R,F>> s);

    R evaluate(const R& x, const R& y, const R& z) const;
    R evaluate(const std::vector<R>& vec) const;


    void print(std::ostream& os) const;

    std::vector<R> isotropic_vector(std::shared_ptr<Prime<R,F>> pR) const;

    static std::shared_ptr<QuadForm<R,F>> reduce(const QuadForm<R,F>& q,
                                                 bool saveIsometry = true);
private:
    /* The coefficients of the quadratic form. The quadratic form has shape:
     *      Q(x,y,z) = ax^2 + by^2 + cz^2 + fyz + gxz + hxy
     */
    R a_, b_, c_, f_, g_, h_;

    /* The discriminant of the quadratic form. */
    R disc_ = -1;

    /* An isometry that transforms this form to the "mother" form defined by
     * the genus to which it belongs. This is unused if this form does not get
     * assigned to a genus as a representative of an isometry class. */
    std::shared_ptr<Isometry<R,F>> s_;

    /* Flag indicating whether this is a reduced form or not. */
    bool reduced_ = false;
};

/* Generic template implementation follows.
 *
 * This code belongs in the header file so that users can make use of this
 * code if they decide to define their own template types, e.g. if taking
 * R = ZZ[\sqrt(2)]. This code will then be compiled for that data type as
 * needed.
 */

template<typename R, typename F>
QuadForm<R,F>::QuadForm(const R& a, const R& b, const R& c,
                        const R& f, const R& g, const R& h)
{
    R disc = a * (4 * b * c - f * f) - b * g * g + h * (f * g - c * h);
    *this = QuadForm<R,F>(disc, a, b, c, f, g, h);
}

template<typename R, typename F>
QuadForm<R,F>::QuadForm(const R& disc,
                        const R& a, const R& b, const R& c,
                        const R& f, const R& g, const R& h)
{
    this->disc_ = disc;
    this->a_ = a; this->b_ = b; this->c_ = c;
    this->f_ = f; this->g_ = g; this->h_ = h;
}

template<typename R, typename F>
const R& QuadForm<R,F>::a(void) const
{
    return this->a_;
}

template<typename R, typename F>
const R& QuadForm<R,F>::b(void) const
{
    return this->b_;
}

template<typename R, typename F>
const R& QuadForm<R,F>::c(void) const
{
    return this->c_;
}

template<typename R, typename F>
const R& QuadForm<R,F>::f(void) const
{
    return this->f_;
}

template<typename R, typename F>
const R& QuadForm<R,F>::g(void) const
{
    return this->g_;
}

template<typename R, typename F>
const R& QuadForm<R,F>::h(void) const
{
    return this->h_;
}

template<typename R, typename F>
const R& QuadForm<R,F>::discriminant(void) const
{
    return this->disc_;
}

template<typename R, typename F>
std::shared_ptr<Isometry<R,F>> QuadForm<R,F>::isometry(void) const
{
    return this->s_;
}

template<typename R, typename F>
void QuadForm<R,F>::isometry(std::shared_ptr<Isometry<R,F>> s)
{
    this->s_ = s;
}

template<typename R, typename F>
R QuadForm<R,F>::evaluate(const R& x, const R& y, const R& z) const
{
    return x * (this->a_ * x + this->h_ * y) +
           y * (this->b_ * y + this->f_ * z) +
           z * (this->c_ * z + this->g_ * x);
}

template<typename R, typename F>
R QuadForm<R,F>::evaluate(const std::vector<R>& vec) const
{
    return this->evaluate(vec[0], vec[1], vec[2]);
}

template<typename R, typename F>
std::ostream& operator<<(std::ostream& os, const QuadForm<R,F>& q)
{
    q.print(os);
    return os;
}

template<typename R, typename F>
std::ostream& operator<<(std::ostream& os, std::shared_ptr<QuadForm<R,F>> qPtr)
{
    qPtr->print(os);
    return os;
}

template<typename R, typename F>
void QuadForm<R,F>::print(std::ostream& os) const
{
    os << this->a_ << " " << this->b_ << " " << this->c_ << " "
       << this->f_ << " " << this->g_ << " " << this->h_;
}

template<typename R, typename F>
bool operator<(const QuadForm<R,F>& q1, const QuadForm<R,F>& q2)
{
    if (q1.a() < q2.a()) { return true; }
    else if (q1.a() > q2.a()) { return false; }
    else if (q1.b() < q2.b()) { return true; }
    else if (q1.b() > q2.b()) { return false; }
    else if (q1.c() < q2.c()) { return true; }
    else if (q1.c() > q2.c()) { return false; }
    else if (q1.f() < q2.f()) { return true; }
    else if (q1.f() > q2.f()) { return false; }
    else if (q1.g() < q2.g()) { return true; }
    else if (q1.g() > q2.g()) { return false; }
    else if (q1.h() < q2.h()) { return true; }
    else if (q1.h() > q2.g()) { return false; }
    else { return false; }
}

#endif // __QUAD_FORM_H_
