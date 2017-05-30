#include <iostream>
#include <gmpxx.h>
#include "QuadForm.h"

template<typename T>
QuadForm<T>::QuadForm(const T& a, const T& b, const T& c,
                      const T& f, const T& g, const T& h)
{
    T disc = a * (4 * b * c - f * f) - b * g * g + h * (f * g - c * h);
    *this = QuadForm<T>(disc, a, b, c, f, g, h);
}

template<typename T>
QuadForm<T>::QuadForm(const T& disc,
                      const T& a, const T& b, const T& c,
                      const T& f, const T& g, const T& h)
{
    this->disc = disc;
    this->a = a; this->b = b; this->c = c;
    this->f = f; this->g = g; this->h = h;
}

template<typename T>
const T& QuadForm<T>::discriminant(void) const
{
    return this->disc;
}

template<typename T>
void QuadForm<T>::print(std::ostream& os) const
{
    os << this->a << " " << this->b << " " << this->c << " "
       << this->f << " " << this->g << " " << this->h;
}

/* Tell the compiler to compile  */
template class QuadForm<mpz_class>;
