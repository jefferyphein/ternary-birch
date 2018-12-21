#ifndef __ISOMETRY_H_
#define __ISOMETRY_H_

#include "birch.h"

template<typename R>
class Isometry
{
public:
    Isometry() : a11(1), a12(0), a13(0), a21(0), a22(1), a23(0), a31(0), a32(0), a33(1) {}

    Isometry(const R& a11, const R& a12, const R& a13,
             const R& a21, const R& a22, const R& a23,
             const R& a31, const R& a32, const R& a33) :
        a11(a11), a12(a12), a13(a13), a21(a21), a22(a22), a23(a23), a31(a31), a32(a32), a33(a33) {}

    void set_values(const R& a11, const R& a12, const R& a13,
                    const R& a21, const R& a22, const R& a23,
                    const R& a31, const R& a32, const R& a33)
    {
        this->a11 = a11; this->a12 = a12; this->a13 = a13;
        this->a21 = a21; this->a22 = a22; this->a23 = a23;
        this->a31 = a31; this->a32 = a32; this->a33 = a33;
    }

    void set_identity(void)
    {
        this->a11 = 1; this->a12 = 0; this->a13 = 0;
        this->a21 = 0; this->a22 = 1; this->a23 = 0;
        this->a31 = 0; this->a32 = 0; this->a33 = 1;
    }

    Isometry<R> inverse(const R& p) const
    {
        Isometry<R> temp;
        temp.a11 = (this->a22 * this->a33 - this->a23 * this->a32) / p;
        temp.a12 = (this->a13 * this->a32 - this->a12 * this->a33) / p;
        temp.a13 = (this->a12 * this->a23 - this->a13 * this->a22) / p;
        temp.a21 = (this->a23 * this->a31 - this->a21 * this->a33) / p;
        temp.a22 = (this->a11 * this->a33 - this->a13 * this->a31) / p;
        temp.a23 = (this->a13 * this->a21 - this->a11 * this->a23) / p;
        temp.a31 = (this->a21 * this->a32 - this->a22 * this->a31) / p;
        temp.a32 = (this->a12 * this->a31 - this->a11 * this->a32) / p;
        temp.a33 = (this->a11 * this->a22 - this->a12 * this->a21) / p;
        return temp;
    }

    Isometry<R> operator*(const Isometry<R>& s) const
    {
        Isometry<R> temp;
        temp.a11 = this->a11*s.a11 + this->a12*s.a21 + this->a13*s.a31;
        temp.a12 = this->a11*s.a12 + this->a12*s.a22 + this->a13*s.a32;
        temp.a13 = this->a11*s.a13 + this->a12*s.a23 + this->a13*s.a33;
        temp.a21 = this->a21*s.a11 + this->a22*s.a21 + this->a23*s.a31;
        temp.a22 = this->a21*s.a12 + this->a22*s.a22 + this->a23*s.a32;
        temp.a23 = this->a21*s.a13 + this->a22*s.a23 + this->a23*s.a33;
        temp.a31 = this->a31*s.a11 + this->a32*s.a21 + this->a33*s.a31;
        temp.a32 = this->a31*s.a12 + this->a32*s.a22 + this->a33*s.a32;
        temp.a33 = this->a31*s.a13 + this->a32*s.a23 + this->a33*s.a33;
        return temp;
    }

    Vector3<R> operator*(const Vector3<R>& vec) const
    {
        Vector3<R> temp;
        temp.x = this->a11 * vec.x + this->a12 * vec.y + this->a13 * vec.z;
        temp.y = this->a21 * vec.x + this->a22 * vec.y + this->a23 * vec.z;
        temp.z = this->a31 * vec.x + this->a32 * vec.y + this->a33 * vec.z;
        return temp;
    }

    bool is_isometry(const QuadForm<R>& from, const QuadForm<R>& to, R scalar)
    {
         if ((from.a()*this->a11*this->a11 +
             from.b()*this->a21*this->a21 +
             from.c()*this->a31*this->a31 +
             from.f()*this->a21*this->a31 +
             from.g()*this->a11*this->a31 +
             from.h()*this->a11*this->a21) != to.a() * scalar) { return false; }
    
        if ((from.a()*this->a12*this->a12 +
             from.b()*this->a22*this->a22 +
             from.c()*this->a32*this->a32 +
             from.f()*this->a22*this->a32 +
             from.g()*this->a12*this->a32 +
             from.h()*this->a12*this->a22) != to.b() * scalar) { return false; }
    
        if ((from.a()*this->a13*this->a13 +
             from.b()*this->a23*this->a23 +
             from.c()*this->a33*this->a33 +
             from.f()*this->a23*this->a33 +
             from.g()*this->a13*this->a33 +
             from.h()*this->a13*this->a23) != to.c() * scalar) { return false; }
    
        if ((2*from.a()*this->a12*this->a13 +
             2*from.b()*this->a22*this->a23 +
             2*from.c()*this->a32*this->a33 +
             from.f()*this->a22*this->a33 +
             from.f()*this->a23*this->a32 +
             from.g()*this->a12*this->a33 +
             from.g()*this->a13*this->a32 +
             from.h()*this->a12*this->a23 +
             from.h()*this->a13*this->a22) != to.f() * scalar) { return false; }
    
        if ((2*from.a()*this->a11*this->a13 +
             2*from.b()*this->a21*this->a23 +
             2*from.c()*this->a31*this->a33 +
             from.f()*this->a21*this->a33 +
             from.f()*this->a23*this->a31 +
             from.g()*this->a11*this->a33 +
             from.g()*this->a13*this->a31 +
             from.h()*this->a11*this->a23 +
             from.h()*this->a13*this->a21) != to.g() * scalar) { return false; }
    
        if ((2*from.a()*this->a11*this->a12 +
             2*from.b()*this->a21*this->a22 +
             2*from.c()*this->a31*this->a32 +
             from.f()*this->a21*this->a32 +
             from.f()*this->a22*this->a31 +
             from.g()*this->a11*this->a32 +
             from.g()*this->a12*this->a31 +
             from.h()*this->a11*this->a22 +
             from.h()*this->a12*this->a21) != to.h() * scalar) { return false; }
   
        return true;
    }

    void A101011001()
    {
        this->a13 += (this->a11 + this->a12);
        this->a23 += (this->a21 + this->a22);
        this->a33 += (this->a31 + this->a32);
    }
    
    void A1t0010001(const R& t)
    {
        this->a12 += t * this->a11;
        this->a22 += t * this->a21;
        this->a32 += t * this->a31;
    }
    
    void A10001t001(const R& t)
    {
        this->a13 += t * this->a12;
        this->a23 += t * this->a22;
        this->a33 += t * this->a32;
    }
    
    void A10t010001(const R& t)
    {
        this->a13 += t * this->a11;
        this->a23 += t * this->a21;
        this->a33 += t * this->a31;
    }
    
    void A0n0n0000n()
    {
        R temp;
        temp = -this->a11; this->a11 = -this->a12; this->a12 = temp; this->a13 = -this->a13;
        temp = -this->a21; this->a21 = -this->a22; this->a22 = temp; this->a23 = -this->a23;
        temp = -this->a31; this->a31 = -this->a32; this->a32 = temp; this->a33 = -this->a33;
    }
    
    void An0000n0n0()
    {
        R temp;
        temp = -this->a12; this->a12 = -this->a13; this->a13 = temp; this->a11 = -this->a11;
        temp = -this->a22; this->a22 = -this->a23; this->a23 = temp; this->a21 = -this->a21;
        temp = -this->a32; this->a32 = -this->a33; this->a33 = temp; this->a31 = -this->a31;
    }
    
    void An00010001()
    {
        this->a11 = -this->a11;
        this->a21 = -this->a21;
        this->a31 = -this->a31;
    }
    
    void A1000n0001()
    {
        this->a12 = -this->a12;
        this->a22 = -this->a22;
        this->a32 = -this->a32;
    }
    
    void A10001000n()
    {
        this->a13 = -this->a13;
        this->a23 = -this->a23;
        this->a33 = -this->a33;
    }
    
    void An010n1001()
    {
        this->a13 += (this->a11 + this->a12); this->a11 = -this->a11; this->a12 = -this->a12;
        this->a23 += (this->a21 + this->a22); this->a21 = -this->a21; this->a22 = -this->a22;
        this->a33 += (this->a31 + this->a32); this->a31 = -this->a31; this->a32 = -this->a32;
    }
    
    void Ann00n0001()
    {
        this->a12 += this->a11; this->a12 = -this->a12; this->a11 = -this->a11;
        this->a22 += this->a21; this->a22 = -this->a22; this->a21 = -this->a21;
        this->a32 += this->a31; this->a32 = -this->a32; this->a31 = -this->a31;
    }
    
    void An0n01000n()
    {
        this->a13 += this->a11; this->a13 = -this->a13; this->a11 = -this->a11;
        this->a23 += this->a21; this->a23 = -this->a23; this->a21 = -this->a21;
        this->a33 += this->a31; this->a33 = -this->a33; this->a31 = -this->a31;
    }
    
    void A1000nn00n()
    {
        this->a13 += this->a12; this->a13 = -this->a13; this->a12 = -this->a12;
        this->a23 += this->a22; this->a23 = -this->a23; this->a22 = -this->a22;
        this->a33 += this->a32; this->a33 = -this->a33; this->a32 = -this->a32;
    }
    
    void Ann001000n()
    {
        this->a12 -= this->a11; this->a11 = -this->a11; this->a13 = -this->a13;
        this->a22 -= this->a21; this->a21 = -this->a21; this->a23 = -this->a23;
        this->a32 -= this->a31; this->a31 = -this->a31; this->a33 = -this->a33;
    }
    
    void An0n0n0001()
    {
        this->a13 -= this->a11; this->a11 = -this->a11; this->a12 = -this->a12;
        this->a23 -= this->a21; this->a21 = -this->a21; this->a22 = -this->a22;
        this->a33 -= this->a31; this->a31 = -this->a31; this->a32 = -this->a32;
    }
    
    void An000nn001()
    {
        this->a13 -= this->a12; this->a12 = -this->a12; this->a11 = -this->a11;
        this->a23 -= this->a22; this->a22 = -this->a22; this->a21 = -this->a21;
        this->a33 -= this->a32; this->a32 = -this->a32; this->a31 = -this->a31;
    }

    
    void A1000010n0()
    {
        R temp;
        temp = this->a13; this->a13 = this->a12; this->a12 = -temp;
        temp = this->a23; this->a23 = this->a22; this->a22 = -temp;
        temp = this->a33; this->a33 = this->a32; this->a32 = -temp;
    }
    
    void A1000100t1(const R& t)
    {
        this->a12 += (this->a13 * t);
        this->a22 += (this->a23 * t);
        this->a32 += (this->a33 * t);
    }
    
    void A100010t01(const R& t)
    {
        this->a11 += (this->a13 * t);
        this->a21 += (this->a23 * t);
        this->a31 += (this->a33 * t);
    }

    void A1000p000p2(const R& p, const R& pp)
    {
        this->a12 *= p;  this->a22 *= p;  this->a32 *= p;
        this->a13 *= pp; this->a23 *= pp; this->a33 *= pp;
    }

    friend std::ostream& operator<<(std::ostream& os, const Isometry<R>& s)
    {
        os << "Matrix(Integers(), 3, (";
        os << s.a11 << "," << s.a12 << "," << s.a13 << ",";
        os << s.a21 << "," << s.a22 << "," << s.a23 << ",";
        os << s.a31 << "," << s.a32 << "," << s.a33 << "))";
        return os;
    }

    static std::vector<std::vector<Isometry<R>>> automorphisms;

    R a11, a12, a13;
    R a21, a22, a23;
    R a31, a32, a33;
};

template<typename R>
std::vector<std::vector<Isometry<R>>> Isometry<R>::automorphisms = {
    {
        { -1, -1, -1,  0,  0,  1,  0,  1,  0 },
        { -1, -1,  0,  0,  1,  0,  0,  0, -1 },
        { -1,  0, -1,  0, -1,  0,  0,  0,  1 },
        { -1,  0,  0,  0,  0, -1,  0, -1,  0 },
        {  1,  0,  1,  0,  0, -1,  0,  1,  0 },
        {  1,  1,  0,  0,  0,  1,  0, -1,  0 },
        {  1,  1,  1,  0, -1,  0,  0,  0, -1 }
    },
    {
        { -1, -1,  0,  0,  1,  0,  0,  0, -1 },
        { -1,  0, -1,  0, -1,  0,  0,  0,  1 },
        {  1,  1,  1,  0, -1,  0,  0,  0, -1 }
    },
    {
        { -1, -1,  0,  0,  1,  0,  0,  0, -1 }
    },
    {
        { -1,  0, -1,  0, -1,  0,  0,  0,  1 }
    },
    {
        { -1,  0,  0,  0, -1, -1,  0,  0,  1 }
    },
    {
        { -1,  0,  0, -1,  1,  0,  0,  0, -1 },
        { -1,  0,  0,  0, -1,  0,  0,  0,  1 },
        { -1,  1,  0, -1,  0,  0,  0,  0,  1 },
        { -1,  1,  0,  0,  1,  0,  0,  0, -1 },
        {  0, -1,  0, -1,  0,  0,  0,  0, -1 },
        {  0, -1,  0,  1, -1,  0,  0,  0,  1 },
        {  0,  1,  0, -1,  1,  0,  0,  0,  1 },
        {  0,  1,  0,  1,  0,  0,  0,  0, -1 },
        {  1, -1,  0,  0, -1,  0,  0,  0, -1 },
        {  1, -1,  0,  1,  0,  0,  0,  0,  1 },
        {  1,  0,  0,  1, -1,  0,  0,  0, -1 }
    },
    {
        { -1,  0,  0,  0, -1,  0,  0,  0,  1 },
        { -1,  1,  0,  0,  1,  0,  0,  0, -1 },
        {  1, -1,  0,  0, -1,  0,  0,  0, -1 }
    },
    {
        {  1, -1,  0,  0, -1,  0,  0,  0, -1 }
    },
    {
        { -1,  0,  0,  0,  1, -1,  0,  0, -1 },
        { -1,  0,  1,  0, -1,  1,  0,  0,  1 },
        {  0, -1,  0, -1,  0,  0,  0,  0, -1 },
        {  0, -1,  1,  1,  0,  0,  0,  0,  1 },
        {  0,  1, -1,  1,  0, -1,  0,  0, -1 },
        {  0,  1,  0, -1,  0,  1,  0,  0,  1 },
        {  1,  0, -1,  0, -1,  0,  0,  0, -1 }
    },
    {
        { -1,  0,  0,  0,  1, -1,  0,  0, -1 },
        { -1,  0,  1,  0, -1,  1,  0,  0,  1 },
        {  1,  0, -1,  0, -1,  0,  0,  0, -1 }
    },
    {
        { -1,  0,  0,  0,  1,  0,  0,  0, -1 },
        { -1,  0,  1,  0, -1,  0,  0,  0,  1 },
        {  1,  0, -1,  0, -1,  0,  0,  0, -1 }
    },
    {
        {  1,  0, -1,  0, -1,  0,  0,  0, -1 }
    },
    {
        { -1,  0,  0,  0, -1,  0,  0, -1,  1 },
        { -1,  0,  0,  0, -1,  1,  0,  0,  1 },
        { -1,  0,  0,  0,  0, -1,  0, -1,  0 },
        { -1,  0,  0,  0,  0,  1,  0,  1,  0 },
        { -1,  0,  0,  0,  1, -1,  0,  0, -1 },
        { -1,  0,  0,  0,  1,  0,  0,  1, -1 },
        {  1,  0,  0,  0, -1,  0,  0,  0, -1 },
        {  1,  0,  0,  0, -1,  1,  0, -1,  0 },
        {  1,  0,  0,  0,  0, -1,  0,  1, -1 },
        {  1,  0,  0,  0,  0,  1,  0, -1,  1 },
        {  1,  0,  0,  0,  1, -1,  0,  1,  0 }
    },
    {
        { -1,  0,  0,  0, -1,  1,  0,  0,  1 },
        { -1,  0,  0,  0,  1, -1,  0,  0, -1 },
        {  1,  0,  0,  0, -1,  0,  0,  0, -1 }
    },
    {
        { -1,  0,  0,  0,  1, -1,  0,  0, -1 }
    },
    {
        { -1,  0,  0, -1,  0,  1, -1,  1,  0 },
        { -1,  0,  0,  0,  0, -1,  0, -1,  0 },
        { -1,  0,  1, -1,  1,  0, -1,  0,  0 },
        { -1,  0,  1,  0, -1,  1,  0,  0,  1 },
        { -1,  1,  0, -1,  0,  0, -1,  0,  1 },
        { -1,  1,  0,  0,  1,  0,  0,  1, -1 },
        {  0, -1,  0, -1,  0,  0,  0,  0, -1 },
        {  0, -1,  0,  1, -1,  0,  0, -1,  1 },
        {  0, -1,  1,  0, -1,  0,  1, -1,  0 },
        {  0, -1,  1,  0,  0,  1, -1,  0,  1 },
        {  0,  0, -1,  0, -1,  0, -1,  0,  0 },
        {  0,  0, -1,  0,  1, -1,  1,  0, -1 },
        {  0,  0,  1, -1,  0,  1,  0, -1,  1 },
        {  0,  0,  1,  1,  0,  0,  0,  1,  0 },
        {  0,  1, -1, -1,  1,  0,  0,  1,  0 },
        {  0,  1, -1,  1,  0, -1,  0,  0, -1 },
        {  0,  1,  0,  0,  0,  1,  1,  0,  0 },
        {  0,  1,  0,  0,  1, -1, -1,  1,  0 },
        {  1, -1,  0,  0, -1,  1,  0, -1,  0 },
        {  1, -1,  0,  1,  0, -1,  1,  0,  0 },
        {  1,  0, -1,  0,  0, -1,  0,  1, -1 },
        {  1,  0, -1,  1,  0,  0,  1, -1,  0 },
        {  1,  0,  0,  1, -1,  0,  1,  0, -1 }
    },
    {
        { -1,  0,  0, -1,  0,  1, -1,  1,  0 },
        { -1,  0,  1,  0, -1,  1,  0,  0,  1 },
        {  0, -1,  0, -1,  0,  0,  0,  0, -1 },
        {  0, -1,  1,  0, -1,  0,  1, -1,  0 },
        {  0,  1, -1,  1,  0, -1,  0,  0, -1 },
        {  0,  1,  0,  0,  1, -1, -1,  1,  0 },
        {  1,  0, -1,  1,  0,  0,  1, -1,  0 }
    },
    {
        { -1,  0,  1,  0, -1,  1,  0,  0,  1 },
        {  0, -1,  0, -1,  0,  0,  0,  0, -1 },
        {  0,  1, -1,  1,  0, -1,  0,  0, -1 }
    },
    {
        { -1,  0,  0,  0,  0, -1,  0, -1,  0 },
        { -1,  0,  1,  0, -1,  1,  0,  0,  1 },
        { -1,  1,  0,  0,  1,  0,  0,  1, -1 },
        {  1, -1,  0,  0, -1,  1,  0, -1,  0 },
        {  1,  0, -1,  0,  0, -1,  0,  1, -1 }
    },
    {
        { -1,  0,  1,  0, -1,  1,  0,  0,  1 }
    },
    {
        { -1,  0,  0,  0, -1,  0,  0,  0,  1 },
        { -1,  0,  0,  0,  0, -1,  0, -1,  0 },
        { -1,  0,  0,  0,  0,  1,  0,  1,  0 },
        { -1,  0,  0,  0,  1,  0,  0,  0, -1 },
        {  0, -1,  0, -1,  0,  0,  0,  0, -1 },
        {  0, -1,  0,  0,  0, -1,  1,  0,  0 },
        {  0, -1,  0,  0,  0,  1, -1,  0,  0 },
        {  0, -1,  0,  1,  0,  0,  0,  0,  1 },
        {  0,  0, -1, -1,  0,  0,  0,  1,  0 },
        {  0,  0, -1,  0, -1,  0, -1,  0,  0 },
        {  0,  0, -1,  0,  1,  0,  1,  0,  0 },
        {  0,  0, -1,  1,  0,  0,  0, -1,  0 },
        {  0,  0,  1, -1,  0,  0,  0, -1,  0 },
        {  0,  0,  1,  0, -1,  0,  1,  0,  0 },
        {  0,  0,  1,  0,  1,  0, -1,  0,  0 },
        {  0,  0,  1,  1,  0,  0,  0,  1,  0 },
        {  0,  1,  0, -1,  0,  0,  0,  0,  1 },
        {  0,  1,  0,  0,  0, -1, -1,  0,  0 },
        {  0,  1,  0,  0,  0,  1,  1,  0,  0 },
        {  0,  1,  0,  1,  0,  0,  0,  0, -1 },
        {  1,  0,  0,  0, -1,  0,  0,  0, -1 },
        {  1,  0,  0,  0,  0, -1,  0,  1,  0 },
        {  1,  0,  0,  0,  0,  1,  0, -1,  0 }
    },
    {
        { -1, -1, -1,  0,  0,  1,  0,  1,  0 },
        { -1, -1, -1,  0,  1,  0,  1,  0,  0 },
        { -1, -1, -1,  1,  0,  0,  0,  0,  1 },
        { -1,  0,  0,  0, -1,  0,  1,  1,  1 },
        { -1,  0,  0,  0,  0, -1,  0, -1,  0 },
        { -1,  0,  0,  1,  1,  1,  0,  0, -1 },
        {  0, -1,  0, -1,  0,  0,  0,  0, -1 },
        {  0, -1,  0,  0,  0, -1,  1,  1,  1 },
        {  0, -1,  0,  1,  1,  1, -1,  0,  0 },
        {  0,  0, -1, -1,  0,  0,  1,  1,  1 },
        {  0,  0, -1,  0, -1,  0, -1,  0,  0 },
        {  0,  0, -1,  1,  1,  1,  0, -1,  0 },
        {  0,  0,  1, -1, -1, -1,  1,  0,  0 },
        {  0,  0,  1,  0,  1,  0, -1, -1, -1 },
        {  0,  0,  1,  1,  0,  0,  0,  1,  0 },
        {  0,  1,  0, -1, -1, -1,  0,  0,  1 },
        {  0,  1,  0,  0,  0,  1,  1,  0,  0 },
        {  0,  1,  0,  1,  0,  0, -1, -1, -1 },
        {  1,  0,  0, -1, -1, -1,  0,  1,  0 },
        {  1,  0,  0,  0,  0,  1, -1, -1, -1 },
        {  1,  1,  1, -1,  0,  0,  0, -1,  0 },
        {  1,  1,  1,  0, -1,  0,  0,  0, -1 },
        {  1,  1,  1,  0,  0, -1, -1,  0,  0 }
    },
    {
        { -1,  0,  0,  0,  0, -1,  0, -1,  0 },
        {  0, -1,  0, -1,  0,  0,  0,  0, -1 },
        {  0,  0, -1,  0, -1,  0, -1,  0,  0 },
        {  0,  0,  1,  1,  0,  0,  0,  1,  0 },
        {  0,  1,  0,  0,  0,  1,  1,  0,  0 }
    },
    {
        { -1,  0,  0,  0, -1,  0,  0,  0,  1 },
        { -1,  0,  0,  0,  1,  0,  0,  0, -1 },
        {  0, -1,  0, -1,  0,  0,  0,  0, -1 },
        {  0, -1,  0,  1,  0,  0,  0,  0,  1 },
        {  0,  1,  0, -1,  0,  0,  0,  0,  1 },
        {  0,  1,  0,  1,  0,  0,  0,  0, -1 },
        {  1,  0,  0,  0, -1,  0,  0,  0, -1 }
    },
    {
        { -1,  0,  0,  0, -1,  0,  0,  0,  1 },
        {  0, -1,  0, -1,  0,  0,  0,  0, -1 },
        {  0,  1,  0,  1,  0,  0,  0,  0, -1 }
    },
    {
        { -1, -1, -1,  1,  0,  0,  0,  0,  1 },
        { -1,  0,  0,  1,  1,  1,  0,  0, -1 },
        {  0, -1,  0, -1,  0,  0,  0,  0, -1 },
        {  0,  1,  0, -1, -1, -1,  0,  0,  1 },
        {  1,  1,  1,  0, -1,  0,  0,  0, -1 }
    },
    {
        {  0, -1,  0, -1,  0,  0,  0,  0, -1 }
    },
    {
        { -1,  0,  0,  0, -1,  0,  0,  0,  1 },
        { -1,  0,  0,  0,  0, -1,  0, -1,  0 },
        { -1,  0,  0,  0,  0,  1,  0,  1,  0 },
        { -1,  0,  0,  0,  1,  0,  0,  0, -1 },
        {  1,  0,  0,  0, -1,  0,  0,  0, -1 },
        {  1,  0,  0,  0,  0, -1,  0,  1,  0 },
        {  1,  0,  0,  0,  0,  1,  0, -1,  0 }
    },
    {
        { -1,  0,  0,  0,  0, -1,  0, -1,  0 },
        { -1,  0,  0,  0,  0,  1,  0,  1,  0 },
        {  1,  0,  0,  0, -1,  0,  0,  0, -1 }
    },
    {
        { -1, -1, -1,  0,  0,  1,  0,  1,  0 },
        { -1,  0,  0,  0,  0, -1,  0, -1,  0 },
        {  1,  1,  1,  0, -1,  0,  0,  0, -1 }
    },
    {
        { -1, -1, -1,  0,  0,  1,  0,  1,  0 },
        { -1,  0,  0,  0,  0, -1,  0, -1,  0 },
        {  1,  1,  1,  0, -1,  0,  0,  0, -1 }
    },
    {
        { -1,  0,  0, -1,  0,  1, -1,  1,  0 },
        { -1,  0,  0,  0,  0, -1,  0, -1,  0 },
        {  0, -1,  1,  0, -1,  0,  1, -1,  0 },
        {  0, -1,  1,  0,  0,  1, -1,  0,  1 },
        {  0,  1, -1, -1,  1,  0,  0,  1,  0 },
        {  0,  1, -1,  1,  0, -1,  0,  0, -1 },
        {  1,  0,  0,  1, -1,  0,  1,  0, -1 }
    },
    {
        { -1,  0,  0,  0,  0, -1,  0, -1,  0 }
    },
    {
        { -1,  0,  0,  0, -1,  0,  0,  0,  1 },
        { -1,  0,  0,  0,  1,  0,  0,  0, -1 },
        {  1,  0,  0,  0, -1,  0,  0,  0, -1 }
    },
    {
        { -1,  0,  0,  0, -1,  0,  0,  0,  1 }
    },
    {
        { -1,  0,  0,  0,  1,  0,  0,  0, -1 }
    },
    {
        {  1,  0,  0,  0, -1,  0,  0,  0, -1 }
    },
    {
        {  1,  1,  1,  0, -1,  0,  0,  0, -1 }
    },
    {
        {  1,  1,  1,  0, -1,  0,  0,  0, -1 }
    },
    {
        { -1,  0,  0, -1,  0,  1, -1,  1,  0 },
        {  0, -1,  1,  0, -1,  0,  1, -1,  0 },
        {  0,  1, -1,  1,  0, -1,  0,  0, -1 }
    },
    {
        {  0,  1, -1,  1,  0, -1,  0,  0, -1 }
    },
    {
    }
};

template<>
void Z_Isometry::set_identity(void);

template<>
Z_Isometry Z_Isometry::inverse(const Z& p) const;

template<>
Z_Isometry Z_Isometry::operator*(const Z_Isometry& s) const;

template<>
void Z_Isometry::A101011001();

template<>
void Z_Isometry::A1t0010001(const Z& t);

template<>
void Z_Isometry::A10001t001(const Z& t);

template<>
void Z_Isometry::A10t010001(const Z& t);

template<>
void Z_Isometry::A0n0n0000n();

template<>
void Z_Isometry::An0000n0n0();

template<>
void Z_Isometry::An00010001();

template<>
void Z_Isometry::A1000n0001();

template<>
void Z_Isometry::A10001000n();

template<>
void Z_Isometry::An010n1001();

template<>
void Z_Isometry::Ann00n0001();

template<>
void Z_Isometry::An0n01000n();

template<>
void Z_Isometry::A1000nn00n();

template<>
void Z_Isometry::Ann001000n();

template<>
void Z_Isometry::An0n0n0001();

template<>
void Z_Isometry::An000nn001();

template<>
void Z_Isometry::A1000010n0();

template<>
void Z_Isometry::A1000100t1(const Z& t);

template<>
void Z_Isometry::A100010t01(const Z& t);

template<>
void Z_Isometry::A1000p000p2(const Z& p, const Z& pp);

#endif // __ISOMETRY_H_
