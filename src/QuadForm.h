#ifndef __QUAD_FORM_H_
#define __QUAD_FORM_H_

#include "birch.h"
#include "birch_util.h"

template<typename R>
class QuadForm
{
public:
    QuadForm() = default;

    QuadForm(const R& a, const R& b, const R& c,
             const R& f, const R& g, const R& h)
    {
        this->a_ = a; this->b_ = b; this->c_ = c;
        this->f_ = f; this->g_ = g; this->h_ = h;
    }

    const R& a(void) const { return this->a_; }
    const R& b(void) const { return this->b_; }
    const R& c(void) const { return this->c_; }
    const R& f(void) const { return this->f_; }
    const R& g(void) const { return this->g_; }
    const R& h(void) const { return this->h_; }

    R discriminant(void) const
    {
        return this->a_ * (4 * this->b_ * this->c_ - this->f_ * this->f_) -
            this->b_ * this->g_ * this->g_ +
            this->h_ * (this->f_ * this->g_ - this->c_ * this->h_);
    }

    bool operator==(const QuadForm<R>& q) const
    {
        return this->a_ == q.a_ && this->b_ == q.b_ && this->c_ == q.c_ &&
            this->f_ == q.f_ && this->g_ == q.g_ && this->h_ == q.h_;
    }

    W64 hash_value(void) const;

    R evaluate(const R& x, const R& y, const R& z) const
    {
        return x * (this->a_ * x + this->g_ * z + this->h_ * y) +
            y * (this->b_ * y + this->f_ * z) + z * z * this->c_;
    }

    R evaluate(const Vector3<R>& vec) const
    {
        return this->evaluate(vec.x, vec.y, vec.z);
    }

    template<typename S, typename T>
    QuadFormFp<S,T> mod(std::shared_ptr<Fp<S,T>> GF) const
    {
        QuadFormFp<S,T> q(GF->mod(this->a_), GF->mod(this->b_), GF->mod(this->c_),
                          GF->mod(this->f_), GF->mod(this->g_), GF->mod(this->h_), GF);
        return q;
    }

    static Z_QuadForm get_quad_form(const std::vector<PrimeSymbol<R>>& primes)
    {
        static_assert( std::is_same<R,Z>::value, "Implemented only for arbitrary precision types." );
        return Z_QuadForm(); // Make the compiler happy.
    }

    static int border(const QuadForm<R>& q, int n)
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
                return 0;
        }
    }

    static int num_automorphisms(const QuadForm<R>& q)
    {
        if (border(q, 1))
        {
            if (border(q, 2))
            {
                if (border(q, 14))
                {
                    if (border(q, 9))
                        return 16;
                    else
                        return 8;
                }
            }
            else
                return 4;
        }

        if (border(q, 2))
            return 4;

        if (border(q, 3))
            return 4;

        if (border(q, 4))
        {
            if (border(q, 10))
            {
                if (border(q, 8))
                    return 24;
                else
                    return 8;
            }
            else
                return 4;
        }

        if (border(q, 5))
        {
            if (border(q, 6))
            {
                if (border(q, 7))
                {
                    if (border(q, 8))
                    {
                        if (border(q, 15))
                            return 16;
                    }
                    else
                        return 8;
                }
            }
            else if (border(q, 11))
                return 8;
            else
                return 4;
        }

        if (border(q, 6))
        {
            if (border(q, 12))
            {
                if (border(q, 9))
                    return 24;
                else
                    return 8;
            }
            else
                return 4;
        }

        if (border(q, 7))
        {
            if (border(q, 8) && border(q, 15))
            {
                if (border(q, 16))
                {
                    if (border(q, 9))
                        return 48;
                    else
                        return 16;
                }
                else
                    return 8;
            }
            else if (border(q, 9))
                return 12;
            else
                return 4;
        }

        if (border(q, 8))
        {
            if (border(q, 9))
            {
                if (border(q, 10) && border(q, 11) && border(q, 12))
                    return 48;
                else if (border(q, 13) && border(q, 14))
                    return 48;
                else
                    return 12;
            }
            else if (border(q, 10))
            {
                if (border(q, 11) && border(q, 12))
                    return 16;
                else
                    return 8;
            }
            else if (border(q, 14))
                return 12;
            else
                return 4;
        }

        if (border(q, 9))
        {
            if (border(q, 12))
            {
                if (border(q, 10) && border(q, 11))
                    return 16;
                else
                    return 8;
            }
            else if (border(q, 14))
            {
                if (border(q, 13))
                    return 8;
                else
                    return 8;
            }
            else if (border(q, 15))
                return 16;
            else
                return 4;
        }

        if (border(q, 10))
        {
            if (border(q, 11) && border(q, 12))
                return 8;
            else
                return 4;
        }

        if (border(q, 11))
            return 4;

        if (border(q, 12))
            return 4;

        if (border(q, 13) && border(q, 14))
            return 4;

        if (border(q, 14))
            return 4;

        if (border(q, 15))
        {
            if (border(q, 16))
                return 8;
            else
                return 4;
        }

        return 2;
    }

    static const std::vector<Isometry<R>>& proper_automorphisms(const QuadForm<R>& q)
    {
        if (border(q, 1))
        {
            if (border(q, 2))
            {
                if (border(q, 14))
                {
                    if (border(q, 9))
                    {
                        return Isometry<R>::automorphisms[0];
                    }
                    else
                    {
                        return Isometry<R>::automorphisms[1];
                    }
                }
            }
            else
            {
                return Isometry<R>::automorphisms[2];
            }
        }

        if (border(q, 2))
        {
            return Isometry<R>::automorphisms[3];
        }

        if (border(q, 3))
        {
            return Isometry<R>::automorphisms[4];
        }

        if (border(q, 4))
        {
            if (border(q, 10))
            {
                if (border(q, 8))
                {
                    return Isometry<R>::automorphisms[5];
                }
                else
                {
                    return Isometry<R>::automorphisms[6];
                }
            }
            else
            {
                return Isometry<R>::automorphisms[7];
            }
        }

        if (border(q, 5))
        {
            if (border(q, 6))
            {
                if (border(q, 7))
                {
                    if (border(q, 8))
                    {
                        if (border(q, 15))
                        {
                            return Isometry<R>::automorphisms[8];
                        }
                    }
                    else
                    {
                        return Isometry<R>::automorphisms[9];
                    }
                }
            }
            else if (border(q, 11))
            {
                return Isometry<R>::automorphisms[10];
            }
            else
            {
                return Isometry<R>::automorphisms[11];
            }
        }

        if (border(q, 6))
        {
            if (border(q, 12))
            {
                if (border(q, 9))
                {
                    return Isometry<R>::automorphisms[12];
                }
                else
                {
                    return Isometry<R>::automorphisms[13];
                }
            }
            else
            {
                return Isometry<R>::automorphisms[14];
            }
        }

        if (border(q, 7))
        {
            if (border(q, 8) && border(q, 15))
            {
                if (border(q, 16))
                {
                    if (border(q, 9))
                    {
                        return Isometry<R>::automorphisms[15];
                    }
                    else
                    {
                        return Isometry<R>::automorphisms[16];
                    }
                }
                else
                {
                    return Isometry<R>::automorphisms[17];
                }
            }
            else if (border(q, 9))
            {
                return Isometry<R>::automorphisms[18];
            }
            else
            {
                return Isometry<R>::automorphisms[19];
            }
        }

        if (border(q, 8))
        {
            if (border(q, 9))
            {
                if (border(q, 10) && border(q, 11) && border(q, 12))
                {
                    return Isometry<R>::automorphisms[20];
                }
                else if (border(q, 13) && border(q, 14))
                {
                    return Isometry<R>::automorphisms[21];
                }
                else
                {
                    return Isometry<R>::automorphisms[22];
                }
            }
            else if (border(q, 10))
            {
                if (border(q, 11) && border(q, 12))
                {
                    return Isometry<R>::automorphisms[23];
                }
                else
                {
                    return Isometry<R>::automorphisms[24];
                }
            }
            else if (border(q, 14))
            {
                return Isometry<R>::automorphisms[25];
            }
            else
            {
                return Isometry<R>::automorphisms[26];
            }
        }

        if (border(q, 9))
        {
            if (border(q, 12))
            {
                if (border(q, 10) && border(q, 11))
                {
                    return Isometry<R>::automorphisms[27];
                }
                else
                {
                    return Isometry<R>::automorphisms[28];
                }
            }
            else if (border(q, 14))
            {
                if (border(q, 13))
                {
                    return Isometry<R>::automorphisms[29];
                }
                else
                {
                    return Isometry<R>::automorphisms[30];
                }
            }
            else if (border(q, 15))
            {
                return Isometry<R>::automorphisms[31];
            }
            else
            {
                return Isometry<R>::automorphisms[32];
            }
        }

        if (border(q, 10))
        {
            if (border(q, 11) && border(q, 12))
            {
                return Isometry<R>::automorphisms[33];
            }
            else
            {
                return Isometry<R>::automorphisms[34];
            }
        }

        if (border(q, 11))
        {
            return Isometry<R>::automorphisms[35];
        }

        if (border(q, 12))
        {
            return Isometry<R>::automorphisms[36];
        }

        if (border(q, 13) && border(q, 14))
        {
            return Isometry<R>::automorphisms[37];
        }

        if (border(q, 14))
        {
            return Isometry<R>::automorphisms[38];
        }

        if (border(q, 15))
        {
            if (border(q, 16))
            {
                return Isometry<R>::automorphisms[39];
            }
            else
            {
                return Isometry<R>::automorphisms[40];
            }
        }

        return Isometry<R>::automorphisms[41];
    }

    static QuadForm<R> reduce(const QuadForm<R>& q, Isometry<R>& s)
    {
        R a = q.a_;
        R b = q.b_;
        R c = q.c_;
        R f = q.f_;
        R g = q.g_;
        R h = q.h_;

        int flag = 1;
        while (flag)
        {
            R t = a + b + f + g + h;
            if (t < 0)
            {
                s.A101011001();
                c += t;
                f += (h + b + b);
                g += (h + a + a);
            }

            if (a >= h)
            {
                t = birch_util::dumb_div<R>(a-h, a+a);
            }
            else
            {
                t = -birch_util::dumb_div<R>(a+h-1, a+a);
            }
            if (t != 0)
            {
                s.A1t0010001(t);
                R temp = a * t;
                h += temp;
                b += h * t;
                f += g * t;
                h += temp;
            }

            if (b >= f)
            {
                t = birch_util::dumb_div<R>(b-f, b+b);
            }
            else
            {
                t = -birch_util::dumb_div<R>(b+f-1, b+b);
            }
            if (t != 0)
            {
                s.A10001t001(t);
                R temp = b * t;
                f += temp;
                c += f * t;
                g += h * t;
                f += temp;
            }

            if (a >= g)
            {
                t = birch_util::dumb_div<R>(a-g, a+a);
            }
            else
            {
                t = -birch_util::dumb_div<R>(a+g-1, a+a);
            }
            if (t != 0)
            {
                s.A10t010001(t);
                R temp = a * t;
                g += temp;
                c += g * t;
                f += h * t;
                g += temp;
            }

            if (a > b || (a == b && abs(f) > abs(g)))
            {
                s.A0n0n0000n();
                t = a; a = b; b = t;
                t = f; f = g; g = t;
            }

            if (b > c || (b == c && abs(g) > abs(h)))
            {
                s.An0000n0n0();
                t = b; b = c; c = t;
                t = g; g = h; h = t;
            }

            if (a > b || (a == b && abs(f) > abs(g)))
            {
                s.A0n0n0000n();
                t = a; a = b; b = t;
                t = f; f = g; g = t;
            }

            int fgh = (f != 0 && g != 0 && h != 0);
            if (fgh)
            {
                if (f < 0) fgh = !fgh;
                if (g < 0) fgh = !fgh;
                if (h < 0) fgh = !fgh;
            }

            if (fgh)
            {
                if (f < 0)
                {
                    s.An00010001();
                    f = -f;
                }

                if (g < 0)
                {
                    s.A1000n0001();
                    g = -g;
                }

                if (h < 0)
                {
                    s.A10001000n();
                    h = -h;
                }
            }
            else
            {
                int s1 = f > 0;
                int s2 = g > 0;
                int s3 = h > 0;

                if ((s1+s2+s3) % 2 == 1)
                {
                    if (f == 0) s1 = 1;
                    else
                    {
                        if (g == 0) s2 = 1;
                        else if (h == 0) s3 = 1;
                    }
                }

                if (s1 == 1)
                {
                    s.An00010001();
                    f = -f;
                }

                if (s2 == 1)
                {
                    s.A1000n0001();
                    g = -g;
                }

                if (s3 == 1)
                {
                    s.A10001000n();
                    h = -h;
                }
            }

            flag = !(abs(f) <= b && abs(g) <= a &&
                    abs(h) <= a && a+b+f+g+h >= 0);
        }

        if (a + b + f + g + h == 0 &&
            a + a + g + g + h > 0)
        {
            s.An010n1001();
            c += a + b + f + g + h;
            f += h + b + b; f = -f;
            g += h + a + a; g = -g;
        }

        if (a == -h && g != 0)
        {
            s.Ann00n0001();
            f += g; f = -f;
            g = -g;
            h = -h;
        }

        if (a == -g && h != 0)
        {
            s.An0n01000n();
            f += h; f = -f;
            h = -h;
            g += (2*a);
        }

        if (b == -f && h != 0)
        {
            s.A1000nn00n();
            g += h; g = -g;
            h = -h;
            f += (2*b);
        }

        if (a == h && g > f + f)
        {
            s.Ann001000n();
            f = g - f;
        }

        if (a == g && h > f + f)
        {
            s.An0n0n0001();
            f = h - f;
        }

        if (b == f && h > g + g)
        {
            s.An000nn001();
            g = h - g;
        }

        if (a == b && abs(f) > abs(g))
        {
            s.A0n0n0000n();
            R t;
            t = a; a = b; b = t;
            t = g; g = f; f = t;
        }

        if (b == c && abs(g) > abs(h))
        {
            s.An0000n0n0();
            R t;
            t = g; g = h; h = t;
        }

        if (a == b && abs(f) > abs(g))
        {
            s.A0n0n0000n();
            R t;
            t = g; g = f; f = t;
        }

        return QuadForm<R>(a, b, c, f, g, h);
    }

    friend std::ostream& operator<<(std::ostream& os, const QuadForm<R>& q)
    {
        os << "QuadForm(" << q.a_ << "," << q.b_ << "," << q.c_ << ","
            << q.f_ << "," << q.g_ << "," << q.h_ << ")";
        return os;
    }

protected:
    R a_, b_, c_, f_, g_, h_;
};

template<typename R, typename S>
class QuadFormFp : public QuadForm<R>
{
public:
    QuadFormFp(const R& a, const R& b, const R& c,
               const R& f, const R& g, const R& h,
               std::shared_ptr<Fp<R,S>> GF) :
        QuadForm<R>(GF->mod(a), GF->mod(b), GF->mod(c),
                    GF->mod(f), GF->mod(g), GF->mod(h))
    {
        this->GF = GF;
    }

    const std::shared_ptr<Fp<R,S>>& field(void) const
    {
        return this->GF;
    }

    R discriminant(void) const
    {
        R res = GF->mul(this->b_, this->c_);    // bc
        res = GF->add(res, res);                // 2bc
        res = GF->add(res, res);                // 4bc
        R temp = GF->mul(this->f_, this->f_);   // ff
        res = GF->sub(res, temp);               // 4bc-ff
        res = GF->mul(this->a_, res);           // a(4bc-ff)

        temp = GF->mul(this->g_, this->g_);     // gg
        temp = GF->mul(this->b_, temp);         // bgg
        res = GF->sub(res, temp);               // a(4bc-ff)-bgg

        temp = GF->mul(this->f_, this->g_);     // fg
        R temp2 = GF->mul(this->c_, this->h_);  // ch
        temp = GF->sub(temp, temp2);            // fg-ch
        temp = GF->mul(this->h_, temp);         // h(fg-ch)
        res = GF->add(res, temp);               // a(4bc-ff)-bgg+h(fg-ch)

        return res;
    }

    R evaluate(const R& x, const R& y, const R& z) const
    {
        R res = GF->mul(this->a_, x);   // ax
        R temp = GF->mul(this->g_, z);  // gz
        res = GF->add(res, temp);       // ax+gz
        temp = GF->mul(this->h_, y);    // hy
        res = GF->add(res, temp);       // ax+gz+hy
        res = GF->mul(x, res);          // x(ax+gz+hy)

        temp = GF->mul(this->b_, y);    // by
        R temp2 = GF->mul(this->f_, z); // fz
        temp = GF->add(temp, temp2);    // by+fz
        temp = GF->mul(y, temp);        // y(by+fz)
        res = GF->add(res, temp);       // x(ax+gz+hy)+y(by+fz)

        temp = GF->mul(this->c_, z);    // cz
        temp = GF->mul(temp, z);        // czz
        res = GF->add(res, temp);       // x(ax+gz+hy)+y(by+fz)+czz

        return res;
    }

    R evaluate(const Vector3<R>& vec) const
    {
        return this->evaluate(vec.x, vec.y, vec.z);
    }

    Vector3<R> isotropic_vector(void) const
    {
        Vector3<R> vec = {0,0,0};

        if (GF->prime() == 2) return this->isotropic_vector_p2();

        while (1)
        {
            R r = 0;
            R alpha = 0;
            R beta = 0;

            while (alpha == 0)
            {
                r = GF->random();                   // r
                beta = GF->mul(this->b_, r);        // br
                alpha = GF->add(beta, this->h_);    // br+h
                alpha = GF->mul(alpha, r);          // (br+h)r
                alpha = GF->add(alpha, this->a_);   // (br+h)r+a = Q(1,r,0)
                alpha = GF->mod(alpha);
            }

            R s = GF->random();

            beta = GF->add(beta, beta);         // 2br
            beta = GF->add(beta, this->h_);     // 2br+h
            beta = GF->mul(beta, s);            // (2br+h)s
            R temp = GF->mul(this->f_, r);      // fs
            beta = GF->add(beta, temp);         // (2br+h)s+fr = (dQ/dy)(s,rs,r)
            beta = GF->add(beta, this->g_);     // (2br+h)s+fr+g

            R gamma = GF->mul(this->b_, s);     // bs
            gamma = GF->add(gamma, this->f_);   // bs+f
            gamma = GF->mul(gamma, s);          // (bs+f)s
            gamma = GF->add(gamma, this->c_);   // (bs+f)s+c = Q(0,s,1)

            R disc = GF->mul(beta, beta);
            gamma = GF->mul(gamma, alpha);
            gamma = GF->add(gamma, gamma);
            gamma = GF->add(gamma, gamma);
            disc = GF->sub(disc, gamma);

            if (GF->legendre(disc) >= 0)
            {
                R root = GF->sqrt(disc);

                root = GF->sub(root, beta);
                alpha = GF->add(alpha, alpha);
                alpha = GF->mod(alpha);

                vec.x = GF->mul(root, GF->inverse(alpha));
                vec.y = GF->mul(r, vec.x);
                vec.y = GF->add(vec.y, s);
                vec.z = 1;

                return vec;
            }
        }

        return vec;
    }

private:
    std::shared_ptr<Fp<R,S>> GF;

    // To avoid unnecessary computation, we encode each of the three 2-isotropic
    // vectors as a coordinate of the return vector. Special care must be taken
    // to obtain the actual isotropic vectors when needed.
    Vector3<R> isotropic_vector_p2(void) const
    {
        Vector3<R> vec = {0,0,0};
        R temp[3] = {0,0,0};

        int index = 0;

        if (this->c_ == 0) temp[index++] = 1;
        if (this->b_ == 0) temp[index++] = 2;
        if ((this->b_ ^ this->f_ ^ this->c_) == 0) temp[index++] = 3;
        if (this->a_ == 0) temp[index++] = 4;
        if ((this->a_ ^ this->g_ ^ this->c_) == 0) temp[index++] = 5;
        if ((this->a_ ^ this->h_ ^ this->b_) == 0) temp[index++] = 6;
        if ((this->a_ ^ this->b_ ^ this->c_ ^ this->f_ ^ this->g_ ^ this->h_) == 0) temp[index++] = 7;

        vec.x = temp[0];
        vec.y = temp[1];
        vec.z = temp[2];

        return vec;
    }
};

namespace std
{
    template<typename R>
    struct hash<QuadForm<R>>
    {
        Z64 operator()(const QuadForm<R>& q) const
        {
            return q.hash_value();
        }
    };

    template<typename R>
    struct hash<GenusRep<R>>
    {
        Z64 operator()(const GenusRep<R>& rep) const
        {
            return rep.q.hash_value();
        }
    };
}

template<typename R>
std::ostream& operator<<(std::ostream& os, const Vector3<R>& vec)
{
    os << "Vector(" << vec.x << "," << vec.y << "," << vec.z << ")";
    return os;
}

template<typename R>
bool operator==(const Vector3<R>& vec1, const Vector3<R>& vec2)
{
    return vec1.x == vec2.x && vec1.y == vec2.y && vec1.z == vec2.z;
}

template<typename R>
Vector3<R> operator+(const Vector3<R>& a, const Vector3<R>& b)
{
    Vector3<R> res;
    res.x = a.x + b.x;
    res.y = a.y + b.y;
    res.z = a.z + b.z;
    return res;
}

template<>
Z_QuadForm Z_QuadForm::get_quad_form(const std::vector<Z_PrimeSymbol>& primes);

#endif // __QUAD_FORM_H_
