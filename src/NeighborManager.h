#ifndef __NEIGHBOR_MANAGER_H_
#define __NEIGHBOR_MANAGER_H_

#include "birch.h"
#include "QuadForm.h"
#include "Isometry.h"
#include "Fp.h"

template<typename R, typename S, typename T>
class NeighborManager
{
public:
    NeighborManager(const QuadForm<T>& q, std::shared_ptr<Fp<R,S>> GF)
    {
        this->q = q;
        this->disc = q.discriminant();

        QuadFormFp<R,S> qp = q.mod(GF);

        this->a = qp.a();
        this->b = qp.b();
        this->c = qp.c();
        this->f = qp.f();
        this->g = qp.g();
        this->h = qp.h();
        this->GF = GF;
        this->vec = qp.isotropic_vector();

        #ifdef DEBUG
        R prime = GF->prime();
        if (prime != 2) assert( qp.evaluate(vec) % prime == 0 );
        #endif

        R temp = GF->mul(a, vec.x);     // ax
        temp = GF->add(temp, temp);     // 2ax
        a0 = GF->sub(a, temp);          // a-2ax
        a0 = GF->sub(a0, g);            // a-2ax-g
        temp = GF->mul(h, vec.y);       // hy
        a0 = GF->sub(a0, temp);         // a-2ax-g-hy

        temp = GF->mul(b, vec.y);       // by
        temp = GF->add(temp, temp);     // 2by
        delta = GF->sub(b, temp);       // b-2by
        delta = GF->sub(delta, f);      // b-2by-f
        delta = GF->add(delta, h);      // b-2by-f+h
        temp = GF->mul(h, vec.x);       // hx
        delta = GF->sub(delta, temp);   // b-2by-f+h-hx
    }

    Vector3<R> isotropic_vector(R t) const
    {
        Vector3<R> res;

        R p = GF->prime();

        if (p == 2) return this->isotropic_vector_p2(t);

        if (t == p)
        {
            R at = GF->sub(delta, h);
            if (at >= p) at = GF->mod(at);
            if (at == 0)
            {
                res.x = vec.x;
                res.y = vec.y ? vec.y-1 : p-1;
                res.z = 1;
            }
            else
            {
                if (b == 0)
                {
                    res.x = 0;
                    res.y = 1;
                    res.z = 0;
                }
                else
                {
                    R inv = GF->mul(GF->inverse(b), at);
                    res.x = vec.x;
                    res.y = GF->add(vec.y, inv-1);
                    res.z = 1;
                }
            }
        }
        else
        {
            R at = a0;
            if (t == 1) at = GF->add(at, delta);
            else if (t >= 2)
            {
                R temp = GF->mul(t-1, b);
                temp = GF->add(temp, delta);
                temp = GF->mul(temp, t);
                at = GF->add(temp, at);
            }
            if (at >= p) at = GF->mod(at);

            if (at == 0)
            {
                res.x = vec.x ? vec.x-1 : p-1;
                res.y = GF->sub(vec.y, t);
                res.z = 1;
            }
            else
            {
                R ct = GF->mul(b, t);
                ct = GF->add(ct, h);
                ct = GF->mul(ct, t);
                ct = GF->add(ct, a);
                if (ct >= p) ct = GF->mod(ct);

                if (ct == 0)
                {
                    if (t == 0)
                    {
                        res.x = 1;
                        res.y = 0;
                        res.z = 0;
                    }
                    else
                    {
                        res.x = GF->inverse(t);
                        res.y = 1;
                        res.z = 0;
                    }
                }
                else
                {
                    R inv = GF->mul(at, GF->inverse(ct));
                    res.x = GF->add(vec.x, inv-1);
                    res.y = GF->sub(vec.y, t);
                    res.y = GF->add(res.y, GF->mul(t, inv));
                    res.z = 1;
                }
            }
        }

        #ifdef DEBUG
        QuadFormFp<R,S> qp = this->q.mod(GF);
        assert( qp.evaluate(res) % this->GF->prime() == 0 );
        #endif

        return res;
    }

    inline GenusRep<T> get_reduced_neighbor_rep(R t) const
    {
        GenusRep<T> rep;
        rep.q = this->get_neighbor(t, rep.s);
        rep.q = QuadForm<T>::reduce(rep.q, rep.s);
        return rep;
    }

    Vector3<R> transform_vector(const GenusRep<T>& dst, Vector3<R> src)
    {
        Vector3<T> temp;
        temp.x = GF->mod(src.x);
        temp.y = GF->mod(src.y);
        temp.z = GF->mod(src.z);

        R p = GF->prime();

        Isometry<T> sinv = dst.s.inverse(p);
        temp = sinv * temp;

        #ifdef DEBUG
        assert( temp.x % p == 0 );
        assert( temp.y % p == 0 );
        assert( temp.z % p == 0 );
        #endif

        temp.x /= p;
        temp.y /= p;
        temp.z /= p;

        Vector3<R> vec;
        vec.x = GF->mod(temp.x);
        vec.y = GF->mod(temp.y);
        vec.z = GF->mod(temp.z);
        if (vec.z != 0)
        {
            R inv = GF->inverse(vec.z);
            vec.x = GF->mod(GF->mul(vec.x, inv));
            vec.y = GF->mod(GF->mul(vec.y, inv));
            vec.z = 1;
        }
        else if (vec.y != 0)
        {
            R inv = GF->inverse(vec.y);
            vec.x = GF->mod(GF->mul(vec.x, inv));
            vec.y = 1;
        }
        else if (vec.z != 0)
        {
            vec.x = 1;
        }

        return vec;
    }

    QuadForm<T> get_neighbor(R t, Isometry<T>& s) const
    {
        Vector3<R> vec = this->isotropic_vector(t);
        return build_neighbor(vec, s);
    }

    QuadForm<T> build_neighbor(Vector3<R>& vec2, Isometry<T>& s) const
    {
        T p = GF->prime();
        T pp = p*p;
        T aa, bb, cc, ff, gg, hh;

        // Convert isotropic vector into the correct domain.
        Vector3<T> vec;
        vec.x = GF->mod(vec2.x);
        vec.y = GF->mod(vec2.y);
        vec.z = GF->mod(vec2.z);

        #ifdef DEBUG
        assert( q.evaluate(vec) % p == 0 );
        #endif

        if (vec.x > (p>>1)) vec.x -= p;
        if (vec.y > (p>>1)) vec.y -= p;

        if (vec.z == 1)
        {
            T u = vec.x;
            T v = vec.y;

            s.set_values(u, 0, -1, v, 1, 0, 1, 0, 0);

            T au = q.a() * u;
            T hu = q.h() * u;
            T hv = q.h() * v;
            T bv = q.b() * v;
            T temp = au + hv + q.g();

            aa = q.c() + u*temp + v * (bv + q.f());
            bb = q.b();
            cc = q.a();
            ff = -q.h();
            gg = -(temp + au);
            hh = q.f() + hu + 2*bv;
        }
        else if (vec.y == 1)
        {
            T u = vec.x;

            s.set_values(u, 0, 1, 1, 0, 0, 0, 1, 0);
            T au = q.a()*u;
            T gu = q.g()*u;
            T temp = au + q.h();

            aa = q.b() + u*temp;
            bb = q.c();
            cc = q.a();
            ff = q.g();
            gg = temp + au;
            hh = q.f() + gu;
        }
        else
        {
            s.set_values(1, 0, 0, 0, 0, -1, 0, 1, 0);

            aa = q.a();
            bb = q.c();
            cc = q.b();
            ff = -q.f();
            gg = -q.h();
            hh = q.g();
        }

        #ifdef DEBUG
        QuadForm<T> test(aa, bb, cc, ff, gg, hh);
        assert( test.discriminant() == q.discriminant() );
        assert( aa % p == 0 );
        #endif

        T gmodp = gg % p;
        if (gmodp == 0)
        {
            s.A1000010n0();

            T temp = bb;
            bb = cc;
            cc = temp;

            temp = gg;
            gg = hh;
            hh = -temp;
            ff = -ff;

            gmodp = gg % p;
        }

        #ifdef DEBUG
        assert( gg % p != 0 );
        assert( gg % p == gmodp );
        #endif

        if (gmodp < 0) gmodp += p;
        T ginv = GF->inverse(gmodp)%p;

        #ifdef DEBUG
        assert( (gg * ginv) % p == 1 || (gg * ginv) % p == 1-p );
        #endif

        T s1 = (-hh * ginv) % p;
        T s2 = (-aa * ginv) % pp;
        if (s1 < -(p>>1)) s1 += p;
        if (s2 < -(pp>>1)) s2 += pp;

        s.A1000100t1(s1);

        T temp1 = cc * s1;
        bb += s1 * (ff + temp1);
        ff += 2*temp1;
        hh += s1 * gg;

        #ifdef DEBUG
        assert( hh % p == 0 );
        #endif

        s.A100010t01(s2);
        T temp2 = cc * s2;
        aa += s2 * (gg + temp2);
        gg += 2 * temp2;
        hh += s2 * ff;

        #ifdef DEBUG
        assert( aa > 0 );
        assert( aa % pp == 0 );
        assert( hh % p == 0 );
        #endif

        s.A1000p000p2(p, pp);
        aa /= pp;
        cc *= pp;
        ff *= p;
        hh /= p;

        QuadForm<T> retval(aa, bb, cc, ff, gg, hh);
        if (std::is_same<T,Z64>::value)
        {
            // If we're not using arbitrary precision, throw an exception if
            // the discriminant of the p-neighbor isn't correct.
            if (retval.discriminant() != this->disc)
            {
                throw std::overflow_error(
                    "An overflow has occurred. The p-neighbor's discriminant "
                    "does not match the original.");
            }
        }
        return retval;
    }

private:
    std::shared_ptr<Fp<R,S>> GF;
    QuadForm<T> q;
    T disc;
    R a, b, c, f, g, h;
    Vector3<R> vec;
    R a0;
    R delta;

    // The 2-isotropic vectors were stored in binary within each of the
    // coordinates of `vec` and so we use this function to unpack them into
    // actual 2-isotropic vectors.
    Vector3<R> isotropic_vector_p2(R t) const
    {
        Vector3<R> res;

        if (t == 0)
        {
            res.z = this->vec.x & 1;
            res.y = !!(this->vec.x & 2);
            res.x = !!(this->vec.x & 4);
        }
        else if (t == 1)
        {
            res.z = this->vec.y & 1;
            res.y = !!(this->vec.y & 2);
            res.x = !!(this->vec.y & 4);
        }
        else
        {
            res.z = this->vec.z & 1;
            res.y = !!(this->vec.z & 2);
            res.x = !!(this->vec.z & 4);
        }

        return res;
    }
};

#endif // __NEIGHBOR_MANAGER_H
