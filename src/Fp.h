#ifndef __FP_H_
#define __FP_H_

#include "birch.h"

template<typename R, typename S>
class Fp
{
public:
    Fp(const R& p, W64 seed, bool use_inverse_lut=false)
    {
        std::random_device rd;
        this->rng = std::unique_ptr<std::mt19937>( new std::mt19937(seed) );
        this->distr = std::unique_ptr<std::uniform_int_distribution<>>(
            new std::uniform_int_distribution<>(0, p-1));

        this->p = p;
        if (this->p != 2)
        {
            this->kp = (((R)-1)/p)*p;
            this->kp_inv = ((S)-1)/kp;
            this->use_inverse_lut = use_inverse_lut;
            if (use_inverse_lut) this->make_inverse_lut();
        }
    }

    const R& prime(void) { return this->p; }

    template<typename T>
    inline R mod(const T& a) const
    {
        static_assert(std::is_integral<T>::value, "Undefined type.");
        T value = (T)a % this->p;
        return (value < 0) ? (R)(value+this->p) : (R)value;
    }

    template<typename T>
    inline Vector3<R> mod(const Vector3<T>& vec) const
    {
        Vector3<R> res;
        res.x = this->mod(vec.x);
        res.y = this->mod(vec.y);
        res.z = this->mod(vec.z);
        return res;
    }

    inline virtual R neg(R a) const
    {
        return this->kp-a;
    }

    inline virtual R mul(R a, R b) const
    {
        S rem = ((S)a)*b;
        R hi = rem >> bits;
        R t = (((S)hi*(S)this->kp_inv) >> bits) + hi;
        rem -= (S)t * this->kp;
        rem = (rem >= this->kp) ? rem-this->kp : rem;
        rem = (rem >= this->kp) ? rem-this->kp : rem;
        rem = (rem >= this->kp) ? rem-this->kp : rem;

        #ifdef DEBUG
        assert( ((S)a*(S)b)%p == rem%p );
        #endif

        return rem;
    }

    inline virtual R add(R a, R b) const
    {
        R neg = this->kp-a;
        return (b >= neg) ? b-neg : this->kp-(neg-b);
    }

    inline virtual R sub(R a, R b) const
    {
        return add(a, this->kp-b);
    }

    inline virtual R pow(R a, Z64 e) const
    {
        if (e == 0) return 1;
        if (e == 1) return a;

        R temp = this->pow(a, e>>1);
        temp = this->mul(temp, temp);
        return (e%2) ? this->mul(temp, a) : temp;
    }

    inline int legendre(R a) const
    {
        Z aa(a);
        Z pp(p);
        return mpz_legendre(aa.get_mpz_t(), pp.get_mpz_t());
    }

    inline virtual R sqrt(R a) const
    {
        a = a % p;
        if (a == 1) return 1;
        if (a == 0) return 0;
        if (this->legendre(a) != 1) return 0;

        R q = p-1;
        R s = 0;
        while(q % 2 == 0) { q >>= 1; s++; }

        if (s == 1) return this->pow(a, (p+1)/4);

        R z = this->random();
        while (z == 0 || this->legendre(z) == 1)
        {
            z = this->random();
        }

        int m = s;
        R c = this->pow(z, q);
        R r = this->pow(a, (q+1)/2);
        R t = this->pow(a, q);
        if (t >= p) t = this->mod(t);

        while (1)
        {
            if (t == 1) return r;
            int i = 0;
            R t1 = t;
            while (t1 != 1)
            {
                t1 = this->mul(t1, t1);
                if (t1 >= p) t1 = this->mod(t1);
                i++;
            }

            int e = 1;
            for (int j=0; j<m-i-1; j++) e <<= 1;
            R b = this->pow(c, e);
            r = this->mul(r, b);
            c = this->mul(b, b);
            t = this->mul(t, c);
            if (t >= p) t %= p;
            m = i;
        }

        return 0;
    }

    inline virtual R inverse(R a) const
    {
        if (this->use_inverse_lut) return this->inverse_lut[a];
        else return this->inv(a);
    }

    inline virtual R inverse(const Z& a) const
    {
        R inv = mpz_get_ui(a.get_mpz_t());
        return this->inverse(inv);
    }

    inline virtual R inverse(const Z64& a) const
    {
        R inv = (R)a;
        return this->inverse(inv);
    }

    inline R random(void) const
    {
        return (R)(*this->distr)(*this->rng);
    }

private:
    R p;
    R kp;
    R kp_inv;
    static constexpr int bits = 8 * sizeof(R);
    bool use_inverse_lut;
    std::vector<R> inverse_lut;

    // Random number generator.
    std::unique_ptr<std::mt19937> rng;
    std::unique_ptr<std::uniform_int_distribution<>> distr;

    inline virtual R inv(R a) const
    {
        if (a == 0) return 0;
        Z aa(a);
        Z pp(p);
        mpz_invert(aa.get_mpz_t(), aa.get_mpz_t(), pp.get_mpz_t());
        R ainv = mpz_get_ui(aa.get_mpz_t());

        #ifdef DEBUG
        assert( ((S)a * (S)ainv) % p == 1 );
        #endif

        return ainv;
    }

    void inverse_lut_populate(Z32 offset, Z32 len)
    {
        Z32 datalen = len<<1;

        // Set the initial values.
        std::iota(this->inverse_lut.begin()+offset,
            this->inverse_lut.begin()+offset+len, offset);

        // Multipley up to the root node...
        for (Z32 i=0, j=len; i<j; i+=2, j++)
        {
            R a = this->inverse_lut[offset+i];
            R b = this->inverse_lut[offset+i+1];
            this->inverse_lut[offset+j] = this->mul(a, b);
        }

        // ...invert the root node...
        R ainv = this->inv(this->inverse_lut[offset+datalen-2]);
        this->inverse_lut[offset+datalen-2] = ainv;

        // ...and then backtrack to the inverse.
        for (Z32 i=datalen-4, j=datalen-2; i>=0; i-=2, --j)
        {
            R temp = this->inverse_lut[offset+i];
            R a = this->inverse_lut[offset+i+1];
            R b = this->inverse_lut[offset+j];
            this->inverse_lut[offset+i] = this->mul(a, b);
            this->inverse_lut[offset+i+1] = this->mul(temp, b);
        }
    }

    void make_inverse_lut(void)
    {
        // TODO: The following could probably be replaced with clz.
        Z32 len = 1;
        R q = p;
        while (q != 1)
        {
            len <<= 1;
            q >>= 1;
        }

        // Make enough room in the lut to grow the tree.
        inverse_lut.resize(1+(len<<1));

        Z32 offset = 1;
        q = p;
        while (q > 1)
        {
            this->inverse_lut_populate(offset, len);

            offset |= len;
            q ^= len;

            R r = q;
            len = 1;
            while (r > 1)
            {
                len <<= 1;
                r >>= 1;
            }
        }

        // Shrink the lut to the appropriate size.
        this->inverse_lut.resize(p);

        #ifdef DEBUG
        for (Z32 i=1; i<p; i++)
        {
            assert( this->mul(i, this->inverse_lut[i]) % p == 1 );
        }
        #endif
    }
};

template<typename R, typename S>
class F2 : public Fp<R,S>
{
public:
    F2(const R& p, W64 seed) : Fp<R,S>(p, seed, false) {}

    inline R mul(R a, R b) const override
    {
        return ((a & b) & 1);
    }

    inline R add(R a, R b) const override
    {
        return ((a ^ b) & 1);
    }

    inline R sub(R a, R b) const override
    {
        return ((a ^ b) & 1);
    }

    inline R pow(R a, Z64 e) const override
    {
        return e == 0 ? 1 : (a & 1);
    }

    inline R sqrt(R a) const override
    {
        return (a & 1);
    }

    inline R inverse(R a) const override
    {
        return (a & 1);
    }

    inline R inverse(const Z& a) const override
    {
        return (mpz_get_ui(a.get_mpz_t()) & 1);
    }

    inline R inverse(const Z64& a) const override
    {
        return (a & 1);
    }
private:
    inline R inv(R a) const override
    {
        return (a & 1);
    }
};

template<>
template<>
W16 W16_Fp::mod(const Z& a) const;

template<>
template<>
W32 W32_Fp::mod(const Z& a) const;

template<>
template<>
W64 W64_Fp::mod(const Z& a) const;

#endif // __FP_H_
