#include "Isometry.h"

typedef QuadForm<mpz_class, mpq_class> QuadFormZZ;
typedef Isometry<mpz_class, mpq_class> IsometryQQ;

template<>
mpq_class IsometryQQ::spinor_norm(const QuadFormZZ& q) const
{
    // Compute the determinant of the isometry.
    mpq_class det = this->a11 * (this->a22 * this->a33 - this->a23 * this->a32) -
                    this->a12 * (this->a21 * this->a33 - this->a23 * this->a31) +
                    this->a13 * (this->a21 * this->a32 - this->a22 * this->a31);

    // We must compute with isometries having determinant +1, and so we
    // multiply the trace by a factor of det to account for this.
    mpq_class tr = det * this->trace();

    // When the trace is not equal to -1, then we have a very simple formula
    // for the spinor norm.
    if (tr != -1)
    {
        return tr + 1;
    }

    // All of these computations are derived from computations with quaternion
    // algebras. This isometry corresponds to the action of a quaternion on the
    // trace zero subspace associated to a ternary quadratic space. This
    // quaternion is denoted:
    //      alpha = t + x*i + y*j + z*k
    // for appropriately chosen i, j, k, x, y, z. In this particular case, we
    // have t = 0.

    // Handle the case where x != 0, as above. We need to multiply all the
    // elements of the isometry by det to account for the fact that we must have
    // determinant +1.
    mpq_class m22m33 = det * (this->a22 + this->a33 -
                             (q.g() * this->a31 + q.h() * this->a21) / (2 * q.a()));
    if (m22m33 != 0)
    {
        return -q.discriminant() / (2 * q.a() * m22m33);
    }

    // Handle the case where x == 0 and y != 0, as above.
    mpq_class abh = 4 * q.a() * q.b() - q.h() * q.h();
    mpq_class m33 = det * (this->a33 + ((q.g() * q.h() - 2 * q.a() * q.f()) *
                           this->a32 + (q.f() * q.h() - 2 * q.b() * q.g()) *
                           this->a31) / abh);
    if (m33 != 1)
    {
        return -2 * q.a() * q.discriminant() * abh / (m33 - 1);
    }

    // Handle the case where x == 0 and y == 0 and z != 0, as above.
    return abh;
}
