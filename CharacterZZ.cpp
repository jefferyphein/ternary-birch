#include "Character.h"
#include "Isometry.h"
#include "QuadForm.h"
#include "Math.h"

typedef Character<mpz_class, mpq_class> CharacterZZ;
typedef Isometry<mpz_class, mpq_class> IsometryQQ;
typedef QuadForm<mpz_class, mpq_class> QuadFormZZ;
typedef Math<mpz_class, mpq_class> MathZZ;

template<>
int64_t CharacterZZ::rho(const IsometryQQ& s, const QuadFormZZ& q) const
{
    int64_t value = 1;

    if (this->ps_.size() != 0)
    {
        mpq_class spin = s.spinor_norm(q);
        for (auto& p : this->ps_)
        {
            int64_t v = abs(MathZZ::valuation(spin, p) % 2);
            value *= (1 - 2 * v);
        }
    }

    return value;
}
