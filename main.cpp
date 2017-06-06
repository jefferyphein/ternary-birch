#include <iostream>
#include <cstdlib>
#include <gmpxx.h>
#include "QuadForm.h"
#include "Genus.h"
#include "Character.h"
#include "Math.h"

typedef Isometry<mpz_class, mpq_class> IsometryZZ;
typedef Genus<mpz_class, mpq_class> GenusZZ;
typedef QuadForm<mpz_class, mpq_class> QuadFormZZ;
typedef Character<mpz_class, mpq_class> CharacterZZ;

int main(int argc, char** argv)
{
    mpz_class a, b, c, f, g, h;

    if (argc != 7)
    {
        std::cerr << "Usage:" << std::endl;
        std::cerr << "   ./birch a b c f g h" << std::endl;
    }
    else
    {
        a = atoi(argv[1]);
        b = atoi(argv[2]);
        c = atoi(argv[3]);
        f = atoi(argv[4]);
        g = atoi(argv[5]);
        h = atoi(argv[6]);
    }

    QuadFormZZ q(a, b, c, f, g, h);
    GenusZZ genus(q);

    auto divs = Math::squarefree_divisors(q.discriminant());
    for (auto& d : divs)
    {
        CharacterZZ chi(d);
        genus.add_character(chi);
    }

    genus.compute_genus(0);

    genus.compute_hecke_operators(7, 0);

    genus.print();
    
    return EXIT_SUCCESS;
}
