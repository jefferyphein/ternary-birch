#include <iostream>
#include <cstdlib>
#include <gmpxx.h>
#include "QuadForm.h"
#include "Genus.h"
#include "Character.h"

typedef Isometry<mpz_class, mpq_class> IsometryZZ;
typedef Genus<mpz_class, mpq_class> GenusZZ;
typedef QuadForm<mpz_class, mpq_class> QuadFormZZ;
typedef Character<mpz_class, mpq_class> CharacterZZ;

int main(int argc, char** argv)
{
    mpz_class a, b, c, f, g, h;

    if (argc != 7)
    {
        // discriminant = 85085
        //mpz_class a(1);
        //mpz_class b(130);
        //mpz_class c(172);
        //mpz_class f(65);
        //mpz_class g(1);
        //mpz_class h(0);

        // discriminant = 1062347
        a = 1;
        b = 187;
        c = 1467;
        f = -187;
        g = 0;
        h = 0;
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
    std::cout << q.discriminant() << std::endl;

    GenusZZ genus(q);
    std::vector<mpz_class> vec = {11, 13};
    CharacterZZ repr(vec);
    genus.add_character(repr);
    genus.print();

    return EXIT_SUCCESS;
}
