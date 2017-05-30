#include <iostream>
#include <cstdlib>
#include <gmpxx.h>
#include "QuadForm.h"
#include "Genus.h"

typedef Isometry<mpz_class, mpq_class> IsometryZZ;
typedef Genus<mpz_class, mpq_class> GenusZZ;
typedef QuadForm<mpz_class, mpq_class> QuadFormZZ;

int main(int, char**)
{
    // discriminant = 85085
    mpz_class a(1);
    mpz_class b(130);
    mpz_class c(172);
    mpz_class f(65);
    mpz_class g(1);
    mpz_class h(0);

    // discriminant = 1062347
    //mpz_class a(1);
    //mpz_class b(187);
    //mpz_class c(1467);
    //mpz_class f(-187);
    //mpz_class g(0);
    //mpz_class h(0);
    QuadFormZZ q(a, b, c, f, g, h);

    //std::cout << q.discriminant() << std::endl;
    //std::cout << q << std::endl;

    GenusZZ genus(q);
    
    return EXIT_SUCCESS;
}
