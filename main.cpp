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
    std::string genfilename;
    std::string eigfilename;

    if (argc < 9)
    {
        std::cerr << "Usage:" << std::endl;
        std::cerr << "   ./birch a b c f g h genus eigs" << std::endl;
    }
    else
    {
        a = atoi(argv[1]);
        b = atoi(argv[2]);
        c = atoi(argv[3]);
        f = atoi(argv[4]);
        g = atoi(argv[5]);
        h = atoi(argv[6]);

        genfilename = argv[7];
        eigfilename = argv[8];
    }

    QuadFormZZ q(a, b, c, f, g, h);
    GenusZZ genus(q);

    auto divs = Math::squarefree_divisors(q.discriminant());
    for (auto& d : divs)
    {
        CharacterZZ chi(d);
        genus.add_character(chi);
    }

    //genus.compute_genus(0);
    genus.import_genus(genfilename);

    genus.import_eigenvectors(eigfilename);

    std::vector<mpz_class> ps = {
          2,   3,   5,   7,  29,  31,  37,  41,  43,  47,
         53,  59,  61,  67,  71,
         73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
        127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
        179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
        233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
        283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
        353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
        419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
        467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
        547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
        607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
        661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
        739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
        811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
        877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
        947, 953, 967, 971, 977, 983, 991, 997,
1009 ,   1013,
1019 ,   1021,    1031,    1033,    1039,    1049,    1051,    1061,    1063,    1069,
1087 ,   1091,    1093,    1097,    1103,    1109,    1117,    1123,    1129,    1151,
1153 ,   1163,    1171,    1181,    1187,    1193,    1201,    1213,    1217,    1223 };

    genus.compute_eigenvalues(ps, 8);

    //for (const mpz_class& p : ps)
    //{
    //    genus.compute_eigenvalues(p, 4);
    //}

    //genus.compute_hecke_operators(2, 8);
    //genus.compute_hecke_operators(3, 8);
    //genus.compute_hecke_operators(5, 8);
    //genus.compute_hecke_operators(7, 8);

    //genus.print();
    
    return EXIT_SUCCESS;
}
