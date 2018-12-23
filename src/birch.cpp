#include "birch.h"
#include "Genus.h"

int main(int argc, char **argv)
{
    std::vector<Z_PrimeSymbol> symbols;
    Z_PrimeSymbol p;

    p.p = 11;
    p.power = 1;
    p.ramified = true;
    symbols.push_back(p);

    p.p = 13;
    p.power = 1;
    p.ramified = true;
    symbols.push_back(p);

    p.p = 17;
    p.power = 1;
    p.ramified = true;
    symbols.push_back(p);

    p.p = 19;
    p.power = 1;
    p.ramified = true;
    symbols.push_back(p);

    p.p = 23;
    p.power = 1;
    p.ramified = true;
    symbols.push_back(p);

    Z_Isometry s;
    Z_QuadForm q = Z_QuadForm::get_quad_form(symbols);

    Z_Genus genus1(q, symbols);
    Z64_Genus genus2 = Z_Genus::convert<Z64>(genus1);

    //Z prime;
    //prime = 2;
    //genus1.hecke_matrix_sparse(prime);
    //prime = 3;
    //genus1.hecke_matrix_sparse(prime);
    //prime = 5;
    //genus1.hecke_matrix_sparse(prime);
    //prime = 7;
    //genus1.hecke_matrix_sparse(prime);
    //prime = 31;
    //genus1.hecke_matrix_sparse(prime);
    //prime = 101;
    //genus1.hecke_matrix_sparse(prime);
    //prime = 373;
    //genus1.hecke_matrix_sparse(prime);
    //prime = 8191;
    //genus1.hecke_matrix_dense(prime);


    //std::vector<Z64_PrimeSymbol> symbols64(symbols.size());
    //std::transform(symbols.begin(), symbols.end(), symbols64.begin(), birch_util::convert_PrimeSymbol<Z,Z64>);

    //Z64_QuadForm q64 = birch_util::convert_QuadForm<Z,Z64>(q);
    //Z64_Genus genus2(q64, symbols64);

    //genus2.hecke_matrix_sparse(2);
    //genus2.hecke_matrix_sparse(3);
    //genus2.hecke_matrix_sparse(5);
    //genus2.hecke_matrix_sparse(7);
    //genus2.hecke_matrix_sparse(101);
    //genus2.hecke_matrix_sparse(373);
    //genus2.hecke_matrix_dense(373);
    genus2.hecke_matrix_dense(8191);

    return EXIT_SUCCESS;
}
