#ifndef __AUTOMORPHISM_ZZ_H
#define __AUTOMORPHISM_ZZ_H

#include <memory>
#include "Isometry.h"

// This class defines a number of integral automorphisms with determinant +1
//  for convenient reference throughout various other classes and objects.
//  The goal being that we wish to avoid consuming too much memory if we need
//  to access automorphisms across various forms with the same automorphism
//  matrices.
//
// The naming convention is row major, so that we have:
//
//  Aabcdefghi = | a  b  c |
//               | d  e  f |
//               | g  h  i |
//
//  where 1 -->  1
//        0 -->  0
//        N --> -1
class AutomorphismZZ {
public:
    typedef Isometry<mpz_class, mpq_class> IsometryQQ;
    typedef std::shared_ptr<IsometryQQ> IsometryQQPtr;

    static IsometryQQPtr A100010001;
    static IsometryQQPtr ANNN001010;
    static IsometryQQPtr ANN001000N;
    static IsometryQQPtr AN0N0N0001;
    static IsometryQQPtr AN0000N0N0;
    static IsometryQQPtr A10100N010;
    static IsometryQQPtr A1100010N0;
    static IsometryQQPtr A1110N000N;
    static IsometryQQPtr AN000NN001;
    static IsometryQQPtr AN00N1000N;
    static IsometryQQPtr AN000N0001;
    static IsometryQQPtr AN10N00001;
    static IsometryQQPtr AN1001000N;
    static IsometryQQPtr A0N0N0000N;
    static IsometryQQPtr A0N01N0001;
    static IsometryQQPtr A010N10001;
    static IsometryQQPtr A01010000N;
    static IsometryQQPtr A1N00N000N;
    static IsometryQQPtr A1N0100001;
    static IsometryQQPtr A1001N000N;
    static IsometryQQPtr AN0001N00N;
    static IsometryQQPtr AN010N1001;
    static IsometryQQPtr A0N1100001;
    static IsometryQQPtr A01N10N00N;
    static IsometryQQPtr A010N01001;
    static IsometryQQPtr A10N0N000N;
    static IsometryQQPtr AN0001000N;
    static IsometryQQPtr AN010N0001;
    static IsometryQQPtr AN000N00N1;
    static IsometryQQPtr AN000N1001;
    static IsometryQQPtr AN00001010;
    static IsometryQQPtr AN0001001N;
    static IsometryQQPtr A1000N000N;
    static IsometryQQPtr A1000N10N0;
    static IsometryQQPtr A10000N01N;
    static IsometryQQPtr A1000010N1;
    static IsometryQQPtr A10001N010;
    static IsometryQQPtr AN00N01N10;
    static IsometryQQPtr AN01N10N00;
    static IsometryQQPtr AN10N00N01;
    static IsometryQQPtr AN1001001N;
    static IsometryQQPtr A0N01N00N1;
    static IsometryQQPtr A0N10N01N0;
    static IsometryQQPtr A0N1001N01;
    static IsometryQQPtr A00N0N0N00;
    static IsometryQQPtr A00N01N10N;
    static IsometryQQPtr A001N010N1;
    static IsometryQQPtr A001100010;
    static IsometryQQPtr A01NN10010;
    static IsometryQQPtr A010001100;
    static IsometryQQPtr A01001NN10;
    static IsometryQQPtr A1N00N10N0;
    static IsometryQQPtr A1N010N100;
    static IsometryQQPtr A10N00N01N;
    static IsometryQQPtr A10N1001N0;
    static IsometryQQPtr A1001N010N;
    static IsometryQQPtr A0N000N100;
    static IsometryQQPtr A0N0001N00;
    static IsometryQQPtr A0N0100001;
    static IsometryQQPtr A00NN00010;
    static IsometryQQPtr A00N010100;
    static IsometryQQPtr A00N1000N0;
    static IsometryQQPtr A001N000N0;
    static IsometryQQPtr A0010N0100;
    static IsometryQQPtr A001010N00;
    static IsometryQQPtr A010N00001;
    static IsometryQQPtr A01000NN00;
    static IsometryQQPtr A10000N010;
    static IsometryQQPtr A1000010N0;
    static IsometryQQPtr ANNN010100;
    static IsometryQQPtr ANNN100001;
    static IsometryQQPtr AN000N0111;
    static IsometryQQPtr AN0011100N;
    static IsometryQQPtr A0N000N111;
    static IsometryQQPtr A0N0111N00;
    static IsometryQQPtr A00NN00111;
    static IsometryQQPtr A00N1110N0;
    static IsometryQQPtr A001NNN100;
    static IsometryQQPtr A001010NNN;
    static IsometryQQPtr A010NNN001;
    static IsometryQQPtr A010100NNN;
    static IsometryQQPtr A100NNN010;
    static IsometryQQPtr A100001NNN;
    static IsometryQQPtr A111N000N0;
    static IsometryQQPtr A11100NN00;
};

#endif /* AUTOMORPHISM_ZZ_H */
