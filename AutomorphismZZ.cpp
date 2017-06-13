#include <memory>
#include "AutomorphismZZ.h"

typedef Isometry<mpz_class, mpq_class> IsometryQQ;
typedef std::shared_ptr<IsometryQQ> IsometryQQPtr;

// Assign all automorphisms used by the QuadForm::automorphisms method.
//  The naming mechanism here is such that each automorphism is assigned a
//  name starting with "A" and followed by nine characters given by 1, 0, or N
//  assigned according to the row major enumeration of the coefficients of the
//  automorphism. In this case, an "N" corresponds to the coefficient -1, where
//  "0" and "1" are assigned consistently with the value of the coefficient.
// This is not a complete list of automorphisms, as each automorphism listed
//  here has determinant +1. This is intentional, since automorphism groups on
//  O(3) always contain -I where I is the identity matrix. As such, when
//  computing the automorphism group of a quadratic form, not only does each
//  corresponding automorphism listed here belong to the group, but so does its
//  negative.
IsometryQQPtr AutomorphismZZ::A100010001 = std::make_shared<IsometryQQ>(1, 0, 0, 0, 1, 0, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::ANNN001010 = std::make_shared<IsometryQQ>(-1, -1, -1, 0, 0, 1, 0, 1, 0);
IsometryQQPtr AutomorphismZZ::ANN001000N = std::make_shared<IsometryQQ>(-1, -1, 0, 0, 1, 0, 0, 0, -1);
IsometryQQPtr AutomorphismZZ::AN0N0N0001 = std::make_shared<IsometryQQ>(-1, 0, -1, 0, -1, 0, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::AN0000N0N0 = std::make_shared<IsometryQQ>(-1, 0, 0, 0, 0, -1, 0, -1, 0);
IsometryQQPtr AutomorphismZZ::A10100N010 = std::make_shared<IsometryQQ>(1, 0, 1, 0, 0, -1, 0, 1, 0);
IsometryQQPtr AutomorphismZZ::A1100010N0 = std::make_shared<IsometryQQ>(1, 1, 0, 0, 0, 1, 0, -1, 0);
IsometryQQPtr AutomorphismZZ::A1110N000N = std::make_shared<IsometryQQ>(1, 1, 1, 0, -1, 0, 0, 0, -1);
IsometryQQPtr AutomorphismZZ::AN000NN001 = std::make_shared<IsometryQQ>(-1, 0, 0, 0, -1, -1, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::AN00N1000N = std::make_shared<IsometryQQ>(-1, 0, 0, -1, 1, 0, 0, 0, -1);
IsometryQQPtr AutomorphismZZ::AN000N0001 = std::make_shared<IsometryQQ>(-1, 0, 0, 0, -1, 0, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::AN10N00001 = std::make_shared<IsometryQQ>(-1, 1, 0, -1, 0, 0, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::AN1001000N = std::make_shared<IsometryQQ>(-1, 1, 0, 0, 1, 0, 0, 0, -1);
IsometryQQPtr AutomorphismZZ::A0N0N0000N = std::make_shared<IsometryQQ>(0, -1, 0, -1, 0, 0, 0, 0, -1);
IsometryQQPtr AutomorphismZZ::A0N01N0001 = std::make_shared<IsometryQQ>(0, -1, 0, 1, -1, 0, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::A010N10001 = std::make_shared<IsometryQQ>(0, 1, 0, -1, 1, 0, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::A01010000N = std::make_shared<IsometryQQ>(0, 1, 0, 1, 0, 0, 0, 0, -1);
IsometryQQPtr AutomorphismZZ::A1N00N000N = std::make_shared<IsometryQQ>(1, -1, 0, 0, -1, 0, 0, 0, -1);
IsometryQQPtr AutomorphismZZ::A1N0100001 = std::make_shared<IsometryQQ>(1, -1, 0, 1, 0, 0, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::A1001N000N = std::make_shared<IsometryQQ>(1, 0, 0, 1, -1, 0, 0, 0, -1);
IsometryQQPtr AutomorphismZZ::AN0001N00N = std::make_shared<IsometryQQ>(-1, 0, 0, 0, 1, -1, 0, 0, -1);
IsometryQQPtr AutomorphismZZ::AN010N1001 = std::make_shared<IsometryQQ>(-1, 0, 1, 0, -1, 1, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::A0N1100001 = std::make_shared<IsometryQQ>(0, -1, 1, 1, 0, 0, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::A01N10N00N = std::make_shared<IsometryQQ>(0, 1, -1, 1, 0, -1, 0, 0, -1);
IsometryQQPtr AutomorphismZZ::A010N01001 = std::make_shared<IsometryQQ>(0, 1, 0, -1, 0, 1, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::A10N0N000N = std::make_shared<IsometryQQ>(1, 0, -1, 0, -1, 0, 0, 0, -1);
IsometryQQPtr AutomorphismZZ::AN0001000N = std::make_shared<IsometryQQ>(-1, 0, 0, 0, 1, 0, 0, 0, -1);
IsometryQQPtr AutomorphismZZ::AN010N0001 = std::make_shared<IsometryQQ>(-1, 0, 1, 0, -1, 0, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::AN000N00N1 = std::make_shared<IsometryQQ>(-1, 0, 0, 0, -1, 0, 0, -1, 1);
IsometryQQPtr AutomorphismZZ::AN000N1001 = std::make_shared<IsometryQQ>(-1, 0, 0, 0, -1, 1, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::AN00001010 = std::make_shared<IsometryQQ>(-1, 0, 0, 0, 0, 1, 0, 1, 0);
IsometryQQPtr AutomorphismZZ::AN0001001N = std::make_shared<IsometryQQ>(-1, 0, 0, 0, 1, 0, 0, 1, -1);
IsometryQQPtr AutomorphismZZ::A1000N000N = std::make_shared<IsometryQQ>(1, 0, 0, 0, -1, 0, 0, 0, -1);
IsometryQQPtr AutomorphismZZ::A1000N10N0 = std::make_shared<IsometryQQ>(1, 0, 0, 0, -1, 1, 0, -1, 0);
IsometryQQPtr AutomorphismZZ::A10000N01N = std::make_shared<IsometryQQ>(1, 0, 0, 0, 0, -1, 0, 1, -1);
IsometryQQPtr AutomorphismZZ::A1000010N1 = std::make_shared<IsometryQQ>(1, 0, 0, 0, 0, 1, 0, -1, 1);
IsometryQQPtr AutomorphismZZ::A10001N010 = std::make_shared<IsometryQQ>(1, 0, 0, 0, 1, -1, 0, 1, 0);
IsometryQQPtr AutomorphismZZ::AN00N01N10 = std::make_shared<IsometryQQ>(-1, 0, 0, -1, 0, 1, -1, 1, 0);
IsometryQQPtr AutomorphismZZ::AN01N10N00 = std::make_shared<IsometryQQ>(-1, 0, 1, -1, 1, 0, -1, 0, 0);
IsometryQQPtr AutomorphismZZ::AN10N00N01 = std::make_shared<IsometryQQ>(-1, 1, 0, -1, 0, 0, -1, 0, 1);
IsometryQQPtr AutomorphismZZ::AN1001001N = std::make_shared<IsometryQQ>(-1, 1, 0, 0, 1, 0, 0, 1, -1);
IsometryQQPtr AutomorphismZZ::A0N01N00N1 = std::make_shared<IsometryQQ>(0, -1, 0, 1, -1, 0, 0, -1, 1);
IsometryQQPtr AutomorphismZZ::A0N10N01N0 = std::make_shared<IsometryQQ>(0, -1, 1, 0, -1, 0, 1, -1, 0);
IsometryQQPtr AutomorphismZZ::A0N1001N01 = std::make_shared<IsometryQQ>(0, -1, 1, 0, 0, 1, -1, 0, 1);
IsometryQQPtr AutomorphismZZ::A00N0N0N00 = std::make_shared<IsometryQQ>(0, 0, -1, 0, -1, 0, -1, 0, 0);
IsometryQQPtr AutomorphismZZ::A00N01N10N = std::make_shared<IsometryQQ>(0, 0, -1, 0, 1, -1, 1, 0, -1);
IsometryQQPtr AutomorphismZZ::A001N010N1 = std::make_shared<IsometryQQ>(0, 0, 1, -1, 0, 1, 0, -1, 1);
IsometryQQPtr AutomorphismZZ::A001100010 = std::make_shared<IsometryQQ>(0, 0, 1, 1, 0, 0, 0, 1, 0);
IsometryQQPtr AutomorphismZZ::A01NN10010 = std::make_shared<IsometryQQ>(0, 1, -1, -1, 1, 0, 0, 1, 0);
IsometryQQPtr AutomorphismZZ::A010001100 = std::make_shared<IsometryQQ>(0, 1, 0, 0, 0, 1, 1, 0, 0);
IsometryQQPtr AutomorphismZZ::A01001NN10 = std::make_shared<IsometryQQ>(0, 1, 0, 0, 1, -1, -1, 1, 0);
IsometryQQPtr AutomorphismZZ::A1N00N10N0 = std::make_shared<IsometryQQ>(1, -1, 0, 0, -1, 1, 0, -1, 0);
IsometryQQPtr AutomorphismZZ::A1N010N100 = std::make_shared<IsometryQQ>(1, -1, 0, 1, 0, -1, 1, 0, 0);
IsometryQQPtr AutomorphismZZ::A10N00N01N = std::make_shared<IsometryQQ>(1, 0, -1, 0, 0, -1, 0, 1, -1);
IsometryQQPtr AutomorphismZZ::A10N1001N0 = std::make_shared<IsometryQQ>(1, 0, -1, 1, 0, 0, 1, -1, 0);
IsometryQQPtr AutomorphismZZ::A1001N010N = std::make_shared<IsometryQQ>(1, 0, 0, 1, -1, 0, 1, 0, -1);
IsometryQQPtr AutomorphismZZ::A0N000N100 = std::make_shared<IsometryQQ>(0, -1, 0, 0, 0, -1, 1, 0, 0);
IsometryQQPtr AutomorphismZZ::A0N0001N00 = std::make_shared<IsometryQQ>(0, -1, 0, 0, 0, 1, -1, 0, 0);
IsometryQQPtr AutomorphismZZ::A0N0100001 = std::make_shared<IsometryQQ>(0, -1, 0, 1, 0, 0, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::A00NN00010 = std::make_shared<IsometryQQ>(0, 0, -1, -1, 0, 0, 0, 1, 0);
IsometryQQPtr AutomorphismZZ::A00N010100 = std::make_shared<IsometryQQ>(0, 0, -1, 0, 1, 0, 1, 0, 0);
IsometryQQPtr AutomorphismZZ::A00N1000N0 = std::make_shared<IsometryQQ>(0, 0, -1, 1, 0, 0, 0, -1, 0);
IsometryQQPtr AutomorphismZZ::A001N000N0 = std::make_shared<IsometryQQ>(0, 0, 1, -1, 0, 0, 0, -1, 0);
IsometryQQPtr AutomorphismZZ::A0010N0100 = std::make_shared<IsometryQQ>(0, 0, 1, 0, -1, 0, 1, 0, 0);
IsometryQQPtr AutomorphismZZ::A001010N00 = std::make_shared<IsometryQQ>(0, 0, 1, 0, 1, 0, -1, 0, 0);
IsometryQQPtr AutomorphismZZ::A010N00001 = std::make_shared<IsometryQQ>(0, 1, 0, -1, 0, 0, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::A01000NN00 = std::make_shared<IsometryQQ>(0, 1, 0, 0, 0, -1, -1, 0, 0);
IsometryQQPtr AutomorphismZZ::A10000N010 = std::make_shared<IsometryQQ>(1, 0, 0, 0, 0, -1, 0, 1, 0);
IsometryQQPtr AutomorphismZZ::A1000010N0 = std::make_shared<IsometryQQ>(1, 0, 0, 0, 0, 1, 0, -1, 0);
IsometryQQPtr AutomorphismZZ::ANNN010100 = std::make_shared<IsometryQQ>(-1, -1, -1, 0, 1, 0, 1, 0, 0);
IsometryQQPtr AutomorphismZZ::ANNN100001 = std::make_shared<IsometryQQ>(-1, -1, -1, 1, 0, 0, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::AN000N0111 = std::make_shared<IsometryQQ>(-1, 0, 0, 0, -1, 0, 1, 1, 1);
IsometryQQPtr AutomorphismZZ::AN0011100N = std::make_shared<IsometryQQ>(-1, 0, 0, 1, 1, 1, 0, 0, -1);
IsometryQQPtr AutomorphismZZ::A0N000N111 = std::make_shared<IsometryQQ>(0, -1, 0, 0, 0, -1, 1, 1, 1);
IsometryQQPtr AutomorphismZZ::A0N0111N00 = std::make_shared<IsometryQQ>(0, -1, 0, 1, 1, 1, -1, 0, 0);
IsometryQQPtr AutomorphismZZ::A00NN00111 = std::make_shared<IsometryQQ>(0, 0, -1, -1, 0, 0, 1, 1, 1);
IsometryQQPtr AutomorphismZZ::A00N1110N0 = std::make_shared<IsometryQQ>(0, 0, -1, 1, 1, 1, 0, -1, 0);
IsometryQQPtr AutomorphismZZ::A001NNN100 = std::make_shared<IsometryQQ>(0, 0, 1, -1, -1, -1, 1, 0, 0);
IsometryQQPtr AutomorphismZZ::A001010NNN = std::make_shared<IsometryQQ>(0, 0, 1, 0, 1, 0, -1, -1, -1);
IsometryQQPtr AutomorphismZZ::A010NNN001 = std::make_shared<IsometryQQ>(0, 1, 0, -1, -1, -1, 0, 0, 1);
IsometryQQPtr AutomorphismZZ::A010100NNN = std::make_shared<IsometryQQ>(0, 1, 0, 1, 0, 0, -1, -1, -1);
IsometryQQPtr AutomorphismZZ::A100NNN010 = std::make_shared<IsometryQQ>(1, 0, 0, -1, -1, -1, 0, 1, 0);
IsometryQQPtr AutomorphismZZ::A100001NNN = std::make_shared<IsometryQQ>(1, 0, 0, 0, 0, 1, -1, -1, -1);
IsometryQQPtr AutomorphismZZ::A111N000N0 = std::make_shared<IsometryQQ>(1, 1, 1, -1, 0, 0, 0, -1, 0);
IsometryQQPtr AutomorphismZZ::A11100NN00 = std::make_shared<IsometryQQ>(1, 1, 1, 0, 0, -1, -1, 0, 0);
