#ifndef __BIRCH_H_
#define __BIRCH_H_

#include <stdint.h>
#include <iostream>
#include <chrono>
#include <cassert>
#include <vector>
#include <map>
#include <array>
#include <random>
#include <memory>
#include <gmpxx.h>

/* Type definitions */

typedef mpz_class W;
typedef uint16_t W16;
typedef uint32_t W32;
typedef uint64_t W64;
typedef __uint128_t W128;

typedef mpz_class Z;
typedef int16_t Z16;
typedef int32_t Z32;
typedef int64_t Z64;
typedef __int128_t Z128;

/* Builtins */

#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

// Constants
constexpr W64 FNV_OFFSET = 0x811c9dc5;
constexpr W64 FNV_PRIME  = 0x01000193;

/* Class declarations */

template<typename R>
class Math;

template<typename R>
class Isometry;

template<typename R>
class QuadForm;

template<typename R, typename S>
class QuadFormFp;

template<typename R, typename S>
class Fp;

template<typename R, typename S>
class F2;

template<typename R>
class Eigenvector;

template<typename R>
class EigenvectorContainer;

template<typename R, typename S, typename T>
class NeighborManager;

template<typename Key>
class HashMap;

template<typename R>
class Genus;

template<typename R>
class Spinor;

template<typename R>
class GenusRep;

template<typename R, typename S, typename T>
class IsometrySequence;

class SetCover;

/* Struct definitions */

template<typename R>
struct Vector3 {
    R x; R y; R z;
};

namespace std
{
    template<typename R>
    struct hash<Vector3<R>>
    {
        Z64 operator()(const Vector3<R>& vec) const
        {
            Z64 fnv = FNV_OFFSET;
            fnv = (fnv ^ vec.x) * FNV_PRIME;
            fnv = (fnv ^ vec.y) * FNV_PRIME;
            fnv = (fnv ^ vec.z) * FNV_PRIME;
            return fnv;
        }
    };
}

template<typename R>
struct PrimeSymbol {
    R p;
    int power;
    bool ramified;
};

/* Templated class types */

// Isometries.
typedef Isometry<Z> Z_Isometry;
typedef Isometry<Z64> Z64_Isometry;

// Quadratic forms over the integers.
typedef QuadForm<Z64> Z64_QuadForm;
typedef QuadForm<Z>   Z_QuadForm;

// Quadratic forms over a finite field.
typedef QuadFormFp<W16,W32>  W16_QuadForm;
typedef QuadFormFp<W32,W64>  W32_QuadForm;
typedef QuadFormFp<W64,W128> W64_QuadForm;

// Vectors.
typedef Vector3<W16> W16_Vector3;
typedef Vector3<W32> W32_Vector3;
typedef Vector3<W64> W64_Vector3;
typedef Vector3<Z64> Z64_Vector3;
typedef Vector3<Z>   Z_Vector3;

// Finite fields.
typedef Fp<W16,W32>  W16_Fp;
typedef Fp<W32,W64>  W32_Fp;
typedef Fp<W64,W128> W64_Fp;
typedef F2<W16,W32>  W16_F2;

// Prime symbols
typedef PrimeSymbol<Z>   Z_PrimeSymbol;
typedef PrimeSymbol<Z64> Z64_PrimeSymbol;

// Math.
typedef Math<Z> Z_Math;
typedef Math<Z64> Z64_Math;

// Neighbor managers.
typedef NeighborManager<W16,W32,Z>  Z_W16_NeighborManager;
typedef NeighborManager<W32,W64,Z>  Z_W32_NeighborManager;
typedef NeighborManager<W64,W128,Z> Z_W64_NeighborManager;
typedef NeighborManager<W16,W32,Z64>  Z64_W16_NeighborManager;
typedef NeighborManager<W32,W64,Z64>  Z64_W32_NeighborManager;
typedef NeighborManager<W64,W128,Z64> Z64_W64_NeighborManager;

// Genus
typedef Genus<Z64> Z64_Genus;
typedef Genus<Z>   Z_Genus;

// Genus representatives
typedef GenusRep<Z64> Z64_GenusRep;
typedef GenusRep<Z> Z_GenusRep;

#endif // __BIRCH_H_
