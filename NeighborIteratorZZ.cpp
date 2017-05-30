#include <memory>
#include <iostream>
#include <iomanip>
#include "NeighborIterator.h"
#include "QuadForm.h"
#include "Math.h"

typedef NeighborIterator<mpz_class, mpq_class> NeighborIteratorZZ;
typedef QuadForm<mpz_class, mpq_class> QuadFormZZ;
typedef Prime<mpz_class, mpq_class> PrimeZZ;
typedef Isometry<mpz_class, mpq_class> IsometryZZ;

template<>
NeighborIteratorZZ::NeighborIterator(std::shared_ptr<QuadFormZZ> q, std::shared_ptr<PrimeZZ> pR) :
    q_(q),
    p_(pR),
    numNeighbors_(mpz_class(pR->norm() + 1)),
    pos_(mpz_class(0)),
    isotropicVector_(q->isotropic_vector(pR)),
    s_(std::make_shared<IsometryZZ>(false))
{
//    mpz_class p = this->p_->principal_generator();
//
//    if (this->isotropicVector_[0] != 0)
//    {
//        mpz_class inv = Math::modinv(this->isotropicVector_[0], p);
//        this->isotropicVector_[0] = 1;
//        this->isotropicVector_[1] = (this->isotropicVector_[1] * inv) % p;
//        this->isotropicVector_[2] = (this->isotropicVector_[2] * inv) % p;
//
//        mpz_class u(this->isotropicVector_[1]);
//        mpz_class v(this->isotropicVector_[2]);
//
//        this->s_->a11 = 1;
//        this->s_->a21 = u;
//        this->s_->a22 = 1;
//        this->s_->a31 = v;
//        this->s_->a33 = 1;
//    }
//    else if (this->isotropicVector_[1] != 0)
//    {
//        mpz_class inv = Math::modinv(this->isotropicVector_[1], p);
//        this->isotropicVector_[1] = 1;
//        this->isotropicVector_[2] = (this->isotropicVector_[1] * inv) % p;
//
//        mpz_class u = this->isotropicVector_[2];
//
//        this->s_->a12 = 1;
//        this->s_->a21 = 1;
//        this->s_->a31 = u;
//        this->s_->a33 = 1;
//    }
//    else
//    {
//        this->isotropicVector_[2] = 1;
//
//        this->s_->a12 = 1;
//        this->s_->a23 = 1;
//        this->s_->a31 = 1;
//    }
}

template<>
std::shared_ptr<QuadFormZZ> NeighborIteratorZZ::next_neighbor(void)
{
    // Return nullptr if no more neighbors to compute.
    if (this->pos_ >= this->numNeighbors_) { return nullptr; }

    mpz_class p = this->p_->principal_generator();

    std::vector<mpz_class> vec(this->isotropicVector_);

    std::cout << "Q_{" << p << "}(" << vec[0] << ","
              << vec[1] << "," 
              << vec[2] << ") = " << (this->q_->evaluate(vec) % p) << std::endl;

    ++this->pos_;

    return this->q_;
}
