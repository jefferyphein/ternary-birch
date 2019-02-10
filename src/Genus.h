#ifndef __GENUS_H_
#define __GENUS_H_

#include "birch.h"
#include "Math.h"
#include "HashMap.h"
#include "Spinor.h"
#include "Isometry.h"
#include "NeighborManager.h"
#include "Eigenvector.h"

template<typename R>
class GenusRep
{
public:
    GenusRep() = default;
    GenusRep(const GenusRep<R>& genus) = default;
    GenusRep(GenusRep<R>&& genus) = default;

    QuadForm<R> q;
    Isometry<R> s;
    Isometry<R> sinv;
    Z64 parent;
    R p;
    std::map<R,int> es;
};

template<typename R>
class Genus
{
    template<typename T>
    friend class Genus;

public:
    Genus() = default;

    Genus(const QuadForm<R>& q, const std::vector<PrimeSymbol<R>>& symbols, W64 seed=0)
    {
        if (seed == 0)
        {
            std::random_device rd;
            seed = rd();
        }

        this->disc = q.discriminant();
        this->seed_ = seed;

        this->prime_divisors.reserve(symbols.size());
        for (const PrimeSymbol<R>& symb : symbols)
        {
            this->prime_divisors.push_back(symb.p);
        }

        Spinor<R> *spin = new Spinor<R>(this->prime_divisors);
        this->spinor = std::unique_ptr<Spinor<R>>(spin);

        if (symbols.size() > 63)
        {
            throw std::domain_error("Must have 63 or fewer prime divisors.");
        }

        size_t num_conductors = 1LL << symbols.size();

        this->conductors.reserve(num_conductors);
        this->conductors.push_back(1);

        size_t bits = 0;
        size_t mask = 1;
        for (size_t n=1; n<num_conductors; n++)
        {
            if (n == 2*mask)
            {
                ++bits;
                mask = 1LL << bits;
            }
            R value = this->prime_divisors[bits] * this->conductors[n ^ mask];
            this->conductors.push_back(value);
        }

        GenusRep<R> rep;
        rep.q = q;
        rep.p = 1;
        rep.parent = -1;

        // Set the mass as a multiple of 24, as this is the largest integer
        // that can appear in its denominator. This value is used to determine
        // when the genus has been fully populated.
        this->mass_x24 = this->get_mass(q, symbols);

        // The mass provides a reasonable estimate for the size of the genus
        // since most isometry classes typically have trivial automorphism
        // group.
        Z64 estimated_size = ceil(mpz_get_d(this->mass_x24.get_mpz_t()) / 24.0);
        auto *ptr = new HashMap<GenusRep<R>>(estimated_size);
        this->hash = std::unique_ptr<HashMap<GenusRep<R>>>(ptr);
        this->hash->add(rep);

        // The spinor primes hash table, used to identify the primes used in
        // constructing the genus representatives.
        auto *ptr2 = new HashMap<W16>();
        this->spinor_primes = std::unique_ptr<HashMap<W16>>(ptr2);

        Z sum_mass_x24 = (48 / QuadForm<R>::num_automorphisms(q));

        Z p = 1;
        W16 prime = 1;

        // A temporary placeholder for the genus representatives before they
        // are fully built.
        GenusRep<R> foo;

        bool done = (sum_mass_x24 == this->mass_x24);
        while (!done)
        {
            // Get the next good prime and build the appropriate finite field.
            do
            {
                mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
                prime = mpz_get_ui(p.get_mpz_t());
            }
            while (this->disc % prime == 0);
            std::shared_ptr<W16_Fp> GF;
            if (prime == 2)
                GF = std::make_shared<W16_F2>(prime, this->seed_);
            else
                GF = std::make_shared<W16_Fp>(prime, this->seed_, true);

            size_t current = 0;
            while (!done && current < this->hash->size())
            {
                // Get the current quadratic form and build the neighbor manager.
                const QuadForm<R>& mother = this->hash->get(current).q;
                NeighborManager<W16,W32,R> manager(mother, GF);

                #ifdef DEBUG
                // Build the affine quadratic form for debugging purposes.
                W16_QuadForm qp = mother.mod(GF);
                #endif

                for (W16 t=0; !done && t<=prime; t++)
                {
                    #ifdef DEBUG
                    // Verify that the appropriate vector is isotropic.
                    W16_Vector3 vec = manager.isotropic_vector(t);
                    assert( qp.evaluate(vec) % prime == 0 );
                    #endif

                    // Construct the neighbor, the isometry is stored in s.
                    foo.s.set_identity();
                    foo.q = manager.get_neighbor(t, foo.s);

                    #ifdef DEBUG
                    // Verify neighbor discriminant matches.
                    assert( rep.q.discriminant() == mother.discriminant() );
                    #endif

                    // Reduce the neighbor to its Eisenstein form and add it to
                    // the hash table.
                    foo.q = QuadForm<R>::reduce(foo.q, foo.s);
                    foo.p = prime;
                    foo.parent = current;

                    bool added = this->hash->add(foo);
                    if (added)
                    {
                        const GenusRep<R>& temp = this->hash->last();
                        sum_mass_x24 += 48 / QuadForm<R>::num_automorphisms(temp.q);
                        done = (sum_mass_x24 == this->mass_x24);
                        this->spinor_primes->add(prime);
                    }
                }

                ++current;
            }
        }

        // Initialize the dimensions to zero, we will compute these values below.
        this->dims.resize(num_conductors, 0);

        // Create the lookup table values for each genus rep at each conductor.
        size_t genus_size = this->hash->size();
        this->lut_positions.resize(num_conductors, std::vector<int>(genus_size, -1));
        this->num_auts.resize(num_conductors);

        // The genus rep isometries were initialized only to contain the
        // isometry between the parent and its child, we now want to update
        // these isometries so that they are rational isometries between the
        // "mother" quadratic form and the genus rep.
        for (size_t n=0; n<this->hash->size(); n++)
        {
            GenusRep<R>& rep = this->hash->at(n);

            // Only compute composite isometries if we are not considering the
            // mother form.
            if (n)
            {
                GenusRep<R>& parent = this->hash->at(rep.parent);

                // Construct the isometries to/from the mother quadratic form.
                rep.sinv = rep.s.inverse(rep.p);
                rep.sinv = rep.sinv * parent.sinv;
                rep.s = parent.s * rep.s;

                // Copy the numerators, and increment the genus rep prime.
                rep.es = parent.es;
                ++rep.es[rep.p];

                #ifdef DEBUG
                R scalar = birch_util::my_pow(rep.es);
                scalar *= scalar;

                // Verify that s is an isometry from the mother form to the rep,
                // and that sinv is an isometry from the rep to the mother form.
                assert( rep.s.is_isometry(q, rep.q, scalar) );
                assert( rep.sinv.is_isometry(rep.q, q, scalar) );
                #endif
            }

            // Determine which subspaces this representative contributes.
            const std::vector<Isometry<R>>& auts = QuadForm<R>::proper_automorphisms(rep.q);
            std::vector<bool> ignore(this->conductors.size(), false);
            for (const Isometry<R>& s : auts)
            {
                Z64 vals = this->spinor->norm(rep.q, s, 1);

                for (size_t k=0; k<num_conductors; k++)
                {
                    if (!ignore[k] && (birch_util::popcnt(vals & k) & 1))
                    {
                        ignore[k] = true;
                    }
                }
            }

            int num = QuadForm<R>::num_automorphisms(rep.q);
            for (size_t k=0; k<num_conductors; k++)
            {
                if (!ignore[k])
                {
                    this->lut_positions[k][n] = this->dims[k];
                    this->num_auts[k].push_back(num);
                }
                this->dims[k] += (ignore[k] ? 0 : 1);
            }
        }
    }

    template<typename T>
    Genus(const Genus<T>& src)
    {
        // Convert the discriminant.
        this->disc = birch_util::convert_Integer<T,R>(src.disc);

        // Convert the prime divisors.
        for (const T& p : src.prime_divisors)
        {
            this->prime_divisors.push_back(birch_util::convert_Integer<T,R>(p));
        }

        // Convert the conductors.
        for (const T& cond : src.conductors)
        {
            this->conductors.push_back(birch_util::convert_Integer<T,R>(cond));
        }

        // Copy dimensions.
        this->dims = src.dims;

        // Copy automorphisms counts.
        this->num_auts = src.num_auts;

        // Copy lookup table dimensions.
        this->lut_positions = src.lut_positions;

        // Copy mass.
        this->mass_x24 = src.mass_x24;

        // Build a copy of the spinor primes hash table.
        this->spinor_primes = std::unique_ptr<HashMap<W16>>(new HashMap<W16>(src.spinor_primes->size()));
        for (W16 x : src.spinor_primes->keys())
        {
            this->spinor_primes->add(x);
        }

        // Build a copy of the genus representatives hash table.
        this->hash = std::unique_ptr<HashMap<GenusRep<R>>>(new HashMap<GenusRep<R>>(src.hash->size()));
        for (const GenusRep<T>& rep : src.hash->keys())
        {
            this->hash->add(birch_util::convert_GenusRep<T,R>(rep));
        }

        // Create Spinor class.
        std::vector<R> primes;
        primes.reserve(src.spinor->primes().size());
        for (const T& p : src.spinor->primes())
        {
            primes.push_back(birch_util::convert_Integer<T,R>(p));
        }
        this->spinor = std::unique_ptr<Spinor<R>>(new Spinor<R>(primes));

        // Copy seed.
        this->seed_ = src.seed_;
    }

    template<typename T>
    static Genus<T> convert(const Genus<R>& src)
    {
        return Genus<T>(src);
    }

    size_t size(void) const
    {
        return this->hash->keys().size();
    }

    W64 seed(void) const
    {
        return this->seed_;
    }

    std::map<R,size_t> dimension_map(void) const
    {
        std::map<R,size_t> temp;
        size_t num_conductors = this->conductors.size();
        for (size_t k=0; k<num_conductors; k++)
        {
            temp[this->conductors[k]] = this->dims[k];
        }
        return temp;
    }

    std::map<R,std::vector<int>> hecke_matrix_dense(const R& p) const
    {
        if (this->disc % p == 0)
        {
            throw std::invalid_argument("Prime must not divide the discriminant.");
        }
        return this->hecke_matrix_dense_internal(p);
    }

    std::map<R,std::vector<std::vector<int>>> hecke_matrix_sparse(const R& p) const
    {
        if (this->disc % p == 0)
        {
            throw std::invalid_argument("Prime must not divide the discriminant.");
        }
        return this->hecke_matrix_sparse_internal(p);
    }

    Eigenvector<R> eigenvector(const std::vector<Z32>& vec, const R& conductor) const
    {
        size_t num_conductors = this->conductors.size();
        bool found = false;

        size_t k;
        for (k=0; k<num_conductors; k++)
        {
            if (this->conductors[k] == conductor)
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            throw std::invalid_argument("Invalid conductor.");
        }

        size_t dim = this->dims[k];
        if (dim != vec.size())
        {
            throw std::invalid_argument("Eigenvector has incorrect dimension.");
        }

        size_t fulldim = this->size();

        std::vector<Z32> temp(this->size());
        const std::vector<int>& lut = this->lut_positions[k];

        for (size_t n=0; n<fulldim; n++)
        {
            if (lut[n] != -1)
            {
                temp[n] = vec[lut[n]];
            }
        }

        return Eigenvector<R>(std::move(temp), k);
    }

    std::vector<Z32> eigenvalues(EigenvectorManager<R>& vector_manager, const R& p) const
    {
        R bits16 = birch_util::convert_Integer<Z64,R>(1LL << 16);
        R bits32 = birch_util::convert_Integer<Z64,R>(1LL << 32);

        if (p == 2)
        {
            W16 prime = 2;
            std::shared_ptr<W16_F2> GF = std::make_shared<W16_F2>(prime, this->seed());
            return this->_eigenvectors<W16,W32>(vector_manager, GF, p);
        }
        else if (p < bits16)
        {
            W16 prime = birch_util::convert_Integer<R,W16>(p);
            std::shared_ptr<W16_Fp> GF = std::make_shared<W16_Fp>(prime, this->seed(), true);
            return this->_eigenvectors<W16,W32>(vector_manager, GF, p);
        }
        else if (p < bits32)
        {
            W32 prime = birch_util::convert_Integer<R,W32>(p);
            std::shared_ptr<W32_Fp> GF = std::make_shared<W32_Fp>(prime, this->seed(), false);
            return this->_eigenvectors<W32,W64>(vector_manager, GF, p);
        }
        else
        {
            W64 prime = birch_util::convert_Integer<R,W64>(p);
            std::shared_ptr<W64_Fp> GF = std::make_shared<W64_Fp>(prime, this->seed(), false);
            return this->_eigenvectors<W64,W128>(vector_manager, GF, p);
        }
    }

    const GenusRep<R>& representative(size_t n) const
    {
        return this->hash->get(n);
    }

    size_t indexof(const GenusRep<R>& rep) const
    {
        return this->hash->indexof(rep);
    }

private:
    R disc;
    std::vector<R> prime_divisors;
    std::vector<R> conductors;
    std::vector<size_t> dims;
    std::vector<std::vector<size_t>> num_auts;
    std::vector<std::vector<int>> lut_positions;
    Z mass_x24;
    std::unique_ptr<HashMap<W16>> spinor_primes;
    std::unique_ptr<HashMap<GenusRep<R>>> hash;
    std::unique_ptr<Spinor<R>> spinor;
    W64 seed_;

    template<typename S, typename T>
    std::vector<Z32> _eigenvectors(EigenvectorManager<R>& vector_manager, std::shared_ptr<Fp<S,T>> GF, const R& p) const
    {
        std::vector<Z32> eigenvalues(vector_manager.size());

        S prime = GF->prime();

        const GenusRep<R>& mother = this->hash->get(0);

        size_t num_indices = vector_manager.indices.size();
        for (size_t index=0; index<num_indices; index++)
        {
            Z64 npos = vector_manager.indices[index];
            const GenusRep<R>& cur = this->hash->get(npos);
            NeighborManager<S,T,R> neighbor_manager(cur.q, GF);
            for (W64 t=0; t<=prime; t++)
            {
                GenusRep<R> foo = neighbor_manager.get_reduced_neighbor_rep((S)t);

                size_t rpos = this->hash->indexof(foo);
                size_t offset = vector_manager.stride * rpos;
                __builtin_prefetch(&vector_manager.strided_eigenvectors[offset], 0, 0);

                W64 spin_vals;
                if (rpos == npos)
                {
                    spin_vals = this->spinor->norm(foo.q, foo.s, p);
                }
                else
                {
                    const GenusRep<R>& rep = this->hash->get(rpos);
                    foo.s = cur.s * foo.s;
                    R scalar = p;

                    foo.s = foo.s * rep.sinv;

                    scalar *= birch_util::my_pow(cur.es);
                    scalar *= birch_util::my_pow(rep.es);

                    spin_vals = this->spinor->norm(mother.q, foo.s, scalar);
                }

                for (Z64 vpos : vector_manager.position_lut[index])
                {
                    W64 cond = vector_manager.conductors[vpos];
                    Z32 value = birch_util::char_val(spin_vals & cond);
                    Z32 coord = vector_manager.strided_eigenvectors[offset + vpos];
                    if (coord)
                    {
                        eigenvalues[vpos] += (value * coord);
                    }
                }
            }

            // Divide out the coordinate associated to the eigenvector to
            // recover the actual eigenvalue.
            for (Z64 vpos : vector_manager.position_lut[index])
            {
                size_t offset = vector_manager.stride * npos;
                Z32 coord = vector_manager.strided_eigenvectors[offset + vpos];
                assert( eigenvalues[vpos] % coord == 0 );
                eigenvalues[vpos] /= coord;
            }
        }

        return eigenvalues;
    }

    // TODO: Add the actual mass formula here for reference.
    Z get_mass(const QuadForm<R>& q, const std::vector<PrimeSymbol<R>>& symbols)
    {
        Z mass = 2 * this->disc;
        Z a = q.h() * q.h() - 4 * q.a() * q.b();
        Z b = -q.a() * this->disc;

        for (const PrimeSymbol<R>& symb : symbols)
        {
            mass *= (symb.p + Math<Z>::hilbert_symbol(a, b, symb.p));
            mass /= 2;
            mass /= symb.p;
        }

        return mass;
    }

    std::map<R,std::vector<std::vector<int>>> hecke_matrix_sparse_internal(const R& p) const
    {
        size_t num_conductors = this->conductors.size();
        size_t num_primes = this->prime_divisors.size();

        std::vector<std::vector<int>> data(num_conductors);
        std::vector<std::vector<int>> indptr;
        std::vector<std::vector<int>> indices(num_conductors);

        W16 prime = birch_util::convert_Integer<R,W16>(p);

        std::shared_ptr<W16_Fp> GF;
        if (prime == 2)
            GF = std::make_shared<W16_F2>(2, this->seed());
        else
            GF = std::make_shared<W16_Fp>((W16)prime, this->seed(), true);

        std::vector<W64> all_spin_vals;
        all_spin_vals.reserve(prime+1);

        std::vector<std::vector<int>> rowdata;
        for (int dim : this->dims)
        {
            rowdata.push_back(std::vector<int>(dim));
            indptr.push_back(std::vector<int>(dim+1, 0));
        }

        const GenusRep<R>& mother = this->hash->keys()[0];
        size_t num_reps = this->size();
        for (size_t n=0; n<num_reps; n++)
        {
            const GenusRep<R>& cur = this->hash->get(n);
            NeighborManager<W16,W32,R> manager(cur.q, GF);

            for (W16 t=0; t<=prime; t++)
            {
                GenusRep<R> foo = manager.get_reduced_neighbor_rep(t);

                #ifdef DEBUG
                assert( foo.s.is_isometry(cur.q, foo.q, p*p) );
                #endif

                size_t r = this->hash->indexof(foo);

                #ifdef DEBUG
                assert( r < this->size() );
                #endif

                W64 spin_vals;
                if (r == n)
                {
                    spin_vals = this->spinor->norm(foo.q, foo.s, p);
                }
                else
                {
                    const GenusRep<R>& rep = this->hash->get(r);
                    foo.s = cur.s * foo.s;
                    R scalar = p;

                    #ifdef DEBUG
                    R temp_scalar = p*p;
                    R temp = birch_util::my_pow(cur.es);
                    temp_scalar *= temp * temp;
                    assert( foo.s.is_isometry(mother.q, foo.q, temp_scalar) );
                    #endif

                    foo.s = foo.s * rep.sinv;

                    #ifdef DEBUG
                    temp = birch_util::my_pow(rep.es);
                    temp_scalar *= temp * temp;
                    assert( foo.s.is_isometry(mother.q, mother.q, temp_scalar) );
                    #endif

                    scalar *= birch_util::my_pow(cur.es);
                    scalar *= birch_util::my_pow(rep.es);

                    #ifdef DEBUG
                    assert( scalar*scalar == temp_scalar );
                    #endif

                    spin_vals = this->spinor->norm(mother.q, foo.s, scalar);
                }

                all_spin_vals.push_back((r << num_primes) | spin_vals);
            }

            for (size_t k=0; k<num_conductors; k++)
            {
                const std::vector<int>& lut = this->lut_positions[k];
                int npos = lut[n];
                if (npos == -1) continue;

                // Populate the row data.
                std::vector<int>& row = rowdata[k];
                for (W64 x : all_spin_vals)
                {
                    int r = x >> num_primes;
                    int rpos = lut[r];
                    if (rpos == -1) continue;

                    int value = birch_util::char_val(x & k);
                    row[rpos] += value;
                }

                // Update data and indices with the nonzero values.
                size_t nnz = 0;
                size_t pos = 0;
                std::vector<int>& data_k = data[k];
                std::vector<int>& indices_k = indices[k];
                for (int x : row)
                {
                    if (x)
                    {
                        data_k.push_back(x);
                        indices_k.push_back(pos);
                        row[pos] = 0; // Clear the nonzero entry.
                        ++nnz;
                    }
                    ++pos;
                }

                // Update indptr
                indptr[k][npos+1] = indptr[k][npos] + nnz;
            }

            all_spin_vals.clear();
        }

        std::map<R,std::vector<std::vector<int>>> csr_matrices;
        for (size_t k=0; k<num_conductors; k++)
        {
            const R& cond = this->conductors[k];
            csr_matrices[cond] = std::vector<std::vector<int>>();
            csr_matrices[cond].push_back(data[k]);
            csr_matrices[cond].push_back(indices[k]);
            csr_matrices[cond].push_back(indptr[k]);
        }
        return csr_matrices;
    }

    std::map<R,std::vector<int>> hecke_matrix_dense_internal(const R& p) const
    {
        size_t num_conductors = this->conductors.size();
        size_t num_primes = this->prime_divisors.size();

        // Allocate memory for the Hecke matrices and create a vector to store
        // pointers to the raw matrix data.
        std::vector<int*> hecke_ptr;
        hecke_ptr.reserve(num_conductors);
        std::vector<std::vector<int>> hecke_matrices;
        for (size_t k=0; k<num_conductors; k++)
        {
            size_t dim = this->dims[k];
            hecke_matrices.push_back(std::vector<int>(dim * dim));
            hecke_ptr.push_back(hecke_matrices.back().data());
        }

        W16 prime = birch_util::convert_Integer<R,W16>(p);
        std::vector<W64> all_spin_vals;
        all_spin_vals.reserve(prime+1);

        std::shared_ptr<W16_Fp> GF;
        if (prime == 2)
            GF = std::make_shared<W16_F2>(2, this->seed());
        else
            GF = std::make_shared<W16_Fp>((W16)prime, this->seed(), true);

        const GenusRep<R>& mother = this->hash->keys()[0];
        size_t num_reps = this->size();

        // Create hash tables for storing isotropic vectors to be skipped
        // at later iterations.
        std::vector<HashMap<W16_Vector3>> vector_hash(num_reps);

        for (size_t n=0; n<num_reps; n++)
        {
            const GenusRep<R>& cur = this->hash->get(n);
            NeighborManager<W16,W32,R> manager(cur.q, GF);

            for (W16 t=0; t<=prime; t++)
            {
                GenusRep<R> foo;
                W16_Vector3 vec = manager.isotropic_vector(t);
                vec.x = GF->mod(vec.x);
                vec.y = GF->mod(vec.y);
                vec.z = GF->mod(vec.z);

                // If this vector has already been identified, skip it. This
                // avoids unnecessarily computing the neighbor, reducing it,
                // and testing for isometry. The Hermitian symmetry property
                // of the Hecke matrix will account for this once we finish
                // processing neighbors.
                if (vector_hash[n].exists(vec)) continue;

                // Build the neighbor and reduce it.
                foo.q = manager.build_neighbor(vec, foo.s);
                foo.q = QuadForm<R>::reduce(foo.q, foo.s);

                #ifdef DEBUG
                assert( foo.s.is_isometry(cur.q, foo.q, p*p) );
                #endif

                size_t r = this->hash->indexof(foo);

                #ifdef DEBUG
                assert( r < this->size() );
                #endif

                W64 spin_vals;
                if (unlikely(r == n))
                {
                    spin_vals = this->spinor->norm(foo.q, foo.s, p);
                }
                else if (r > n)
                {
                    W16_Vector3 result = manager.transform_vector(cur, foo, vec);
                    vector_hash[r].add(result);

                    const GenusRep<R>& rep = this->hash->get(r);
                    foo.s = cur.s * foo.s;
                    R scalar = p;

                    #ifdef DEBUG
                    R temp_scalar = p*p;
                    R temp = birch_util::my_pow(cur.es);
                    temp_scalar *= temp * temp;
                    assert( foo.s.is_isometry(mother.q, foo.q, temp_scalar) );
                    #endif

                    foo.s = foo.s * rep.sinv;

                    #ifdef DEBUG
                    temp = birch_util::my_pow(rep.es);
                    temp_scalar *= temp * temp;
                    assert( foo.s.is_isometry(mother.q, mother.q, temp_scalar) );
                    #endif

                    scalar *= birch_util::my_pow(cur.es);
                    scalar *= birch_util::my_pow(rep.es);

                    #ifdef DEBUG
                    assert( scalar*scalar == temp_scalar );
                    #endif

                    spin_vals = this->spinor->norm(mother.q, foo.s, scalar);
                }
                else continue;

                all_spin_vals.push_back((r << num_primes) | spin_vals);
            }

            for (size_t k=0; k<num_conductors; k++)
            {
                const std::vector<int>& lut = this->lut_positions[k];
                int npos = lut[n];
                if (unlikely(npos == -1)) continue;

                int *row = hecke_ptr[k];

                for (W64 x : all_spin_vals)
                {
                    int r = x >> num_primes;
                    int rpos = lut[r];
                    if (unlikely(rpos == -1)) continue;

                    row[rpos] += birch_util::char_val(x & k);
                }

                hecke_ptr[k] += this->dims[k];
            }

            all_spin_vals.clear();
        }

        // Copy the upper diagonal entries to the lower diagonal using the
        // Hermitian symmetry property and then move the matrix into an
        // associatively map before returning.
        std::map<R,std::vector<int>> matrices;
        for (size_t k=0; k<num_conductors; k++)
        {
            std::vector<int>& matrix = hecke_matrices[k];
            size_t dim = this->dims[k];
            size_t dim2 = dim * dim;
            const std::vector<size_t>& auts = this->num_auts[k];

            // Copy upper diagonal matrix to the lower diagonal.
            for (size_t start=0, row=0; start<dim2; start+=dim+1, row++)
            {
                int row_auts = auts[row];
                for (size_t dst=start+dim, src=start+1, col=row+1; col<dim; src++, col++, dst+=dim)
                {
                    if (matrix[src])
                    {
                        int col_auts = auts[col];
                        if (col_auts == row_auts)
                        {
                            matrix[dst] = matrix[src];
                        }
                        else
                        {
                            matrix[dst] = matrix[src] * col_auts / row_auts;
                        }
                    }
                }
            }

            // Move the matrix in the corresponding entry in the map.
            matrices[this->conductors[k]] = std::move(hecke_matrices[k]);
        }
        return matrices;
    }
};

template<typename R>
bool operator==(const GenusRep<R>& a, const GenusRep<R>& b)
{
    return a.q == b.q;
}

#endif // __GENUS_H_
