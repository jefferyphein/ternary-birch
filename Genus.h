#ifndef __GENUS_H_
#define __GENUS_H_

#include <iostream>
#include <memory>
#include <cassert>
#include <unordered_set>
#include <set>
#include <map>
#include <queue>
#include <thread>
#include <mutex>
#include <chrono>
#include <fstream>
#include <future>
#include <atomic>
#include <random>
#include "QuadForm.h"
#include "NeighborIterator.h"
#include "Character.h"
#include "HeckeOperator.h"
#include "Eigenvector.h"
#include "Math.h"

template<typename R, typename F>
class GenusRep
{
public:
    typedef std::shared_ptr<QuadForm<R,F>> QuadFormPtr;
    typedef std::shared_ptr<Isometry<R,F>> IsometryPtr;

    GenusRep(QuadFormPtr qPtr);

    int64_t pos(const Character<R,F>& chi) const;
    void set_position(const Character<R,F>& chi, int64_t x);
    void set_absolute_position(int64_t x);
    int64_t position(const Character<R,F>& chi) const;
    int64_t position(const R& cond) const;
    int64_t absolute_position(void) const;
    void set_inverse(void);
    QuadFormPtr quad_form(void) const;
    IsometryPtr inverse(void) const;

    void set_dimension(const Character<R,F>& chi, int64_t dim);
    int64_t dimension(const Character<R,F>& chi) const;

    bool operator==(const GenusRep<R,F>& rep) const;
    bool operator<(const GenusRep<R,F>& rep) const;
private:
    QuadFormPtr q_;
    int64_t pos_;
    IsometryPtr inv_;
    std::map<R, int64_t> positionMap_;
    std::map<R, int64_t> dimensionMap_;
};

template<typename R, typename F>
class Genus
{
public:
    typedef std::shared_ptr<QuadForm<R,F>> QuadFormPtr;
    typedef std::shared_ptr<HeckeOperator<R,F>> HeckePtr;

    Genus(const QuadForm<R,F>& q, int64_t numThreads=0);
    Genus(const std::string& filename);

    bool computed(void) const;

    QuadFormPtr quad_form(void) const;

    const R& discriminant(void) const;

    R smallest_good_prime(void) const;
    std::vector<R> bad_primes(void) const;

    void add_character(const Character<R,F>& chi);

    void set_dimensions(GenusRep<R,F>& rep);

    void print(std::ostream& os) const;

    size_t size(void) const;

    int64_t dimension(const Character<R,F>& chi) const;
    int64_t dimension(const R& cond) const;

    void compute_hecke_operators(const R& p, int64_t numThreads=0);

    HeckePtr hecke_operator(const R& p, const Character<R,F>& chi);

    void import_eigenvectors(const std::string& filename);
    void export_eigenvectors(const std::string& filename);

    void compute_eigenvalues(const R& p, int64_t numThreads=0);
    void compute_eigenvalues(const std::vector<R>& ps, int64_t numThreads=0);

private:
    ////////// GENERAL MEMBER VARIABLES

    /* A pointer to the originating form. */
    QuadFormPtr q_;

    /* The discriminant. */
    R disc_;

    /* Flag indicating whether the genus has been computed in full. */
    bool computed_ = false;

    /* Character dimensions. */
    std::map<R, int64_t> dimensionMap_;

    /* Attempts to add a new genus representative. */
    void add_genus_rep(QuadFormPtr q, QuadFormPtr neighbor=nullptr, bool force=false);

    /* Finds a genus representative. Will throw an exception if not found. */
    const GenusRep<R,F>& find_genus_rep(QuadFormPtr q) const;

    /* A filename used for importing the genus representatives. */
    std::string genusFilename;

    /* The number of threads to use when computing the genus. */
    int64_t genusNumThreads = 0;

    /* Creates the genus from a file. */
    void import_genus(const std::string& filename);

    /* Computes the full genus. */
    void compute_genus(int64_t numThreads=0);

    ////////// THREADING RESOURCES

    /* This function is passed to threads, pulls a genus representative off of
     * a queue and computes its p-neighbors. */
    void threaded_compute_genus(const R& p, int64_t threadid);

    /* This function is passed to threads, pulls a genus representative out of
     * a vector based on the value of genusVecIndex_. */
    void threaded_compute_hecke_operators(const R& p, int64_t threadid);

    /* This function is passed to threads, accesses p-neighbors, and computes eigenvalues. */
    void threaded_compute_eigenvalues(std::promise<std::map<int64_t, std::map<R, int64_t>>>& aps);

    /* A vector of booleans used to determine whether a thread is blocked due
     * to there being no objects in the genusQueue_. */
    std::vector<bool> threadBlocked_;

    /* A queue of shared QuadForm pointers used for multithreading. */
    std::queue<QuadFormPtr> genusQueue_;

    /* A thread-safe value used to access GenusReps when computing Hecke
     * operators. */
    std::atomic<int64_t> genusVecIndex_;

    /* A mutex used for accessing genus reps when multithreading. */
    std::mutex genusMutex_;

    /* A mutex used for updating Hecke operator values when multithreading. */
    std::mutex heckeMutex_;

    /* A neighbor iterator object to be shared amongst threads when computing
     * eigenvalues. */
    std::vector<NeighborIterator<R,F>> eigenvalueNeighborVec_;

    /* A thread-safe value used to iterate over neighbors when computing
     * Hecke eigenvalues. */
    std::atomic<int64_t> neighborIndex_;

    ////////// GENUS DATA STRUCTURES

    /* An unordered set which stores distinct genus representatives. */
    std::unordered_set<GenusRep<R,F>> genusReps_;

    /* A vector of shared QuadForm pointers. An unordered list of the genus. */
    std::vector<QuadFormPtr> genusVec_;

    /* A set of characters to compute. */
    std::set<Character<R,F>> charSet_;

    /* A map from conductor to character. */
    std::map<R, Character<R,F>> condToChar_;

    /* A set of prime characters to compute. When representation
     * values are computed for characters in charSet_, these will first be
     * computed and then multiplied together to obtain the desired
     * representation value. */
    std::set<Character<R,F>> primeCharSet_;

    ////////// HECKE OPERATOR RESOURCES

    /* A map structure used to organize Hecke operators. */
    std::map<R, std::map<R, HeckePtr>> heckeMap_;

    /* A helper function that updates Hecke operators. */
    void update_hecke_operators(const GenusRep<R,F>& rep,
                                QuadFormPtr neighbor,
                                std::map<R, std::map<int64_t, int64_t>>& rowMap);

    ////////// EIGENVECTOR RESOURCES

    /* A function used to compute the pivots used for computing eigenvalues.
     * This is an implementation of a greedy algorithm. */
    void assign_eigenvector_pivots_greedy(void);

    /* A map of Eigenvectors, associated by conductors. */
    std::map<Character<R,F>, std::vector<Eigenvector>> eigenvectorMap_;

    /* A map of lookup positions used to associate eigenvector positions to
     * their associated genus representative. */
    std::map<R, std::vector<int64_t>> absolutePosition_;

    /* A vector of pivot positions used to compute Hecke eigenvalues. */
    std::vector<int64_t> eigenvectorPivots_;

    /* A two-dimensional map that stores eigenvalues associated to Eigenvector
     * pointers and prime values. */
    std::map<int64_t, std::map<R, int64_t>> eigenvalueMap_;

    /* A map of pivots used for each eigenvector. */
    std::map<int64_t, int64_t> pivotMap_;
};

template<typename R, typename F>
GenusRep<R,F>::GenusRep(QuadFormPtr q)
{
    this->q_ = q;
    this->inv_ = nullptr;
}

template<typename R, typename F>
void GenusRep<R,F>::set_inverse(void)
{
    this->inv_ = this->q_->isometry()->inverse();
}

template<typename R, typename F>
void GenusRep<R,F>::set_position(const Character<R,F>& chi, int64_t x)
{
    this->positionMap_[chi.conductor()] = x;
}

template<typename R, typename F>
void GenusRep<R,F>::set_absolute_position(int64_t x)
{
    this->pos_ = x;
}

template<typename R, typename F>
int64_t GenusRep<R,F>::position(const Character<R,F>& chi) const
{
    return this->positionMap_.find(chi.conductor())->second;
}

template<typename R, typename F>
int64_t GenusRep<R,F>::position(const R& cond) const
{
    return this->positionMap_.find(cond)->second;
}

template<typename R, typename F>
int64_t GenusRep<R,F>::absolute_position(void) const
{
    return this->pos_;
}

template<typename R, typename F>
std::shared_ptr<QuadForm<R,F>> GenusRep<R,F>::quad_form(void) const
{
    return this->q_;
}

template<typename R, typename F>
bool GenusRep<R,F>::operator==(const GenusRep<R,F>& rep) const
{
    return *(rep.quad_form()) == *this->q_;
}

template<typename R, typename F>
std::shared_ptr<Isometry<R,F>> GenusRep<R,F>::inverse(void) const
{
    return this->inv_;
}

template<typename R, typename F>
bool GenusRep<R,F>::operator<(const GenusRep<R,F>& rep) const
{
    return *this->q_ < *(rep.quad_form());
}

template<typename R, typename F>
int64_t GenusRep<R,F>::dimension(const Character<R,F>& chi) const
{
    const R& cond = chi.conductor();
    if (this->dimensionMap_.count(cond) > 0)
    {
        return this->dimensionMap_.find(cond)->second;
    }
    else
    {
        return 0;
    }
}

namespace std
{
    template<typename R, typename F>
    struct hash<GenusRep<R,F>>
    {
        int64_t operator()(const GenusRep<R,F>& rep) const
        {
            return rep.quad_form()->hash_value();
        }
    };
}

template<typename R, typename F>
Genus<R,F>::Genus(const QuadForm<R,F>& q, int64_t numThreads)
{
    this->q_ = QuadForm<R,F>::reduce(q, false);
    this->disc_ = this->q_->discriminant();

    if (!Math<R,F>::is_squarefree(this->disc_))
    {
        throw std::runtime_error("Discriminant must be squarefree.");
    }

    if (!Math<R,F>::is_positive(this->disc_))
    {
        throw std::runtime_error("Discriminant must be positive.");
    }

    this->genusNumThreads = numThreads;
}

template<typename R, typename F>
Genus<R,F>::Genus(const std::string& filename)
{
    // Set genus input filename.
    this->genusFilename = filename;

    // Open the genus file, so that we can read the initial quadratic form.
    std::ifstream infile(filename.c_str(), std::ios::in);
    if (!infile.is_open())
    {
        throw std::runtime_error("Unable to open input file.");
    }

    // Read the dimension. We don't need this now.
    int64_t temp;
    infile >> temp;

    // Read the coefficients of the initial form.
    R a, b, c, f, g, h;
    infile >> a >> b >> c >> f >> g >> h;
    infile.close();

    // Set the quadratic form and discriminant.
    this->q_ = std::make_shared<QuadForm<R,F>>(a, b, c, f, g, h, true);
    this->disc_ = this->q_->discriminant();
}

template<typename R, typename F>
bool Genus<R,F>::computed(void) const
{
    return this->computed_;
}

template<typename R, typename F>
std::shared_ptr<QuadForm<R,F>> Genus<R,F>::quad_form(void) const
{
    return this->q_;
}

template<typename R, typename F>
const R& Genus<R,F>::discriminant(void) const
{
    return this->disc_;
}

template<typename R, typename F>
std::vector<R> Genus<R,F>::bad_primes(void) const
{
    return Math<R,F>::prime_divisors_naive(this->disc_);
}

template<typename R, typename F>
size_t Genus<R,F>::size(void) const
{
    return this->genusReps_.size();
}

template<typename R, typename F>
void Genus<R,F>::threaded_compute_genus(const R& p, int64_t threadid)
{
    while (true)
    {
        // Obtain the mutex.
        this->genusMutex_.lock();

        // Check the size of the queue to determine whether we are blocked.
        if (this->genusQueue_.size() == 0)
        {
            // If there is nothing in the queue, we're blocked.
            this->threadBlocked_[threadid] = true;

            // Determine whether all threads are blocked.
            bool done = true;
            for (bool value : this->threadBlocked_)
            {
                done = done && value;
            }

            // If all threads are blocked, release the lock and terminate.
            // This can only happen if there are no more elements in the queue
            // to process.
            if (done) 
            {
                this->genusMutex_.unlock();
                return;
            }

            // Release the lock and sleep for a random number of milliseconds.
            this->genusMutex_.unlock();
            std::default_random_engine dre(threadid);
            std::uniform_int_distribution<int64_t> id(0, 20);
            std::this_thread::sleep_for(std::chrono::milliseconds(id(dre)));
            continue;
        }
        else
        {
            // Indicate that the thread is not blocked. This will ensure that
            // other threads will not terminate if they obtain the lock.
            this->threadBlocked_[threadid] = false;

            // Pop the next genus representative off the queue, and unlock the
            // mutex.
            QuadFormPtr cur = this->genusQueue_.front();
            this->genusQueue_.pop();
            this->genusMutex_.unlock();

            // The neighbor iterator.
            NeighborIterator<R,F> it(cur, p);

            // The first neighbor.
            QuadFormPtr pn = it.next_neighbor();
#ifdef DEBUG
            assert( pn->isometry()->is_isometry(*this->q_, *pn) );
#endif

            // Loop over all p-neighbors.
            while (pn != nullptr)
            {
#ifdef DEBUG
                assert( pn->isometry()->is_isometry(*this->q_, *pn) );
#endif

                // Reduce the p-neighbor.
                QuadFormPtr qq = QuadForm<R,F>::reduce(*pn);

                // Obtain the lock, add this genus representative to the queue
                // and release the lock.
                this->genusMutex_.lock();
                this->add_genus_rep(qq, pn);
                this->genusMutex_.unlock();

                // Get the next p-neighbor.
                pn = it.next_neighbor();
            }
        }
    }
}

template<typename R, typename F>
void Genus<R,F>::compute_genus(int64_t numThreads)
{
    // Do nothing if the genus has already been computed.
    if (this->computed_)
    {
        return;
    }

    // Initialize the dimensions for each character to zero.
    for (auto& chi : this->charSet_)
    {
        this->dimensionMap_[chi.conductor()] = 0;
    }

    // Read genus from file, if provided.
    if (!this->genusFilename.empty())
    {
        this->import_genus(this->genusFilename);
        return;
    }

    // Add the defining quadratic form as the initial genus rep.
    this->add_genus_rep(this->q_);

    // Determine the smallest good prime.
    R p = this->smallest_good_prime();

    // The number of p-neighbors this prime will produce.
    int64_t numNeighbors = NeighborIterator<R,F>::num_neighbors(p);

    // Determine the number of threads to use, up to the specified value. This
    // allows for faster genus computation in some cases where the prime only
    // provides a small number of p-neighbors.
    if (numNeighbors <= 6)
    {
        // When there are only a few p-neighbors, do not multithread.
        numThreads = 0;
    }

    // Should we use threads?
    if (numThreads > 0)
    {
        // Initialize a vector of booleans to false. Each thread will modify
        // this vector to indicate that it is blocked on the queue. Once all
        // threads have indicated that they are blocked, each thread will
        // terminate execution.
        this->threadBlocked_ = std::vector<bool>(numThreads, false);

        // Build the threads.
        std::vector<std::thread> threads;
        for (int64_t n = 0; n < numThreads; n++)
        {
            threads.push_back(std::thread(
                &Genus<R,F>::threaded_compute_genus, this, p, n
            ));
            this->threadBlocked_[n] = false;
        }

        // Wait for all threads to terminate.
        for (auto& t : threads)
        {
            t.join();
        }
    }
    else
    {
        // Loop over all genus representatives as the genus is built.
        size_t index = 0;
        while (index < this->genusVec_.size())
        {
            // The current genus representative.
            QuadFormPtr cur = this->genusVec_[index++];

            // The neighbor iterator.
            NeighborIterator<R,F> it(cur, p);

            // The first neighbor.
            QuadFormPtr pn = it.next_neighbor();

#ifdef DEBUG
            assert( pn->isometry()->is_isometry(*this->q_, *pn) );
#endif

            // Loop over all p-neighbors.
            while (pn != nullptr)
            {
#ifdef DEBUG
                assert( pn->isometry()->is_isometry(*this->q_, *pn) );
#endif

                // Reduce the p-neighbor.
                QuadFormPtr qq = QuadForm<R,F>::reduce(*pn);

                // Add this reduce form to the genus, if necessary.
                this->add_genus_rep(qq, pn);

                // Get the next p-neighbor.
                pn = it.next_neighbor();
            }
        }
    }

    // Flag the genus as computed.
    this->computed_ = true;
}

template<typename R, typename F>
const GenusRep<R,F>& Genus<R,F>::find_genus_rep(QuadFormPtr q) const
{
    GenusRep<R,F> rep(q);
    return *(this->genusReps_.find(rep));
}

template<typename R, typename F>
void Genus<R,F>::add_genus_rep(QuadFormPtr q, QuadFormPtr neighbor, bool force)
{
    // Create a GenusRep object from the shared pointer.
    GenusRep<R,F> rep(q);

    // Is this quadratic form already in the genus? If not, add it.
    if (force || this->genusReps_.count(rep) == 0)
    {
        // Compose to obtain a global isometry between this genus rep and
        // the original quadratic form.
        if (neighbor != nullptr)
        {
            q->isometry()->multiply_on_left_by(neighbor->isometry());

#ifdef DEBUG
            assert( q->isometry()->is_isometry(*this->q_, *q) );
#endif
        }

        // Set the dimensions for this genus representative.
        this->set_dimensions(rep);

        for (auto& chi : this->charSet_)
        {
            const R& cond = chi.conductor();

            // Update the dimensions for each character.
            this->dimensionMap_[cond] += rep.dimension(chi);

            // Set the relative position for this genus rep with respect to
            // each character. A position of -1 indicates that this genus rep
            // does not contribute to the dimension of the space associated to
            // this character.
            if (rep.dimension(chi) > 0)
            {
                rep.set_position(chi, this->dimensionMap_[cond]-1);
            }
            else
            {
                rep.set_position(chi, -1);
            }
        }

        // Set the absolute position of this genus rep.
        rep.set_absolute_position(this->genusVec_.size());

        // Set the inverse value of the genus rep.
        rep.set_inverse();

        // Update data structures.
        this->genusReps_.insert(rep);
        this->genusVec_.push_back(q);
        this->genusQueue_.push(q);
    }
}

template<typename R, typename F>
void GenusRep<R,F>::set_dimension(const Character<R,F>& chi, int64_t dim)
{
    this->dimensionMap_[chi.conductor()] = dim;
}

template<typename R, typename F>
void Genus<R,F>::add_character(const Character<R,F>& chi)
{
    if (this->disc_ % chi.conductor() != 0)
    {
        throw std::runtime_error("Conductor must divide discriminant.");
    }

    // If this is a new character, add it to the set.
    if (this->charSet_.count(chi) == 0)
    {
        this->charSet_.insert(chi);
        this->condToChar_[chi.conductor()] = chi;

        // The list of primes associated to this character.
        std::vector<R> ps = chi.primes();

        // Loop over all of these primes and add any new primes to the prime
        // set.
        for (auto& p : ps)
        {
            Character<R,F> temp(p);
            if (this->primeCharSet_.count(temp) == 0)
            {
                this->primeCharSet_.insert(temp);
            }
        }
    }
}

template<typename R, typename F>
std::shared_ptr<HeckeOperator<R,F>> Genus<R,F>::hecke_operator(const R& p, const Character<R,F>& chi)
{
    auto ptr = std::make_shared<HeckeOperator<R,F>>(*this, chi);
    return ptr;
}

template<typename R, typename F>
void Genus<R,F>::update_hecke_operators(const GenusRep<R,F>& rep,
                                        QuadFormPtr neighbor,
                                        std::map<R, std::map<int64_t, int64_t>>& rowMap)
{
#ifdef DEBUG
    assert( neighbor->isometry()->is_isometry(*this->q_, *neighbor) );
#endif

    // The reduced form of the p-neighbor.
    QuadFormPtr qq = QuadForm<R,F>::reduce(*neighbor);

    // The genus rep for this p-neighbor.
    const GenusRep<R,F>& rrep = this->find_genus_rep(qq);

    // Multiply the p-neighbor isometry on the right by the reduction
    // isometry. This now represents an isometry from the original
    // quadratic form to the genus representative isometric to this
    // p-neighbor.
    neighbor->isometry()->multiply_on_right_by(qq->isometry());

#ifdef DEBUG
    assert( neighbor->isometry()->is_isometry(*this->q_, *rrep.quad_form()) );
#endif

    // Obtain an automorphism of the original quadratic form.
    neighbor->isometry()->multiply_on_right_by(rrep.inverse());

#ifdef DEBUG
    assert( neighbor->isometry()->is_automorphism(*this->q_) );
#endif

    // Convenient reference for the automorphism we'll use.
    std::shared_ptr<Isometry<R,F>> aut = neighbor->isometry();

    // Compute the value of each primitive character.
    std::map<R, int64_t> primeValues;
    for (auto& chi : this->primeCharSet_)
    {
        primeValues[chi.conductor()] = chi.rho(*aut, *this->q_);
    }

    // Build the values of the composite characters from the values of the
    // primitive characters we just computed.
    for (auto& chi : this->charSet_)
    {
        int64_t row = rep.position(chi);
        if (row == -1) { continue; }

        int64_t col = rrep.position(chi);
        if (col == -1) { continue; }

        int64_t value = 1;
        auto facs = chi.primes();
        for (auto& fac : facs)
        {
            value *= primeValues[fac];
        }

        rowMap[chi.conductor()][col] += value;
    }
}

template<typename R, typename F>
void Genus<R,F>::threaded_compute_hecke_operators(const R& p, int64_t)
{
    while (true)
    {
        // Obtain the mutex and determine the index of the genus rep this
        // thread should process.
        size_t index = this->genusVecIndex_++;

        // Have we reached the end of the genus vector?
        if (index >= this->genusVec_.size())
        {
            // If so, release the lock and return.
            return;
        }

        // Obtain a genus representative form to process and release the lock.
        QuadFormPtr cur = this->genusVec_[index];

        // Create a row map for each character. This will store temporary
        // values which will then be updated all at once when we finish
        // processing each p-neighbor.
        std::map<R, std::map<int64_t, int64_t>> rowMap;
        for (auto& chi : this->charSet_)
        {
            rowMap[chi.conductor()].clear();
        }

        // The current genus representative object.
        const GenusRep<R,F>& rep = this->find_genus_rep(cur);

        // The neighbor iterator.
        NeighborIterator<R,F> it(cur, p);

        // The first p-neighbor.
        QuadFormPtr pn = it.next_neighbor();

        // Loop over all p-neighbors.
        while (pn != nullptr)
        {
            // Process the p-neighbor and update all Hecke operators.
            this->update_hecke_operators(rep, pn, rowMap);

            // The next p-neighbor.
            pn = it.next_neighbor();
        }

        // Obtain the Hecke lock, and update this row of the Hecke operator.
        for (auto& chi : this->charSet_)
        {
            const R& cond = chi.conductor();
            int64_t pos = rep.position(chi);
            if (pos != -1)
            {
                this->heckeMutex_.lock();
                this->heckeMap_[cond][p]->update_row(pos, rowMap[cond]);
                this->heckeMutex_.unlock();
            }
        }
    }
}

template<typename R, typename F>
void Genus<R,F>::compute_hecke_operators(const R& p, int64_t numThreads)
{
    if (!this->computed_)
    {
        this->compute_genus(this->genusNumThreads);
    }

    bool needToCompute = false;
    for (auto& chi : this->charSet_)
    {
        const R& cond = chi.conductor();

        // Create map for this conductor, if necessary.
        if (this->heckeMap_.count(cond) == 0)
        {
            this->heckeMap_[cond] = std::move(std::map<R, HeckePtr>());
        }

        // Create Hecke operator for this prime, if necessary.
        if (this->heckeMap_[cond].count(p) == 0)
        {
            this->heckeMap_[cond][p] = std::make_shared<HeckeOperator<R,F>>(*this, chi);
            needToCompute = true;
        }
    }

    if (!needToCompute) { return; }

    if (numThreads > 0)
    {
        // Initialize the shared pointer used by threads to compute the Hecke
        // operators.
        this->genusVecIndex_.store(0);

        // Create a vector of threads.
        std::vector<std::thread> threads;

        for (int64_t n = 0; n < numThreads; n++)
        {
            threads.push_back(std::thread(
                &Genus<R,F>::threaded_compute_hecke_operators, this, p, n
            ));
        }

        // Wait for threads to finish before moving along.
        for (auto& t : threads)
        {
            t.join();
        }
    }
    else
    {
        for (auto& cur : this->genusVec_)
        {
            // The neighbor iterator.
            NeighborIterator<R,F> it(cur, p);

            const GenusRep<R,F>& rep = this->find_genus_rep(cur);

            // Create a row map for each character.
            std::map<R, std::map<int64_t, int64_t>> rowMap;
            for (auto& chi : this->charSet_)
            {
                rowMap[chi.conductor()].clear();
            }

            // The first p-neighbor.
            QuadFormPtr pn = it.next_neighbor();

            // Loop over all p-neighbors.
            while (pn != nullptr)
            {
                // Process the p-neighbor and update all Hecke operators.
                this->update_hecke_operators(rep, pn, rowMap);

                // The next p-neighbor.
                pn = it.next_neighbor();
            }

            // Add the current row to the Hecke operator.
            for (auto& chi : this->charSet_)
            {
                const R& cond = chi.conductor();
                int64_t pos = rep.position(chi);
                if (pos != -1)
                {
                    this->heckeMap_[cond][p]->add_row(rep.position(chi), rowMap[cond]);
                }
            }
        }
    }
}

template<typename R, typename F>
void Genus<R,F>::threaded_compute_eigenvalues(std::promise<std::map<int64_t, std::map<R, int64_t>>>& aps)
{
    // Initialize running eigenvalue tallies.
    std::map<int64_t, std::map<R, int64_t>> values;
    for (auto& it : this->eigenvectorMap_)
    {
        for (Eigenvector& vec : it.second)
        {
            for (NeighborIterator<R,F>& neighborIt : this->eigenvalueNeighborVec_)
            {
                values[vec.index()][neighborIt.prime()] = 0;
            }
        }
    }

    // A NeighborIterator pointer.
    int64_t itIndex = 0;

    // The number of primes we're going to compute.
    int64_t numPrimes = this->eigenvalueNeighborVec_.size();

    // Initialize a vector to keep track of the boundaries for our neighbor
    // index.
    std::vector<int64_t> boundaries;
    boundaries.reserve(numPrimes+1);
    boundaries.push_back(0);

    // Set the boundaries.
    int64_t maxNeighbors = 0;
    for (NeighborIterator<R,F>& temp : this->eigenvalueNeighborVec_)
    {
        maxNeighbors += temp.num_neighbors();
        boundaries.push_back(maxNeighbors);
    }

    while (true)
    {
        // Get a neighbor index.
        int64_t nIndex = this->neighborIndex_++;

        // If we've exceeded the maximum number of neighbors, we're done.
        if (nIndex >= maxNeighbors)
        {
            aps.set_value(std::move(values));
            return;
        }

        // Update the NeighborIterator pointer.
        while (nIndex >= boundaries[itIndex+1])
        {
            ++itIndex;
        }

        // The prime and p-neighbor number that we're going to compute.
        const R& p = this->eigenvalueNeighborVec_[itIndex].prime();
        int64_t N = nIndex - boundaries[itIndex];

        // The source quadratic form.
        QuadFormPtr src = this->eigenvalueNeighborVec_[itIndex].quad_form();

        // Get the p-neighhor, reduce it, and identify its associated
        // genus rep.
        QuadFormPtr neighbor = this->eigenvalueNeighborVec_[itIndex].get_neighbor(N);
        QuadFormPtr reduced = QuadForm<R,F>::reduce(*neighbor);
        const GenusRep<R,F>& rep = this->find_genus_rep(reduced);

        // Determine whether we actually need to compute anything for this
        // p-neighbor.
        bool needed = false;
        for (auto& it : this->eigenvectorMap_)
        {
            const R& cond = it.first.conductor();
            int64_t relPos = rep.position(cond);
            if (relPos == -1) { continue; }

            for (auto& vec : it.second)
            {
                if (vec[relPos] == 0) { continue; }
                int64_t pivot = this->pivotMap_[vec.index()];
                if (*src == *this->genusVec_[pivot])
                {
                    needed = true;
                }
            }
        }

        // If we don't need to compute anything, get the next neighbor and
        // continue.
        if (!needed)
        {
            continue;
        }

        // Multiply the p-neighbor isometry on the right by the reduction
        // isometry. This now represents an isometry from the original
        // quadratic form to the genus representative isometric to this
        // p-neighbor.
        neighbor->isometry()->multiply_on_right_by(reduced->isometry());

#ifdef DEBUG
        assert( neighbor->isometry()->is_isometry(*this->q_, *rep.quad_form()) );
#endif

        // Obtain an automorphism of the original quadratic form.
        neighbor->isometry()->multiply_on_right_by(rep.inverse());

#ifdef DEBUG
        assert( neighbor->isometry()->is_automorphism(*this->q_) );
#endif

        // Convenient reference for the automorphism we'll use.
        std::shared_ptr<Isometry<R,F>> aut = neighbor->isometry();

        // Compute the character value for each of the primitive characters.
        std::map<R,int64_t> primeValues;
        for (const Character<R,F>& chi : this->primeCharSet_)
        {
            primeValues[chi.conductor()] = chi.rho(*aut, *this->q_);
        }

        for (auto& it : this->eigenvectorMap_)
        {
            // Get the character and the relative position of this genus rep
            // with respect to the character. If the p-neighbor we computed
            // doesn't contribute to this eigenvector, skip it.
            const Character<R,F>& chi = it.first;
            int64_t pos = rep.position(chi);
            if (pos == -1) { continue; }

            // Compute the value of this character.
            int64_t value = 1;
            for (const R& prime : chi.primes())
            {
                value *= primeValues[prime];
            }

            // Update local eigenvalues.
            for (Eigenvector& vec : it.second)
            {
                int64_t pivot = this->pivotMap_[vec.index()];
                QuadFormPtr q = this->genusVec_[pivot];
                if (*q == *src)
                {
                    values[vec.index()][p] += (vec[pos] * value);
                }
            }
        }
    }
}

template<typename R, typename F>
void Genus<R,F>::compute_eigenvalues(const R& p, int64_t numThreads)
{
    std::vector<R> ps = {p};
    this->compute_eigenvalues(ps, numThreads);
}

template<typename R, typename F>
void Genus<R,F>::compute_eigenvalues(const std::vector<R>& ps, int64_t numThreads)
{
    // Throw an error if there are no eigenvectors.
    if (this->eigenvectorMap_.size() == 0)
    {
        throw std::runtime_error("No eigenvectors imported.");
    }

    if (this->eigenvectorPivots_.empty())
    {

        // Assign eigenvector pivots. The goal of this method is to identify a
        // set of genus representatives for which we need to compute
        // p-neighbors in order to compute all eigenvalues for all
        // eigenvectors. This is a nontrivial problem and is equivalent to the
        // covering set problem which is NP-complete (thanks to Gonzalo
        // Tornaria for pointing this out). We have separated this out as its
        // own method so that alternative algorithms can be designed and
        // called based upon criteria such as the number of vectors, how many
        // coordinates they have, the number of eigenvalues to be computed,
        // their size, etc.
        this->assign_eigenvector_pivots_greedy();
    }

    // Get the bad primes.
    std::vector<R> badPrimes = std::move(this->bad_primes());

    // A vector of flags indicating whether we actually need to compute its
    // associated prime.
    std::vector<bool> needToCompute(ps.size(), false);

    // Initialize the eigenvalues; this will be updated as we go.
    for (auto& it : this->eigenvectorMap_)
    {
        const Character<R,F>& chi = it.first;
        const R& cond = chi.conductor();
        for (Eigenvector& vec : it.second)
        {
            // Create map of prime/eigenvalue pairs.
            if (this->eigenvalueMap_.count(vec.index()) == 0)
            {
                this->eigenvalueMap_[vec.index()] = std::move(std::map<R, int64_t>());
            }
            std::map<R, int64_t>& tempMap = this->eigenvalueMap_[vec.index()];

            // Initialize eigenvalues at the good primes, if necessary.
            int64_t index = 0;
            for (auto& p : ps)
            {
                if (tempMap.count(p) == 0)
                {
                    tempMap[p] = 0;
                    needToCompute[index] = true;
                }
                ++index;
            }

            // Initialize eigenvalues at the bad primes.
            for (auto& p : badPrimes)
            {
                if (cond % p == 0)
                {
                    tempMap[p] = -1;
                }
                else
                {
                    tempMap[p] = 1;
                }
            }
        }
    }

    if (numThreads > 0)
    {
        // Assign the shared neighbor iterator.
        int64_t index = 0;
        for (auto& p : ps)
        {
            if (needToCompute[index])
            {
                for (int64_t pivot : this->eigenvectorPivots_)
                {
                    QuadFormPtr q = this->genusVec_[pivot];
                    this->eigenvalueNeighborVec_.emplace_back(q, p);
                }
            }
            ++index;
        }

        // Initialize the neighbor index.
        this->neighborIndex_.store(0);

        // Create promises.
        std::vector<std::promise<std::map<int64_t, std::map<R, int64_t>>>> aps(numThreads);

        // Create the threads.
        std::vector<std::thread> threads;
        for (int64_t n = 0; n < numThreads; n++)
        {
            threads.push_back(std::thread(
                &Genus<R,F>::threaded_compute_eigenvalues, this, std::ref(aps[n])
            ));
        }

        // Wait for the threads to terminate.
        for (auto& t : threads)
        {
            t.join();
        }

        // Set eigenvalues.
        for (int64_t n = 0; n < numThreads; n++)
        {
            auto eigenvalues = aps[n].get_future().get();
            for (auto& it1 : eigenvalues)
            {
                int64_t index = it1.first;
                for (auto& it2 : it1.second)
                {
                    const R& p = it2.first;
                    this->eigenvalueMap_[index][p] += it2.second;
                }
            }
        }
    }
    else
    {
        int64_t index = 0;
        for (const R& p : ps)
        {
            // Skip this prime if we don't need to compute it for any of the
            // eigenvectors.
            if (!needToCompute[index++])
            {
                continue;
            }

            for (int64_t pivot : this->eigenvectorPivots_)
            {
                QuadFormPtr q = this->genusVec_[pivot];
                //const GenusRep<R,F>& g = this->find_genus_rep(q);

                // Build the neighbor iterator and let's get started.
                NeighborIterator<R,F> it(q, p);
                QuadFormPtr neighbor = it.next_neighbor();

                while (neighbor != nullptr)
                {
                    QuadFormPtr reduced = QuadForm<R,F>::reduce(*neighbor);
                    const GenusRep<R,F>& rep = this->find_genus_rep(reduced);

                    // Determine whether we actually need to compute anything for this
                    // p-neighbor.
                    bool needed = false;
                    for (auto& it : this->eigenvectorMap_)
                    {
                        const R& cond = it.first.conductor();
                        int64_t relPos = rep.position(cond);
                        if (relPos == -1) { continue; }

                        for (auto& vec : it.second)
                        {
                            if (vec[relPos] == 0) { continue; }
                            int64_t vecPivot = this->pivotMap_[vec.index()];
                            if (*q == *this->genusVec_[vecPivot])
                            {
                                needed = true;
                            }
                        }
                    }

                    // If we don't need to compute anything, get the next neighbor and
                    // continue.
                    if (!needed)
                    {
                        neighbor = it.next_neighbor();
                        continue;
                    }

                    // Multiply the p-neighbor isometry on the right by the reduction
                    // isometry. This now represents an isometry from the original
                    // quadratic form to the genus representative isometric to this
                    // p-neighbor.
                    neighbor->isometry()->multiply_on_right_by(reduced->isometry());

#ifdef DEBUG    
                    assert( neighbor->isometry()->is_isometry(*this->q_, *rep.quad_form()) );
#endif

                    // Obtain an automorphism of the original quadratic form.
                    neighbor->isometry()->multiply_on_right_by(rep.inverse());

#ifdef DEBUG    
                    assert( neighbor->isometry()->is_automorphism(*this->q_) );
#endif

                    // Convenient reference for the automorphism we'll use.
                    std::shared_ptr<Isometry<R,F>> aut = neighbor->isometry();

                    // Compute the character value for each of the primitive characters.
                    std::map<R,int64_t> primeValues;
                    for (const Character<R,F>& chi : this->primeCharSet_)
                    {
                        primeValues[chi.conductor()] = chi.rho(*aut, *this->q_);
                    }

                    for (auto& it : this->eigenvectorMap_)
                    {
                        // Get the character and the relative position of this genus rep
                        // with respect to the character. If the p-neighbor we computed
                        // doesn't contribute to this eigenvector, skip it.
                        const Character<R,F>& chi = it.first;
                        int64_t repPos = rep.position(chi);
                        if (repPos == -1) { continue; }

                        // Compute the value of this character.
                        int64_t value = 1;
                        for (const R& p : chi.primes())
                        {
                            value *= primeValues[p];
                        }

                        // Loop over the eigenvectors associated to this character.
                        for (auto& vec : it.second)
                        {
                            if (vec[repPos] == 0) { continue; }

                            int64_t pivot = this->pivotMap_[vec.index()];
                            if (*q == *this->genusVec_[pivot])
                            {
                                this->eigenvalueMap_[vec.index()][p] += (vec[repPos] * value);
                            }
                        }
                    }

                    neighbor = it.next_neighbor();
                }
            }
        }
    }

    // Recover the eigenvalues for each eigenvector.
    for (auto& it : this->eigenvectorMap_)
    {
        const Character<R,F>& chi = it.first;

        // Update the eigenvalues.
        for (auto& vec : it.second)
        {
            int64_t pivot = this->pivotMap_[vec.index()];
            QuadFormPtr q = this->genusVec_[pivot];
            const GenusRep<R,F>& g = this->find_genus_rep(q);
            int64_t relPos = g.position(chi);
            int64_t index = 0;
            for (const R& p : ps)
            {
                if (needToCompute[index])
                {
#ifdef DEBUG
                    assert( this->eigenvalueMap_[vec.index()][p] % vec[relPos] == 0 );
#endif
                    this->eigenvalueMap_[vec.index()][p] /= vec[relPos];
                }
                ++index;
            }
        }
    }
}

template<typename R, typename F>
int64_t Genus<R,F>::dimension(const Character<R,F>& chi) const
{
    const R& cond = chi.conductor();
    if (this->dimensionMap_.count(cond) > 0)
    {
        return this->dimensionMap_.find(cond)->second;
    }
    else
    {
        return 0;
    }
}

template<typename R, typename F>
int64_t Genus<R,F>::dimension(const R& cond) const
{
    if (this->dimensionMap_.count(cond) > 0)
    {
        return this->dimensionMap_.find(cond)->second;
    }
    else
    {
        return 0;
    }
}

template<typename R, typename F>
void Genus<R,F>::print(std::ostream& os) const
{
    os << this->genusVec_.size() << std::endl;
    for (QuadFormPtr q : this->genusVec_)
    {
        os << q << " ";
        os << q->isometry() << std::endl;
    }
    os << std::endl;

    for (auto& it1 : this->heckeMap_)
    {
        const R& cond = it1.first;
        const std::map<R, HeckePtr>& heckeMap = it1.second;
        for (auto& it2 : heckeMap)
        {
            const R& p = it2.first;
            HeckePtr hecke = it2.second;
            os << cond << " " << p << " ";
            os << *hecke;
            os << std::endl << std::endl;
        }
    }

    os << std::endl;
}

template<typename R, typename F>
void Genus<R,F>::import_genus(const std::string& filename)
{
    // Open file stream.
    std::ifstream infile(filename.c_str(), std::ios::in);

    // If the file isn't open... ¯\_(ツ)_/¯
    if (!infile.is_open())
    {
        throw std::runtime_error("Unable to open input file.");
    }

    int64_t dim = 0;
    infile >> dim;

    for (int64_t n = 0; n < dim; n++)
    {
        // Quadratic form coefficients.
        R a, b, c, f, g, h;
        infile >> a >> b >> c >> f >> g >> h;

        // Isometry coefficients.
        F a11, a12, a13, a21, a22, a23, a31, a32, a33;
        infile >> a11 >> a12 >> a13 >> a21 >> a22 >> a23 >> a31 >> a32 >> a33;

        // Create a quadratic form shared pointer, as well as the isometry
        // shared pointer, assign the isometry to the quadratic form, and add
        // it to the genus.
        QuadFormPtr qq = std::make_shared<QuadForm<R,F>>(this->disc_, a, b, c, f, g, h, true);
        auto s = std::make_shared<Isometry<R,F>>(a11, a12, a13, a21, a22, a23, a31, a32, a33);
        qq->isometry(s);

        // Set the initial quadratic form.
        if (n == 0)
        {
            this->q_ = qq;
            this->disc_ = qq->discriminant();

            if (!Math<R,F>::is_squarefree(this->disc_))
            {
                throw std::runtime_error("Discriminant must be squarefree.");
            }

            if (!Math<R,F>::is_positive(this->disc_))
            {
                throw std::runtime_error("Discriminant must be positive.");
            }
        }

        this->add_genus_rep(qq);
    }

    // Read Hecke operators from file.
    R cond, p;
    infile >> cond >> p;

    // Create character and add it to the genus.
    Character<R,F> chi(cond);
    this->add_character(chi);

    // Create map for this conductor, if one doesn't already exist.
    if (this->heckeMap_.count(cond) == 0)
    {
        this->heckeMap_[cond] = std::move(std::map<R, HeckePtr>());
    }

    // Continue reading file until we reach the end.
    while (!infile.eof())
    {
        // Create shared pointer to Hecke operator.
        std::shared_ptr<HeckeOperator<R,F>> hecke = std::make_shared<HeckeOperator<R,F>>(*this, chi);

        // Import values into Hecke operator.
        hecke->import(infile);

        // Assign Hecke operator.
        this->heckeMap_[cond][p] = hecke;

        // Attempt to read conductor and prime from file.
        infile >> cond >> p;
    }
}

template<typename R, typename F>
void Genus<R,F>::import_eigenvectors(const std::string& filename)
{
    // Open a file stream.
    std::ifstream infile(filename.c_str(), std::ios::in);

    // Can't open the file? No problem. Throw an exception!
    if (!infile.is_open())
    {
        throw std::runtime_error("Unable to open input file.");
    }

    // Determine how many conductors we have, and build the associated characters.
    int64_t numConds;
    infile >> numConds;

    // Continue reading until end of file.
    for (int64_t n = 0; n < numConds; n++)
    {
        R cond;
        infile >> cond;

        // Get the primes dividing this conductor.
        std::vector<R> ps = Math<R,F>::prime_divisors_naive(cond);

        // Construct the character.
        Character<R,F> chi(ps);
        this->add_character(chi);
    }

    // Compute the genus.
    this->compute_genus(this->genusNumThreads);

    // Read eigenvectors from file.
    size_t numEigenvectors = 0;

    R cond;
    infile >> cond;

    while (!infile.eof())
    {
        // Do we have a character with this conductor?
        if (this->dimensionMap_.count(cond) == 0)
        {
            throw std::runtime_error("Eigenvector with unassociated conductor found, cannot proceed.");
        }

        // Okay, let's get it, and figure out the dimension.
        const Character<R,F>& chi = this->condToChar_[cond];
        int64_t dim = this->dimensionMap_[cond];

        // Read the eigenvector coefficients.
        Eigenvector vec(dim, numEigenvectors);
        infile >> vec;

        // Set up data structure for storing eigenvectors.
        if (this->eigenvectorMap_.count(chi) == 0)
        {
            this->eigenvectorMap_[chi] = std::move(std::vector<Eigenvector>());
        }
        this->eigenvectorMap_[chi].push_back(std::move(vec));

        // Set up data structure for storing eigenvalues.
        this->eigenvalueMap_[numEigenvectors] = std::map<R, int64_t>();

        // Read pre-computed eigenvalues.
        int64_t numPrimes;
        infile >> numPrimes;
        for (int64_t n = 0; n < numPrimes; n++)
        {
            R p;
            int64_t ap;
            infile >> p >> ap;
            this->eigenvalueMap_[numEigenvectors][p] = ap;
        }

        ++numEigenvectors;

        // Read the next conductor, or get eof flag.
        infile >> cond;
    }

    // Close file stream.
    infile.close();

    // Populate the absolute position vectors with zeroes.
    for (auto& it : this->eigenvectorMap_)
    {
        const R& cond = it.first.conductor();
        this->absolutePosition_[cond] = std::vector<int64_t>(this->dimension(cond), -1);
    }

    // Populate the absolute position.
    for (auto& rep : this->genusReps_)
    {
        int64_t absPos = rep.absolute_position();
        for (auto& it : this->eigenvectorMap_)
        {
            const R& cond = it.first.conductor();
            int64_t relPos = rep.position(cond);

            if (relPos != -1)
            {
                this->absolutePosition_[cond][relPos] = absPos;
            }
        }
    }
}

template<typename R, typename F>
void Genus<R,F>::assign_eigenvector_pivots_greedy(void)
{
    // A vector which will be populated with eigenvector pivots. These
    // correspond to the row of the Hecke operator which needs to be computed
    // in order to generate all eigenvalue for all eigenvectors.
    this->eigenvectorPivots_.clear();

    // An unordered map of the eigenvector indices which are known to be
    // accounted for by the eigenvectorPivots_ vector.
    std::unordered_set<int64_t> accountedFor;

    // The total number of eigenvectors.
    size_t numEigenvectors = this->eigenvalueMap_.size();

    // The number of eigenvectors we hope to account for during this iteration.
    int64_t goal = numEigenvectors;

    while (accountedFor.size() < numEigenvectors)
    {
        // A vector containing the number of nonzero entries at each absolute
        // position across all eigenvectors. This will be used to identify which
        // genus rep(s) will be used to compute eigenvalues.
        std::vector<int64_t> nonzero(this->genusVec_.size(), 0);

        for (auto& it : this->eigenvectorMap_)
        {
            const R& cond = it.first.conductor();
            const std::vector<int64_t>& absPos = this->absolutePosition_[cond];
            for (auto& vec : it.second)
            {
                if (accountedFor.count(vec.index()) > 0) { continue; }

                const auto& coeffs = vec.coefficients();
                int64_t len = coeffs.size();
                for (int64_t k = 0; k < len; k++)
                {
                    if (coeffs[k] != 0)
                    {
                        ++nonzero[absPos[k]];
                        if (nonzero[absPos[k]] == goal)
                        {
                            // Push the pivot onto the vector.
                            this->eigenvectorPivots_.push_back(absPos[k]);

                            // Assign this pivot to all remaining vectors.
                            for (size_t n = 0; n < numEigenvectors; n++)
                            {
                                if (accountedFor.count(n) == 0)
                                {
                                    this->pivotMap_[n] = absPos[k];
                                }
                            }
                            return;
                        }
                    }
                }
            }
        }

        // Determine the best (greedy) pivot position to include.
        auto maxelt = std::max_element(nonzero.begin(), nonzero.end());
        int64_t bestIndex = std::distance(nonzero.begin(), maxelt);
        QuadFormPtr qq = this->genusVec_[bestIndex];
        const GenusRep<R,F>& rep = this->find_genus_rep(qq);

        for (auto& it : this->eigenvectorMap_)
        {
            const R& cond = it.first.conductor();

            int64_t relPos = rep.position(cond);

            for (auto& vec : it.second)
            {
                // Skip this vector if it has already been assigned a pivot.
                if (accountedFor.count(vec.index()) > 0)
                {
                    continue;
                }

                // Skip this vector if the character does not allow it to be
                // accounted for.
                if (relPos == -1)
                {
                    continue;
                }

                // If the vector has a nonzero value in the best position, add
                // it to the accounted for set.
                if (vec[relPos] != 0)
                {
                    accountedFor.insert(vec.index());
                    this->pivotMap_[vec.index()] = bestIndex;
                }
            }
        }

        // Add this pivot to the list.
        this->eigenvectorPivots_.push_back(bestIndex);

        // Update the goal.
        goal = numEigenvectors - accountedFor.size();
    }
}

template<typename R, typename F>
void Genus<R,F>::export_eigenvectors(const std::string& filename)
{
    std::ofstream outfile(filename.c_str(), std::ios::out);

    if (!outfile.is_open())
    {
        throw std::runtime_error("Unable to open output file for writing.");
    }

    outfile << this->eigenvectorMap_.size();
    for (auto& it1 : this->eigenvectorMap_)
    {
        outfile << " " << it1.first.conductor();
    }
    outfile << std::endl << std::endl;

    for (auto& it1 : this->eigenvectorMap_)
    {
        const Character<R,F>& chi = it1.first;
        std::vector<Eigenvector>& list = it1.second;
        for (Eigenvector& vec : list)
        {
            outfile << chi.conductor() << " " << vec << std::endl;

            auto& eigs = this->eigenvalueMap_[vec.index()];
            outfile << eigs.size();
            for (auto& it2 : eigs)
            {
                outfile << " " << it2.first << " " << it2.second;
            }
            outfile << std::endl << std::endl;
        }
    }

    outfile.close();
}

template<typename R, typename F>
std::ostream& operator<<(std::ostream& os, const Genus<R,F>& genus)
{
    genus.print(os);
    return os;
}

#endif // __GENUS_H_
