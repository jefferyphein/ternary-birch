#ifndef __GENUS_H_
#define __GENUS_H_

#include <memory>
#include <cassert>
#include <unordered_set>
#include <set>
#include <map>
#include <queue>
#include <thread>
#include <mutex>
#include <chrono>
#include "QuadForm.h"
#include "NeighborIterator.h"
#include "Character.h"
#include "HeckeOperator.h"

template<typename R, typename F>
class GenusRep
{
public:
    typedef std::shared_ptr<QuadForm<R,F>> QuadFormPtr;
    typedef std::shared_ptr<Isometry<R,F>> IsometryPtr;

    GenusRep(QuadFormPtr qPtr);

    int64_t pos(const Character<R,F>& chi) const;
    void set_position(const Character<R,F>& chi, int64_t x);
    int64_t position(const Character<R,F>& chi) const;
    void set_inverse(void);
    QuadFormPtr quad_form(void) const;
    IsometryPtr inverse(void) const;

    void set_dimension(const Character<R,F>& chi, int64_t dim);
    int64_t dimension(const Character<R,F>& chi) const;

    bool operator==(const GenusRep<R,F>& rep) const;
    bool operator<(const GenusRep<R,F>& rep) const;
private:
    QuadFormPtr q_;
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

    Genus(const QuadForm<R,F>& q);

    QuadFormPtr quad_form(void) const;

    R smallest_good_prime(void) const;

    void add_character(const Character<R,F>& chi);

    void set_dimensions(GenusRep<R,F>& rep);

    void print(void) const;

    size_t size(void) const;

    int64_t dimension(const Character<R,F>& chi) const;

    void compute_hecke_operators(const R& p, int64_t numThreads=0);

    HeckePtr hecke_operator(const R& p, const Character<R,F>& chi);

    /* Computes the full genus. */
    // TODO: MAKE THIS PRIVATE. THIS IS PUBLIC FOR TESTING PURPSOES ONLY.
    void compute_genus(int64_t numThreads=0);

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
    void add_genus_rep(QuadFormPtr q, QuadFormPtr neighbor=nullptr);

    /* Finds a genus representative. Will throw an exception if not found. */
    const GenusRep<R,F>& find_genus_rep(QuadFormPtr q) const;

    ////////// THREADING RESOURCES

    /* This function is passed to threads, pulls a genus representative off of
     * a queue and computes its p-neighbors. */
    void threaded_compute_genus(const R& p, int64_t threadid);

    /* This function is passed to threads, pulls a genus representative out of
     * a vector based on the value of genusVecIndex_. */
    void threaded_compute_hecke_operators(const R& p, int64_t threadid);

    /* A vector of booleans used to determine whether a thread is blocked due
     * to there being no objects in the genusQueue_. */
    std::vector<bool> threadBlocked_;

    /* A queue of shared QuadForm pointers used for multithreading. */
    std::queue<QuadFormPtr> genusQueue_;

    /* An index used to access GenusReps when computing Hecke operators. */
    int64_t genusVecIndex_;

    /* A mutex used for accessing genus reps when multithreading. */
    std::mutex genusMutex_;

    /* A mutex used for updating Hecke operator values when multithreading. */
    std::mutex heckeMutex_;

    ////////// GENUS DATA STRUCTURES

    /* An unordered set which stores distinct genus representatives. */
    std::unordered_set<GenusRep<R,F>> genusReps_;

    /* A vector of shared QuadForm pointers. An unordered list of the genus. */
    std::vector<QuadFormPtr> genusVec_;

    /* A set of characters to compute. */
    std::set<Character<R,F>> charSet_;

    /* A set of prime characters to compute. When representation
     * values are compuated for characters in charSet_, these will first be
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
int64_t GenusRep<R,F>::position(const Character<R,F>& chi) const
{
    return this->positionMap_.find(chi.conductor())->second;
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
    if (this->dimensionMap_.count(chi.conductor()) > 0)
    {
        return this->dimensionMap_.find(chi.conductor())->second;
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
Genus<R,F>::Genus(const QuadForm<R,F>& q)
{
    this->q_ = QuadForm<R,F>::reduce(q, false);
    this->disc_ = this->q_->discriminant();
}

template<typename R, typename F>
std::shared_ptr<QuadForm<R,F>> Genus<R,F>::quad_form(void) const
{
    return this->q_;
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
    if (this->computed_) { return; }

    // Initialize the dimensions for each character to zero.
    for (auto& chi : this->charSet_)
    {
        this->dimensionMap_[chi.conductor()] = 0;
    }

    // Add the defining quadratic form as the initial genus rep.
    this->add_genus_rep(this->q_);

    // Determine the smallest good prime.
    R p = this->smallest_good_prime();

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
        int64_t index = 0;
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
}

template<typename R, typename F>
const GenusRep<R,F>& Genus<R,F>::find_genus_rep(QuadFormPtr q) const
{
    GenusRep<R,F> rep(q);
    return *(this->genusReps_.find(rep));
}

template<typename R, typename F>
void Genus<R,F>::add_genus_rep(QuadFormPtr q, QuadFormPtr neighbor)
{
    // Create a GenusRep object from the shared pointer.
    GenusRep<R,F> rep(q);

    // Is this quadratic form already in the genus? If not, add it.
    if (this->genusReps_.count(rep) == 0)
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
            // Update the dimensions for each character.
            this->dimensionMap_[chi.conductor()] += rep.dimension(chi);

            // Set the relative position for this genus rep with respect to
            // each character. A position of -1 indicates that this genus rep
            // does not contribute to the dimension of the space associated to
            // this character.
            if (rep.dimension(chi) > 0)
            {
                rep.set_position(chi, this->dimensionMap_[chi.conductor()]-1);
            }
            else
            {
                rep.set_position(chi, -1);
            }
        }

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
    // If this is a new character, add it to the set.
    if (this->charSet_.count(chi) == 0)
    {
        this->charSet_.insert(chi);

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
        this->genusMutex_.lock();
        int64_t index = this->genusVecIndex_++;
        this->genusMutex_.unlock();

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
            int64_t pos = rep.position(chi);
            if (pos != -1)
            {
                this->heckeMutex_.lock();
                this->heckeMap_[p][chi.conductor()]->update_row(pos, rowMap[chi.conductor()]);
                this->heckeMutex_.unlock();
            }
        }
    }
}

template<typename R, typename F>
void Genus<R,F>::compute_hecke_operators(const R& p, int64_t numThreads)
{
    if (this->heckeMap_.count(p) > 0)
    {
        const std::map<R, HeckePtr> hecke = this->heckeMap_.find(p)->second;

        std::cout << "Hecke operators at " << p << " already computed." << std::endl;
        return;
    }

    // Build a Hecke operator map based on character conductors.
    this->heckeMap_[p] = std::map<R, HeckePtr>();
    for (auto& chi : this->charSet_)
    {
        this->heckeMap_[p][chi.conductor()] = std::make_shared<HeckeOperator<R,F>>(*this, chi);
    }

    if (numThreads > 0)
    {
        // Initialize the shared pointer used by threads to compute the Hecke
        // operators.
        this->genusVecIndex_ = 0;

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
                this->heckeMap_[p][cond]->add_row(rep.position(chi), rowMap[cond]);
            }
        }
    }

    //auto h = this->heckeMap_.find(p)->second;
    //for (auto& chi : this->charSet_)
    //{
    //    auto mat = h.find(chi.conductor())->second;
    //    std::cout << chi.conductor() << " " << p << " ";
    //    mat->print();
    //    std::cout << std::endl;
    //}
}

template<typename R, typename F>
int64_t Genus<R,F>::dimension(const Character<R,F>& chi) const
{
    if (this->dimensionMap_.count(chi.conductor()) > 0)
    {
        return this->dimensionMap_.find(chi.conductor())->second;
    }
    else
    {
        return 0;
    }
}

template<typename R, typename F>
void Genus<R,F>::print(void) const
{
    //std::cout << " conductor   dimension" << std::endl;
    //for (auto& chi : this->charSet_)
    //{
    //    std::cout << std::setw(10) << chi.conductor() << "   ";
    //    std::cout << this->dimension(chi) << std::endl;
    //}
    //for (auto& q : this->genusVec_)
    //{
    //    std::cout << q << std::endl;
    //}
}

#endif // __GENUS_H_
