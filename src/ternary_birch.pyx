# distutils: language = c++
# distutils: sources = birch_util.cpp Fp.cpp Isometry.cpp Math.cpp QuadForm.cpp SetCover.cpp
# distutils: extra_compile_args = -g -O3 -Wall -Werror -std=c++11 -fvar-tracking-assignments-toggle

from __future__ import print_function

import logging
import numpy as np

from queue import Queue
from datetime import datetime
from scipy.sparse import csr_matrix

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map as cppmap
from libcpp.memory cimport shared_ptr, make_shared
from libc.stdint cimport int32_t as Z32
from libc.stdint cimport int64_t as Z64
from libc.stdint cimport uint16_t as W16
from libc.stdint cimport uint32_t as W32
from libc.stdint cimport uint64_t as W64
from libc.math cimport sqrt, floor

from operator import itemgetter
from random import randint

from cython.operator import dereference as deref
from cython.operator import postincrement as incr

from sage.rings.integer cimport Integer
from sage.libs.gmp.types cimport mpz_t

from sage.all import QuadraticForm
from sage.all import Integers, Rationals
from sage.all import random_prime
from sage.all import matrix
from sage.all import next_prime
from sage.all import GF

from math import log, exp, lgamma

HASSE_MULTIPLIER = 50

cdef extern from "<utility>" namespace "std" nogil:
    T move[T](T)

cdef extern from "gmpxx.h":
    cdef cppclass mpz_class:
        mpz_class(mpz_t a)
        string get_str(int base)

cdef extern from "QuadForm.h":
    cdef cppclass PrimeSymbol[R]:
        R p
        int power
        int ramified

    cdef cppclass QuadForm[R]:
        QuadForm()
        const R& a() const
        const R& b() const
        const R& c() const
        const R& f() const
        const R& g() const
        const R& h() const

        @staticmethod
        QuadForm[R] get_quad_form(const vector[PrimeSymbol[R]]& primes) except +

cdef extern from "Eigenvector.h":
    cdef cppclass Eigenvector[R]:
        const vector[Z32]& data() const

    cdef cppclass EigenvectorManager[R]:
        void add_eigenvector(Eigenvector[R]&& vector)
        size_t size() const
        const Eigenvector[R]& operator[](size_t index) const
        void finalize()

cdef extern from "Genus.h":
    cdef cppclass Genus[R]:
        Genus()
        Genus(const QuadForm[R]& q, const vector[PrimeSymbol[R]]& symbols, W64 seed)
        cppmap[R,size_t] dimension_map() const
        W64 seed() const
        cppmap[R,vector[int]] hecke_matrix_dense(const R& p) except +
        cppmap[R,vector[vector[int]]] hecke_matrix_sparse(const R& p) except +

        Eigenvector[R] eigenvector(const vector[Z32]& vec, const R& conductor) except +
        vector[Z32] eigenvalues(EigenvectorManager[R]& manager, const R& p) except +

        @staticmethod
        Genus[T] convert[T](const Genus[R]& src)

cdef extern from "Isometry.h":
    cdef cppclass Isometry[R]:
        R a11
        R a12
        R a13
        R a21
        R a22
        R a23
        R a31
        R a32
        R a33
        pass

cdef extern from "IsometrySequence.h":
    cdef cppclass IsometrySequenceData[T]:
        Isometry[T] isometry
        T denominator
        size_t src
        size_t dst
    cdef cppclass IsometrySequence[R,S,T]:
        IsometrySequence(shared_ptr[Genus[T]], const T& p)
        int done() const
        IsometrySequenceData next()

ctypedef mpz_class Z
ctypedef PrimeSymbol[Z] Z_PrimeSymbol
ctypedef QuadForm[Z] Z_QuadForm

cdef class BirchGenus:
    cdef shared_ptr[Genus[Z]] Z_genus
    cdef shared_ptr[Genus[Z64]] Z64_genus
    cdef EigenvectorManager[Z] Z_manager
    cdef EigenvectorManager[Z64] Z64_manager
    cdef Z64_genus_is_set
    cdef level_
    cdef ramified_primes_
    cdef facs
    cdef dims
    cdef seed_
    cdef hecke
    cdef sage_hecke
    cdef eigenvectors

    def __init__(self, level, ramified_primes=None, seed=None):
        """
        EXAMPLES::

            sage: from ternary_birch import BirchGenus
            sage: g = BirchGenus(11*13*17*19*23)
        """
        self.level_ = Integer(level)
        self.facs = self.level_.factor()
        ps = list(map(itemgetter(0), self.facs))
        es = list(map(itemgetter(1), self.facs))

        if ramified_primes is None:
            logging.info("Ramified primes: chosen automatically")
            if len(ps) % 2 == 0:
                self.ramified_primes_ = ps[:-1]
            else:
                self.ramified_primes_ = ps[:]
        else:
            logging.info("Ramified primes: specified by user")
            self.ramified_primes_ = [ p for p in ramified_primes if p in ps ]

        cdef vector[Z_PrimeSymbol] primes
        cdef Z_PrimeSymbol prime
        for n,p in enumerate(ps):
            prime.p = Z(Integer(p).value)
            prime.power = int(es[n])
            prime.ramified = p in self.ramified_primes_
            primes.push_back(prime)
            logging.info("%s at %s", "Ramified" if prime.ramified else "Unramified", p)

        cdef Z_QuadForm q
        try:
            logging.info("Determining desired quadratic form")
            q = Z_QuadForm.get_quad_form(primes)
            a = _Z_to_int(q.a())
            b = _Z_to_int(q.b())
            c = _Z_to_int(q.c())
            f = _Z_to_int(q.f())
            g = _Z_to_int(q.g())
            h = _Z_to_int(q.h())
            S = "Chose Q(x,y,z) = "
            if a:
                S += "{}{}x^2".format("" if a > 0 else "-", abs(a) if abs(a) != 1 else "")
            if b:
                S += " {} {}y^2".format("+" if b > 0 else "-", abs(b) if abs(b) != 1 else "")
            if c:
                S += " {} {}z^2".format("+" if c > 0 else "-", abs(c) if abs(c) != 1 else "")
            if f:
                S += " {} {}yz".format("+" if f > 0 else "-", abs(f) if abs(f) != 1 else "")
            if g:
                S += " {} {}xz".format("+" if g > 0 else "-", abs(g) if abs(g) != 1 else "")
            if h:
                S += " {} {}xy".format("+" if h > 0 else "-", abs(h) if abs(h) != 1 else "")
            logging.info(S)
        except Exception as e:
            raise Exception(e.message)

        seed = seed if seed else 0
        cdef W64 arg_seed = seed

        logging.info("Computing genus representatives...")
        genus_start = datetime.now()
        self.Z_genus = make_shared[Genus[Z]](q, primes, arg_seed)
        genus_stop = datetime.now()
        logging.info("Finished computing genus representatives (time: %s)", genus_stop-genus_start)
        self.seed_ = deref(self.Z_genus).seed()
        logging.info("Seed = %s (%s)", self.seed_, "provided by user" if seed else "set randomly")

        cdef cppmap[Z,size_t] mymap = deref(self.Z_genus).dimension_map()
        cdef cppmap[Z,size_t].iterator it = mymap.begin()
        self.dims = dict()
        while it != mymap.end():
            self.dims[_Z_to_int(deref(it).first)] = deref(it).second
            incr(it)

        self.Z64_genus_is_set = False

        self.hecke = dict()
        self.sage_hecke = dict()
        self.eigenvectors = None

    def dimensions(self):
        return self.dims

    def dimension(self):
        total = 0
        for x in self.dims:
            total += self.dims[x]
        return total

    def level(self):
        return self.level_

    def seed(self):
        return self.seed_

    def ramified_primes(self):
        return self.ramified_primes_

    def next_good_prime(self, p):
        while True:
            p = next_prime(p)
            if self.level_ % p != 0:
                break
        return p

    def rational_eigenvectors(self, precise=True, sparse=None, force=False):
        # Do not duplicate work if we've already computed the eigenvectors
        if not force and self.eigenvectors is not None:
            return self.eigenvectors

        self.eigenvectors = []

        # Find the smallest good prime.
        p = self.next_good_prime(1)

        # The Hasse bound for our starting prime.
        hasse = int(floor(2*sqrt(1.0*p)))

        # Pick a prime considerably larger than the absolute value of the Hasse
        # bound. We do this to avoid exhausting over all possible eigenvalues,
        # since this can get fairly costly for large primes p. Instead, we
        # look for candidate eigenvalues by computing the roots of the
        # characteristic polynomial of our Hecke matrix over some larger
        # finite field q. By choosing a larger finite field, we reduce the
        # probability of picking up "incidental" eigenvalues which do not
        # actually correspond to rational eigenvectors.
        q = next_prime(HASSE_MULTIPLIER * hasse)

        # Produce some initial jobs to get the ball rolling.
        job_queue = Queue()
        for conductor in self.dims:
            A = self.sage_hecke_matrix(p, conductor, precise=precise, sparse=sparse)

            logging.info("Looking for eigenvalues (p=%s, conductor=%s) over GF(%s)...", p, conductor, q)
            start_time = datetime.now()

            # Determine all possible eigenvalues within the Hasse bound.
            roots = A.change_ring(GF(q)).characteristic_polynomial().roots()
            roots = [ Integers()(pair[0]) for pair in roots ]
            roots = [ rt-q if rt > hasse else rt for rt in roots ]
            roots = [ rt for rt in roots if abs(rt) <= hasse ]
            end_time = datetime.now()
            logging.info("  found %s possible eigenvalue(s) (time: %s)", len(roots), end_time-start_time)

            for e in roots:
                job = dict()
                job['p'] = p
                job['e'] = e
                job['aps'] = dict()
                job['matrix'] = A
                job['conductor'] = conductor
                job['subspace'] = None
                job_queue.put(job)

        # Process kernels as long as there are jobs in the queue.
        while job_queue.qsize() > 0:
            job = job_queue.get()

            # Unpack the job.
            p = job['p']
            e = job['e']
            matrix = job['matrix']
            conductor = job['conductor']
            subspace = job['subspace']

            logging.info("Computing kernel of %s-matrix (p=%s, conductor=%s, eigenvalue=%s)...", matrix.dimensions(), p, conductor, e)

            # Find the rational eigenspace with the specified value.
            start_time = datetime.now()
            if subspace is None:
                kernel = (matrix - e).right_kernel()
            else:
                kernel = (matrix - e*subspace.transpose()).right_kernel()
            end_time = datetime.now()

            dim = kernel.dimension()

            # Nothing to do if there is no eigenspace.
            if dim == 0:
                logging.info("  trivial kernel found (time: %s)", end_time-start_time)
                continue

            # Assign the eigenvalue to a new dictionary.
            aps = dict(job['aps'])
            aps[p] = e

            # Get the kernel as a matrix object and return it to the original
            # coordinates.
            S = kernel.matrix()
            if subspace is not None:
                S = S * subspace

            if dim == 1:
                logging.info("  eigenvector found (time: %s)", end_time-start_time)

                x = dict()
                x['vector'] = S[0]
                x['aps'] = aps
                x['conductor'] = conductor
                self.eigenvectors.append(x)
                continue

            logging.info("  %s-dimensional eigenspace found (time: %s)", S.dimensions()[0], end_time-start_time)

            # If we found a multi-dimensional eigenspace, we need to break it
            # apart using the Hecke matrix at the next good prime.
            p = self.next_good_prime(p)

            hasse = int(floor(2*sqrt(1.0*p)))
            q = next_prime(HASSE_MULTIPLIER*hasse)

            # Project the Hecke matrix onto the desired eigenspace.
            A = self.sage_hecke_matrix(p, conductor, precise=precise, sparse=sparse)

            logging.info("Looking for eigenvalues (p=%s, conductor=%s) over GF(%s)...", p, conductor, q)
            start_time = datetime.now()
            roots = A.change_ring(GF(q)).characteristic_polynomial().roots()
            roots = [ Integers()(pair[0]) for pair in roots ]
            roots = [ rt-q if rt > hasse else rt for rt in roots ]
            roots = [ rt for rt in roots if abs(rt) <= hasse ]
            end_time = datetime.now()
            logging.info("  found %s possible eigenvalue(s) (time: %s)", len(roots), end_time-start_time)

            A = A * S.transpose()

            for e in roots:
                newjob = dict()
                newjob['p'] = p
                newjob['e'] = e
                newjob['aps'] = aps
                newjob['matrix'] = A
                newjob['conductor'] = conductor
                newjob['subspace'] = S
                job_queue.put(newjob)

        # Set up the eigenvector manager so we can compute eigenvalues.
        self.reset_eigenvector_manager()

        return self.eigenvectors

    # TODO: Make it so that computing neighbors with 32-bit primes doesn't
    # crash when using precise=False.

    def add_eigenvector(self, eigenvector):
        if self.eigenvectors is None:
            self.eigenvectors = []

        self.eigenvectors.append(dict(eigenvector))

    def compute_eigenvalues(self, p, precise=True):
        prime = Integer(p)

        if not prime.is_prime():
            raise Exception("p is not prime.")

        if self.level_ % prime == 0:
            raise Exception("Cannot compute eigenvalues at primes dividing the level.")

        # If we haven't already computed the eigenvectors, do so now.
        if self.eigenvectors is None:
            self.rational_eigenvectors()

        # TODO: Make this better.
        cdef vector[Z32] data
        if not precise and not self.Z64_genus_is_set:
            self.Z64_genus = make_shared[Genus[Z64]](deref(self.Z_genus))
            self.Z64_genus_is_set = True

            for entry in self.eigenvectors:
                vec = entry['vector']
                cond = entry['conductor']
                dimension = len(vec)
                data = vector[Z32](dimension)
                for n,value in enumerate(vec):
                    data[n] = value
                self.Z64_manager.add_eigenvector(deref(self.Z64_genus).eigenvector(data, Integer(cond)))
            self.Z64_manager.finalize()

        # Compute the eigenvalues.
        cdef vector[Z32] aps
        if precise:
            aps = deref(self.Z_genus).eigenvalues(self.Z_manager, Z(Integer(p).value))
        else:
            aps = deref(self.Z64_genus).eigenvalues(self.Z64_manager, Integer(p))

        # Store eigenvalues in memory.
        for n,vec in enumerate(self.eigenvectors):
            vec['aps'][prime] = aps[n]

        return [ aps[n] for n in range(aps.size()) ]

    def reset_eigenvector_manager(self):
        cdef EigenvectorManager[Z] _Z_manager
        cdef EigenvectorManager[Z64] _Z64_manager

        if not self.Z64_genus_is_set:
            self.Z64_genus = make_shared[Genus[Z64]](deref(self.Z_genus))
            self.Z64_genus_is_set = True

        for entry in self.eigenvectors:
            vec = entry['vector']
            cond = entry['conductor']
            dimension = len(vec)
            data = vector[Z32](dimension)
            for n,value in enumerate(vec):
                data[n] = value
            _Z64_manager.add_eigenvector(deref(self.Z64_genus).eigenvector(data, Integer(cond)))
            _Z_manager.add_eigenvector(deref(self.Z_genus).eigenvector(data, Z(Integer(cond).value)))

        _Z64_manager.finalize()
        _Z_manager.finalize()

        self.Z64_manager = _Z64_manager
        self.Z_manager = _Z_manager

    def compute_eigenvalues_upto(self, upper, precise=True, force=False):
        ps = []
        p = 1
        while True:
            p = self.next_good_prime(p)
            if p <= upper:
                ps.append(p)
            else:
                break

        # If eigenvectors haven't already been computed, do so now.
        if self.eigenvectors is None:
            self.rational_eigenvectors()

        # TODO: Make this better.
        cdef vector[Z32] data
        if not precise and not self.Z64_genus_is_set:
            self.Z64_genus = make_shared[Genus[Z64]](deref(self.Z_genus))
            self.Z64_genus_is_set = True

            for entry in self.eigenvectors:
                vec = entry['vector']
                cond = entry['conductor']
                dimension = len(vec)
                data = vector[Z32](dimension)
                for n,value in enumerate(vec):
                    data[n] = value
                self.Z64_manager.add_eigenvector(deref(self.Z64_genus).eigenvector(data, Integer(cond)))
            self.Z64_manager.finalize()

        cdef vector[Z32] aps
        for p in ps:
            # Skip this prime if all eigenvectors already has its eigenvalue.
            if not force and all(p in vec['aps'] for vec in self.eigenvectors):
                continue

            if precise:
                aps = deref(self.Z_genus).eigenvalues(self.Z_manager, Z(Integer(p).value))
            else:
                aps = deref(self.Z64_genus).eigenvalues(self.Z64_manager, Integer(p))

            for n,vec in enumerate(self.eigenvectors):
                vec['aps'][p] = aps[n]

    def isometry_sequence(self, p, precise=True):
        prime = Integer(p)

        if not prime.is_prime():
            raise Exception("p is not prime.")

        if self.level_ % p == 0:
            raise Exception("Cannot compute Hecke matrix at primes dividing the level.")

        if not precise:
            if not self.Z64_genus_is_set:
                logging.info("Converting arbitrary precision Genus object to fixed-precision Genus object...")
                self.Z64_genus = make_shared[Genus[Z64]](deref(self.Z_genus))
                self.Z64_genus_is_set = True

            return self._isometry_sequence_imprecise(prime)
        else:
            return self._isometry_sequence_precise(prime)

    def _isometry_sequence_imprecise(self, Integer p):
        cdef shared_ptr[IsometrySequence[W16,W32,Z64]] sequence
        cdef Z64 prime = p
        sequence = make_shared[IsometrySequence[W16,W32,Z64]](self.Z64_genus, prime)

        cdef IsometrySequenceData[Z64] data
        while not deref(sequence).done():
            data = deref(sequence).next()

            a11 = Integer(data.isometry.a11)
            a12 = Integer(data.isometry.a12)
            a13 = Integer(data.isometry.a13)
            a21 = Integer(data.isometry.a21)
            a22 = Integer(data.isometry.a22)
            a23 = Integer(data.isometry.a23)
            a31 = Integer(data.isometry.a31)
            a32 = Integer(data.isometry.a32)
            a33 = Integer(data.isometry.a33)
            s = matrix(Integers(), 3, [ a11, a12, a13, a21, a22, a23, a31, a32, a33 ])

            retval = dict()
            retval['isometry'] = s
            retval['denominator'] = Integer(data.denominator)
            retval['src'] = Integer(data.src)
            retval['dst'] = Integer(data.dst)
            yield retval

    def _isometry_sequence_precise(self, Integer p):
        cdef shared_ptr[IsometrySequence[W16,W32,Z]] sequence
        sequence = make_shared[IsometrySequence[W16,W32,Z]](self.Z_genus, Z(p.value))

        cdef IsometrySequenceData[Z] data
        while not deref(sequence).done():
            data = deref(sequence).next()
            a11 = _Z_to_int(data.isometry.a11)
            a12 = _Z_to_int(data.isometry.a12)
            a13 = _Z_to_int(data.isometry.a13)
            a21 = _Z_to_int(data.isometry.a21)
            a22 = _Z_to_int(data.isometry.a22)
            a23 = _Z_to_int(data.isometry.a23)
            a31 = _Z_to_int(data.isometry.a31)
            a32 = _Z_to_int(data.isometry.a32)
            a33 = _Z_to_int(data.isometry.a33)
            s = matrix(Integers(), 3, [ a11, a12, a13, a21, a22, a23, a31, a32, a33 ])

            retval = dict()
            retval['isometry'] = s
            retval['denominator'] = _Z_to_int(data.denominator)
            retval['src'] = Integer(data.src)
            retval['dst'] = Integer(data.dst)
            yield retval

    def hecke_matrix(self, p, conductor, sparse=None, precise=True):
        prime = Integer(p)

        if not prime.is_prime():
            raise Exception("p is not prime.")

        if self.level_ % conductor != 0:
            raise Exception("Conductor must divide the level.")

        if self.level_ % p == 0:
            raise Exception("Cannot compute Hecke matrix at primes dividing the level.")

        if prime in self.hecke:
            if conductor in self.hecke[prime]:
                return self.hecke[prime][conductor]
            else:
                raise Exception("No Hecke matrix associated to this conductor. How did this happen?")

        if sparse is None:
            density = self._estimated_density(p, conductor)
            logging.info("Hecke matrix density estimate at p=%s: %s%%", prime, density * 100.0)
            sparse = density <= 0.3
            logging.info("Based on density estimates, computing %s matrices", "sparse" if sparse else "dense")

        if not precise:
            if not self.Z64_genus_is_set:
                logging.info("Converting arbitrary precision Genus object to fixed-precision Genus object...")
                self.Z64_genus = make_shared[Genus[Z64]](deref(self.Z_genus))
                self.Z64_genus_is_set = True

            start_time = datetime.now()
            if sparse:
                logging.info("Computing Hecke matrices (p=%s, sparse, int64_t)...", prime)
                self._hecke_matrix_sparse_imprecise(prime)
            else:
                logging.info("Computing Hecke matrices (p=%s, dense, int64_t)...", prime)
                self._hecke_matrix_dense_imprecise(prime)
        else:
            start_time = datetime.now()
            if sparse:
                logging.info("Computing Hecke matrices (p=%s, sparse, arbitrary)...", prime)
                self._hecke_matrix_sparse_precise(prime)
            else:
                logging.info("Computing Hecke matrices (p=%s, dense, arbitrary)...", prime)
                self._hecke_matrix_dense_precise(prime)
        end_time = datetime.now()
        logging.info("Finished computing Hecke matrices at p=%s (time: %s)", prime, end_time-start_time)

        if conductor in self.hecke[prime]:
            return self.hecke[prime][conductor]
        else:
            raise Exception("No Hecke matrix associated to this conductor. How did this happen?")

    def sage_hecke_matrix(self, p, conductor, precise=True, sparse=None):
        prime = Integer(p)

        if prime in self.sage_hecke:
            if conductor in self.sage_hecke[prime]:
                return self.sage_hecke[prime][conductor]

        # Get the desired Hecke matrix and convert it into a sage matrix object.
        A = self.hecke_matrix(p, conductor, precise=precise, sparse=sparse)

        logging.info("Converting numpy matrix to sage matrix (p=%s, conductor=%s)...", p, conductor)

        start_time = datetime.now()
        # Check whether A has the attributes of a sparse matrix.
        if hasattr(A, 'indices') and hasattr(A, 'indptr'):
            # Build a sparse matrix over the Integers(). This is considerably
            # faster than calling matrix(Integers(), A.toarray(), sparse=True)
            # as of 20190214.    ¯\_(ツ)_/¯
            dim = self.dims[conductor]
            B = matrix(Integers(), nrows=dim, ncols=dim, sparse=True)
            pos = 0
            for row,ptr in enumerate(A.indptr[1:]):
                while pos < ptr:
                    col = A.indices[pos]
                    B.add_to_entry(row, col, A.data[pos])
                    pos += 1
        else:
            B = matrix(Integers(), A, sparse=False)
        end_time = datetime.now()

        logging.info("  conversion time: %s", end_time-start_time)

        if not prime in self.sage_hecke:
            self.sage_hecke[prime] = dict()

        self.sage_hecke[prime][conductor] = B
        return self.sage_hecke[prime][conductor]

    def _hecke_matrix_dense_precise(self, Integer p):
        cdef cppmap[Z,vector[int]] mymap
        cdef cppmap[Z,vector[int]].iterator it

        try:
            start_time = datetime.now()
            mymap = deref(self.Z_genus).hecke_matrix_dense(Z(p.value))
            end_time = datetime.now()
            logging.info("  call time: %s", end_time-start_time)
            self.hecke[p] = dict()
        except Exception as e:
            raise Exception(e.message)

        start_time = datetime.now()
        it = mymap.begin()
        while it != mymap.end():
            cond = _Z_to_int(deref(it).first)
            self.hecke[p][cond] = _make_matrix(self.dims[cond], deref(it).second)
            incr(it)

        end_time = datetime.now()
        logging.info("  copy time: %s", end_time-start_time)

    def _hecke_matrix_dense_imprecise(self, Integer p):
        cdef cppmap[Z64,vector[int]] mymap
        cdef cppmap[Z64,vector[int]].iterator it

        try:
            start_time = datetime.now()
            mymap = deref(self.Z64_genus).hecke_matrix_dense(p)
            end_time = datetime.now()
            logging.info("  call time: %s", end_time-start_time)
            self.hecke[p] = dict()
        except Exception as e:
            raise Exception(e.message)

        start_time = datetime.now()
        it = mymap.begin()
        while it != mymap.end():
            cond = Integer(deref(it).first)
            self.hecke[p][cond] = _make_matrix(self.dims[cond], deref(it).second)
            incr(it)

        end_time = datetime.now()
        logging.info("  copy time: %s", end_time-start_time)

    def _hecke_matrix_sparse_precise(self, Integer p):
        cdef cppmap[Z,vector[vector[int]]] mymap
        cdef cppmap[Z,vector[vector[int]]].iterator it

        try:
            start_time = datetime.now()
            mymap = deref(self.Z_genus).hecke_matrix_sparse(Z(p.value))
            end_time = datetime.now()
            logging.info("  call time: %s", end_time-start_time)
            self.hecke[p] = dict()
        except Exception as e:
            raise Exception(e.message)

        cdef int data_len
        cdef int indices_len
        cdef int intptr_len

        start_time = datetime.now()
        it = mymap.begin()
        while it != mymap.end():
            cond = _Z_to_int(deref(it).first)
            dim = self.dims[cond]

            mat = csr_matrix((np.array([]), np.array([]), np.zeros(dim+1)), shape=(dim,dim))

            data_len = deref(it).second[0].size()
            indices_len = deref(it).second[1].size()
            indptr_len = deref(it).second[2].size()

            mat.data = _make_array(data_len, deref(it).second[0])
            mat.indices = _make_array(indices_len, deref(it).second[1])
            mat.indptr = _make_array(indptr_len, deref(it).second[2])

            self.hecke[p][cond] = mat

            incr(it)

        end_time = datetime.now()
        logging.info("  copy time: %s", end_time-start_time)

    def _hecke_matrix_sparse_imprecise(self, Integer p):
        cdef cppmap[Z64,vector[vector[int]]] mymap
        cdef cppmap[Z64,vector[vector[int]]].iterator it

        try:
            start_time = datetime.now()
            mymap = deref(self.Z64_genus).hecke_matrix_sparse(p)
            end_time = datetime.now()
            logging.info("  call time: %s", end_time-start_time)
            self.hecke[p] = dict()
        except Exception as e:
            raise Exception(e.message)

        cdef int data_len
        cdef int indices_len
        cdef int intptr_len

        start_time = datetime.now()
        it = mymap.begin()
        while it != mymap.end():
            cond = Integer(deref(it).first)
            dim = self.dims[cond]

            mat = csr_matrix((np.array([]), np.array([]), np.zeros(dim+1)), shape=(dim,dim))

            data_len = deref(it).second[0].size()
            indices_len = deref(it).second[1].size()
            indptr_len = deref(it).second[2].size()

            mat.data = _make_array(data_len, deref(it).second[0])
            mat.indices = _make_array(indices_len, deref(it).second[1])
            mat.indptr = _make_array(indptr_len, deref(it).second[2])

            self.hecke[p][cond] = mat

            incr(it)

        end_time = datetime.now()
        logging.info("  copy time: %s", end_time-start_time)

    def _estimated_density(self, p, conductor):
        if p == 2: p=p+1
        M = (p+1) / 2
        N = self.dims[conductor]

        if N == 0:
            return 0.0

        if N == 1:
            return 1.0

        if conductor == 1:
            return 1.0 - exp(2 * M * (log(N-1) - log(N)))

        lg2M1 = lgamma(2 * M + 1)
        log4 = log(4)

        logN1 = log(N-1)
        logN = log(N)

        density = 0.0
        for i in range(M):
            logprob = 2*(M-i)*logN1 + lg2M1 - (2*M*logN +
                2*lgamma(i+1) + i*log4 + lgamma(2*M-2*i+1))
            e = exp(logprob)
            if e == 0: break
            density += e

        return 1.0 - density

    def __repr__(self):
        return '''Birch genus with level {} = {}
Ramified primes = {}
Dimensions = {}
Seed = {}'''.format(self.level_, self.facs, self.ramified_primes_, self.dims, self.seed_)

    def __reduce__(self):
        return (BirchGenus,
            (self.level_, self.ramified_primes_, self.seed_),
            self.__getstate__())

    def __getstate__(self):
        state = dict()
        if self.eigenvectors is not None:
            state['eigenvectors'] = list(self.eigenvectors)
        else:
            state['eigenvectors'] = None
        return state

    def __setstate__(self, state):
        if state['eigenvectors'] is not None:
            self.eigenvectors = list(state['eigenvectors'])
        else:
            self.eigenvectors = None

cdef _Z_to_int(const Z& x):
    return Integer(x.get_str(10), 10)

cdef class _MatrixWrapper:
    cdef vector[int] vec
    cdef Py_ssize_t shape[2]
    cdef Py_ssize_t strides[2]
    cdef dim

    def __init__(self, dim=0):
        self.dim = dim

    cdef set_data(self, vector[int]& data):
        self.vec = move(data)

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef Py_ssize_t itemsize = sizeof(self.vec[0])

        self.shape[0] = self.dim
        self.shape[1] = self.dim
        self.strides[0] = sizeof(int) * self.dim
        self.strides[1] = sizeof(int)
        buffer.buf = <char *>&(self.vec[0])
        buffer.format = 'i'
        buffer.internal = NULL
        buffer.itemsize = itemsize
        buffer.len = self.dim * self.dim * itemsize
        buffer.ndim = 2
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self.shape
        buffer.strides = self.strides
        buffer.suboffsets = NULL

cdef class _ArrayWrapper:
    cdef vector[int] vec
    cdef Py_ssize_t shape[1]
    cdef Py_ssize_t strides[1]
    cdef dim

    def __init__(self, dim=0):
        self.dim = dim

    cdef set_data(self, vector[int]& data):
        self.vec = move(data)

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef Py_ssize_t itemsize = sizeof(self.vec[0])

        self.shape[0] = self.dim
        self.strides[0] = sizeof(int)
        buffer.buf = <char *>&(self.vec[0])
        buffer.format = 'i'
        buffer.internal = NULL
        buffer.itemsize = itemsize
        buffer.len = self.dim * itemsize
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self.shape
        buffer.strides = self.strides
        buffer.suboffsets = NULL

cdef _make_array(dim, vector[int]& data):
    cdef _ArrayWrapper aw
    aw = _ArrayWrapper(dim)
    aw.set_data(data)
    return np.asarray(aw)

cdef _make_matrix(dim, vector[int]& data):
    cdef _MatrixWrapper mw
    mw = _MatrixWrapper(dim)
    mw.set_data(data)
    return np.asarray(mw)

cdef do_something(const IsometrySequenceData[Z]& data):
    print(data.src, data.dst)
