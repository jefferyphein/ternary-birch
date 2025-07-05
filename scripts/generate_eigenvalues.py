import argparse
import json
import logging
import multiprocessing
import operator
import os
import pathlib
import random
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from functools import partial, reduce

from sage.rings.integer import Integer
from sage.modular.dims import sturm_bound
from ternary_birch import BirchGenus

LOGGER = logging.getLogger(__name__)


@dataclass
class Job:
    level: int
    seed: int
    ramified_primes: list[int]
    ramified_primes_product: int
    unramified_primes: list[int]
    unramified_primes_product: int
    sturm_bound_of_ramified_primes_product: int
    eigenvalues_upto: int
    elapsed_seconds: float | None = None
    elapsed_seconds_genus: float | None = None
    elapsed_seconds_eigenvectors: float | None = None
    elapsed_seconds_eigenvalues: float | None = None


@dataclass
class JobResult:
    job: Job
    dimensions_by_conductor: dict[int, int]
    dimension_total: int
    rational_eigenvectors_count: int
    rational_eigenvectors: list[dict]

    def json(self, **kwargs) -> str:
        return json.dumps(asdict(self), **kwargs)


@dataclass
class Stopwatch:
    elapsed_sum: int = 0
    start_mark: int | None = None

    def start(self):
        if self.start_mark is None:
            self.start_mark = time.monotonic_ns()

    def stop(self):
        if self.start_mark is not None:
            self.elapsed_sum += time.monotonic_ns() - self.start_mark
            self.start_mark = None

    def elapsed(self) -> float:
        """Convert nanoseconds into seconds."""
        total_elapsed = self.elapsed_sum + (
            0 if self.start_mark is None else time.monotonic_ns() - self.start_mark
        )
        return total_elapsed / 10**9

    @staticmethod
    @contextmanager
    def contextmanager():
        stopwatch = Stopwatch()
        stopwatch.start()
        try:
            yield stopwatch
        finally:
            stopwatch.stop()


def is_squarefree(n: int, *, num_primes: int | None = None):
    """Whether provided integer is squarefree with correct number of divisors."""
    primes = Integer(n).prime_divisors()
    # Is n squarefree?
    if reduce(operator.mul, primes) != n:
        return False
    if num_primes is not None:
        return len(primes) == num_primes
    return True


def generate_jobs(level: int, *, seed: int, upto: int | None = None) -> Job:
    """Generate jobs for calculating rational eigenvectors at the provided level.

    Preconditions:
        - level is squarefree

    """
    primes = Integer(level).prime_divisors()
    ramified_primes = sorted(primes)[:-1] if len(primes) % 2 == 0 else sorted(primes)
    ramified_primes_product = reduce(operator.mul, ramified_primes)
    unramified_primes = [sorted(primes)[-1]] if len(primes) % 2 == 0 else []
    unramified_primes_product = reduce(operator.mul, unramified_primes, 1)
    the_sturm_bound = sturm_bound(ramified_primes_product, 2)
    upto = min(upto or the_sturm_bound, the_sturm_bound)
    return Job(
        level=level,
        seed=seed,
        ramified_primes=list(map(int, ramified_primes)),
        ramified_primes_product=int(ramified_primes_product),
        unramified_primes=list(map(int, unramified_primes)),
        unramified_primes_product=int(unramified_primes_product),
        eigenvalues_upto=int(upto),
        sturm_bound_of_ramified_primes_product=int(the_sturm_bound),
    )


def normalize_vectors(vecs: dict):
    """Convert all Sage objects into native Python objects so they can be serialized."""
    for vec in vecs:
        vec["vector"] = tuple(int(coef) for coef in vec["vector"])
        vec["conductor"] = int(vec["conductor"])
        vec["aps"] = {int(p): int(e) for p, e in vec["aps"].items()}
    return vecs


def find_rational_eigenvectors(job: Job) -> JobResult:
    """Find all rational eigenvectors and calculate eigenvalues up to Sturm bound.

    Preconditions:
        - job.level is squarefree

    Calculates as many eigenvalues as needed according to the Sturm bound of the fully
    ramified level. For example, if level = p*q with p<q, calculate all eigenvalues up
    to the Sturm bound required by level p.

    """
    with Stopwatch.contextmanager() as stopwatch:
        with Stopwatch.contextmanager() as stopwatch_genus:
            genus = BirchGenus(job.level, seed=job.seed)
        with Stopwatch.contextmanager() as stopwatch_eigenvectors:
            vecs = genus.rational_eigenvectors()
        with Stopwatch.contextmanager() as stopwatch_eigenvalues:
            genus.compute_eigenvalues_upto(job.eigenvalues_upto)
    job.elapsed_seconds = stopwatch.elapsed()
    job.elapsed_seconds_genus = stopwatch_genus.elapsed()
    job.elapsed_seconds_eigenvectors = stopwatch_eigenvectors.elapsed()
    job.elapsed_seconds_eigenvalues = stopwatch_eigenvalues.elapsed()
    return JobResult(
        job=job,
        dimensions_by_conductor={
            int(conductor): int(dimension)
            for conductor, dimension in genus.dimensions().items()
        },
        dimension_total=int(genus.dimension()),
        rational_eigenvectors_count=len(vecs),
        rational_eigenvectors=normalize_vectors(vecs),
    )


if __name__ == "__main__":
    # Set up logging.
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    # Parse command-line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("savepath", type=pathlib.Path)
    parser.add_argument(
        "--lower",
        type=int,
        required=False,
        default=500000,
        help="Lower bound of integers to check",
    )
    parser.add_argument(
        "--upper",
        type=int,
        required=False,
        default=1000000,
        help="Upper bound of integers to check",
    )
    parser.add_argument("--num-primes", type=int, required=False, default=None)
    parser.add_argument(
        "--upto",
        type=int,
        default=None,
        help=(
            "Compute eigenvalues up to this bound. "
            "Defaults to the Sturm bound for the level."
        ),
    )
    parser.add_argument(
        "--seed", type=int, required=False, default=random.randint(0, 2**32 - 1)
    )
    parser.add_argument(
        "--max-workers", type=int, required=False, default=multiprocessing.cpu_count()
    )
    args = parser.parse_args()

    # Use as many workers as requested.
    pool = ProcessPoolExecutor(max_workers=args.max_workers)

    # Generate list of squarefree primes.
    squarefree_primes = filter(
        partial(is_squarefree, num_primes=args.num_primes),
        range(args.lower, args.upper),
    )
    jobs = map(
        partial(generate_jobs, seed=args.seed, upto=args.upto), squarefree_primes
    )
    tasks = [pool.submit(find_rational_eigenvectors, job=job) for job in jobs]
    LOGGER.info("Submitting %s tasks using %s worker(s)", len(tasks), args.max_workers)
    for n, task in enumerate(as_completed(tasks)):
        try:
            result = task.result()
        except Exception:
            LOGGER.exception("Unexpected error occurred while processing task.")
        else:
            savepath = args.savepath / f"{result.job.level}.json"
            with open(savepath, "w") as f:
                f.write(result.json())
            LOGGER.info(
                "(%s / %s) Finished computing eigenvalues at level %s",
                n + 1,
                len(tasks),
                result.job.level,
            )
