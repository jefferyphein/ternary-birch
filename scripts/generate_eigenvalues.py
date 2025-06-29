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
    sturm_bound_of_ramified_primes_product: int


@dataclass
class SubJob:
    level: int
    seed: int
    sturm_bound: int


@dataclass
class JobResult:
    job: Job
    elapsed_seconds: float
    rational_eigenvectors_count: int
    rational_eigenvectors: list[dict]

    def json(self, **kwargs) -> str:
        return json.dumps(asdict(self), **kwargs)


@dataclass
class SubJobResult:
    subjob: SubJob
    elapsed_seconds: float
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


def is_squarefree_six_prime(n: int) -> bool:
    """Determines whether provided integer is squarefree with six prime divisors."""
    primes = Integer(n).prime_divisors()
    # Does n has six prime divisors?
    if len(primes) != 6:
        return False
    # Is n squarefree?
    return reduce(operator.mul, primes) == n


def generate_jobs(level: int, seed: int) -> Job:
    """Generate jobs for calculating rational eigenvectors at the provided level.

    Preconditions:
        - level is squarefree
        - level has an even number of prime divisors

    """
    primes = Integer(level).prime_divisors()
    ramified_primes = sorted(primes)[:-1]
    ramified_primes_product = reduce(operator.mul, ramified_primes)
    return Job(
        level=level,
        seed=seed,
        ramified_primes=list(map(int, ramified_primes)),
        ramified_primes_product=int(ramified_primes_product),
        sturm_bound_of_ramified_primes_product=int(
            sturm_bound(ramified_primes_product, 2)
        ),
    )


def generate_subjobs(level: int, seed: int) -> SubJob:
    """Generate subjobs for calculating rational eigenvectors at a fully ramified level.

    Preconditions:
        - level is squarefree
        - level has an odd number of prime divisors

    """
    return SubJob(
        level=level,
        seed=seed,
        sturm_bound=int(sturm_bound(level, 2)),
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
        - job.level has an even number of prime divisors

    Calculates as many eigenvalues as needed according to the Sturm bound of the fully
    ramified level. For example, if level = p*q with p<q, calculate all eigenvalues up
    to the Sturm bound required by level p.

    """
    with Stopwatch.contextmanager() as stopwatch:
        genus = BirchGenus(job.level, seed=job.seed)
        vecs = genus.rational_eigenvectors()
        genus.compute_eigenvalues_upto(job.sturm_bound_of_ramified_primes_product)
    return JobResult(
        job=job,
        elapsed_seconds=stopwatch.elapsed(),
        rational_eigenvectors_count=len(vecs),
        rational_eigenvectors=normalize_vectors(vecs),
    )


def find_rational_eigenvectors_at_fully_ramified_level(subjob: SubJob) -> SubJobResult:
    """Find all rational eigenvectors and calculate eigenvalues up to Sturm bound.

    Preconditions:
        - subjob.level is squarefree
        - subjob.level has an odd number of prime divisors

    Calculates as many eigenvalues as needed according to the Sturm bound of the level.

    """
    with Stopwatch.contextmanager() as stopwatch:
        genus = BirchGenus(subjob.level, seed=subjob.seed)
        vecs = genus.rational_eigenvectors()
        the_sturm_bound = sturm_bound(subjob.level, 2)
        genus.compute_eigenvalues_upto(the_sturm_bound)
    return SubJobResult(
        subjob=subjob,
        elapsed_seconds=stopwatch.elapsed(),
        rational_eigenvectors_count=len(vecs),
        rational_eigenvectors=normalize_vectors(vecs),
    )


if __name__ == "__main__":
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
    parser.add_argument(
        "--seed", type=int, required=False, default=random.randint(0, 2**32 - 1)
    )
    parser.add_argument(
        "--max-workers", type=int, required=False, default=multiprocessing.cpu_count()
    )
    args = parser.parse_args()

    # Attempt to create the path for saving oldforms.
    oldforms_path = args.savepath / "oldforms"
    os.makedirs(oldforms_path, exist_ok=True)

    # Use as many workers as requested.
    pool = ProcessPoolExecutor(max_workers=args.max_workers)

    # Generate list of jobs for calculating eigenforms for squarefree 6-primes.
    six_primes = filter(is_squarefree_six_prime, range(args.lower, args.upper))
    jobs = list(map(partial(generate_jobs, seed=args.seed), six_primes))

    # Deduplicate all observed fully ramified levels.
    fully_ramified_levels = set(job.ramified_primes_product for job in jobs)

    # Generate old forms.
    subjobs = map(partial(generate_subjobs, seed=args.seed), fully_ramified_levels)
    subtasks = [
        pool.submit(find_rational_eigenvectors_at_fully_ramified_level, subjob)
        for subjob in subjobs
    ]

    # Process each subtask as it completes and save results to disk.
    for subtask in as_completed(subtasks):
        try:
            result = subtask.result()
        except Exception:
            LOGGER.exception("Unexpected error occurred while processing task.")
        else:
            savepath = oldforms_path / f"{result.subjob.level}.json"
            with open(savepath, "w") as f:
                f.write(result.json(indent=2))

    # Now that the subtasks are completed, process the primary jobs.
    tasks = [pool.submit(find_rational_eigenvectors, job=job) for job in jobs]
    for task in as_completed(tasks):
        try:
            result = task.result()
        except Exception:
            LOGGER.exception("Unexpected error occurred while processing task.")
        else:
            savepath = args.savepath / f"{result.job.level}.json"
            with open(savepath, "w") as f:
                f.write(result.json(indent=2))
