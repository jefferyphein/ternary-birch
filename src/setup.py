from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension("ternary_birch", ["ternary_birch.pyx"],
        libraries=['gmp', 'gmpxx'],
    )
]

setup(
    name = "ternary-birch",
    version = "2.0.0-beta.3",
    author = "Jeffery Hein",
    author_email = "jefferyphein@gmail.com",
    maintainer = "Jeffery Hein",
    maintainer_email = "jefferyphein@gmail.com",
    url = "https://github.com/jefferyphein/ternary-birch",
    description = "A library for computing Hecke matrices, eigenvectors, and eigenvalues for positive definite rational ternary quadratic forms.",
    license = "GNU GPL",
    ext_modules = cythonize(extensions)
)
