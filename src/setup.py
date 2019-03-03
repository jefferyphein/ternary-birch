from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "ternary-birch",
    version = "2.0.0-beta.3",
    ext_modules = cythonize("ternary_birch.pyx")
)
