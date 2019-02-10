from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "ternary-birch",
    version = "2.0.0-beta.2",
    ext_modules = cythonize("py_birch.pyx")
)
