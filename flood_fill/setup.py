from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules = cythonize("flood_fill.pyx"),
    include_dirs=[numpy.get_include()]
)

# there are several ways of building a cython code on http://docs.cython.org/src/quickstart/build.html and here is the distutils one. include_dirs here for give the correct of numpy loading problem.