from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
extensions = [Extension('scqo', ['scqo.pyx', 'ode.c', 'pool.c', 'operators.c', 'scquantumoptics.c'])]
setup(
    ext_modules = cythonize(extensions)
)
