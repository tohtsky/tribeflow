from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os
import sys
import setuptools

__version__ = '0.0.1'

def get_eigen3_include():
    EIGEN3_INCLUDE_DIR = os.environ.get('EIGEN3_INCLUDE_DIR', None)
    if EIGEN3_INCLUDE_DIR is None:
        raise RuntimeError('''Specify Eigen3 include directory by the Env variable
    ``EIGEN3_INCLUDE_DIR``
''')
    return EIGEN3_INCLUDE_DIR

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


ext_modules = [
    Extension(
        'tribeflowpp',
        ['tribeflow/plearn.cpp', 'tribeflow/dataio.cpp',
            'tribeflow/learn_body.cpp', 'tribeflow/kernels/base.cpp'],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True),
            get_eigen3_include()
        ],
        language='c++'
    ),
]


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        opts.append('-O3') 
        opts.append('-std=c++17')
        opts.append('-march=native')
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)

setup(
    name='tribeflowpp',
    version=__version__,
    description='A tribeflow re-implementation using C++11',
    long_description='',
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.2'],
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
)

