#! /usr/bin/env python3

from distutils.core import setup
from Cython.Build import cythonize

# Usage:
# python3 compile_cython.py build_ext --inplace

if __name__ == '__main__':
    setup(ext_modules = cythonize('motif_finder_tool_modules.pyx'))