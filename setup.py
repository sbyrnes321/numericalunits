# -*- coding: utf-8 -*-

import os
from setuptools import setup
from numericalunits import __version__

# Utility function to read the README file.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

descrip = ("A package that lets you define quantities with units, which can "
           "then be used in almost any numerical calculation in any "
           "programming language. Checks that calculations pass dimensional "
           "analysis, performs unit conversions, and defines physical "
           "constants.")

setup(
    name = "numericalunits",
    version = __version__,
    author = "Steven Byrnes",
    author_email = "steven.byrnes@gmail.com",
    description = descrip,
    license = "MIT",
    keywords = "units, quantities, physical constants, dimensional analysis",
    url = "http://pypi.python.org/pypi/numericalunits",
    py_modules=['numericalunits'],
    long_description=read('README.rst'),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3"]
)
