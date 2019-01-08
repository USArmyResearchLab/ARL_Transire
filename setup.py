#!/usr/bin/env python

import sys
from distutils.core import setup
from setuptools import find_packages

long_description = """\
Transire is a python package built on the ASE package 
in combination with methods including CP2K."""

if sys.version_info < (2, 6, 0, 'final', 0):
    raise SystemExit('Python 2.6 or later is required!')

setup(name='transire',
      version='1.0',
      description="Interface and electron transport tools",
      author="Caleb Carlin",
      author_email='caleb.m.carlin.ctr@mail.mil',
      package_dir={'transire': ''},
      packages=find_packages(),
      long_description = long_description,
      scripts = ['transire.py'])
