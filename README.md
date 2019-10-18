Transire generates and manipulates crystalline interface structures built 
on the Atomic Simulation Environment (ASE).   A number of search methods 
are available to minimize the interface energy, including a random walk 
Markov Chain search and rotations of one slab relative to the other.  
Energy calculations are performed using CP2K or LAMMPS, and any method 
available in either can be used for calculating the interface energy.
The electron transport properties across the interface can be calculated
using the non-equilibrium Green's function (NEGF) method implemented in ASE
(requires CP2K calculation).

Transire is Latin for "to cross over" in reference to transport
across the interface.

Requirements
==============

Python 2.7

ASE 3.12 or newer

pycp2k (if using cp2k)
        
        available from project page: https://github.com/SINGROUP/pycp2k

cp2k 4.1 or newer (optional)
        
        If not using cp2k, make sure all energy related flags are set to
        False in the input file 

LAMMPS (26 Jan 2017) or newer (optional)
        
        install serial version with mode=shlib and install-python options

lammpslib
        
        https://svn.fysik.dtu.dk/projects/ase-extra/trunk/ase/calculators
        
        copy the file into your local library folder, usually
        ~/.local/lib/python2.7/site-packages which is where python installs
        the other libraries

Installation
==============

Installation performed using setup.py in the main directory:

```
>>>  python setup.py install
```

To install into the local user space:

```
>>>  python setup.py install --user
```

If transire is installed in local space, ensure that the path environmental variables are set.
The $PYTHONPATH variable needs to include ~/.local/lib/pythonX/site-packages where pythonX is
the version of python on your system.  The $PATH variable needs to include ~/.local/bin.
Example commands to set environment variables (cshell):

```
>>>  setenv PYTHONPATH $(PYTHONPATH):~/.local/lib/pythonX/site-packages
>>>  setenv PATH $(PATH):~/.local/bin
```

To test for successful installation use:

```
>>>  transire.py -h
```  

or

```
>>>  python -m transire.py -h
```

If installed correctly, the help description of input options for transire will be printed.


Documentation
=============

Documentation is designed for sphinx v1.7.0 and can be compiled by running:

>>> make html

in the 'doc' directory.  Documentation output defaults to 'doc/html/'.


Citations
=========

If you find ARL Transire useful, please cite its [![DOI](https://zenodo.org/badge/164665839.svg)](https://zenodo.org/badge/latestdoi/164665839)

Legal
=====

Transire is available under the GNU Lesser General Public License (LGPL).
A copy of the LGPL is provided in the LICENSE file.

Copyright (2018) Caleb Carlin.
