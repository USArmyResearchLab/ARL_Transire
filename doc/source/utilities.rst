.. module:: source.utilities

Utility Routines
================

This is a collection of routines and classes that are either general
purpose or used throughout the rest of the code.

____________

.. autofunction:: lcm

====

.. autofunction:: write_xyz

====

.. autoclass:: InterfaceConfigStorage
    :members:

====

.. autofunction:: generate_surfaces_list

====

.. autofunction:: printx

====

.. autofunction:: angle_between

====

.. autofunction:: unit_vector

====

In version 3.13 of ASE changed the behavior of functions in ase.build
so that the length of the third lattice  vector is not conserved.
As this breaks functionality in Transire, the two following functions
are taken from ASE 3.12.

.. autofunction:: surface_from_ase

.. autofunction:: build_from_ase

====

.. autofunction:: ext_gcd

====

.. autofunction:: almost_zero
