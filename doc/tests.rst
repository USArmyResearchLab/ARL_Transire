.. _tests:

Tests
=====

Test inputs for ARL_Transire are provided for using CP2K, LAMMPS, and running
ARL_Transire without energy calculations.  The tests are located in the tests
folder and are further divided into three folders.

1) CP2K_INPUTS

2) LAMMPS_INPUTS

3) NON_ENERGY_INPUTS

Each folder contains a set of test inputs to cover the major functionality of
ARL_Transire features applicable to the energy method being used as well as the
required additional input files.  The test inputs can be run individually or
as a batch using the respective scripts in the tests folder.

CP2K
____

It is recommended that the CP2K tests be run in parallel as running the entire
set of tests can take a while, despite limiting the number of atoms used in
the tests.  ARL_Transire output files for the tests follow the pattern:

cp2k_test_X.out where X matches the corresponding input cp2k_X.in

**Tests include**:

*cp2k_angle.in*
        calculates the energy of bulk 3c SiC and the symmetry identical
        90 degree rotation.

*cp2k_mc.in*
        calculates the energy of bulk 3c SiC and then performs 20 steps
        of the Markov Chain random walk.  The small max number of allowed
        atoms and the low energy starting configuration means that the
        output should show only steps rejected for higher energy and for
        exceeding the size constraint.

*cp2k_et.in*
        calculates the energy of the 001/101 3c SiC grain boundary and
        reads in the results to perform the electron transmission
        calculation.

*cp2k_insert.in*
        inserts two Pt atoms inbetween bulk 3c SiC and calculates the
        energy.

*cp2k_separation.in*
        performs a separtion optimization on the bulk 3c SiC interface.

LAMMPS
______

LAMMPS output files follow the same naming convention as for cp2k
with lammps\_ substituted for cp2k\_.

*lammps_angle.in*
        generates interface structures of bulk Cu for 9 twist rotations with a
        step size of 10 degrees.  A RAS search is further performed and the
        energy of the result structures are calculated.

*lammps_mc.in*
        calculates the energy of bulk Cu and then performs 100 steps of
        the Markov Chain random walk.  The low energy of the starting
        interface structure means that the output should report that most
        steps were rejected for being higher in energy.

*lammps_surface.in*
        calculates the energy of bulk Cu and then repeats the calculation
        for grain boundaries involving the 001/001, 001/100, and 001/101
        sets of Miller Indices.  Symmetric sets of Miller Planes must
        report the same energy.

*lammps_insert.in*
        calculates the energy of bulk Cu and then inserts two Cu atoms at the
        interface and calculates the inserted energy.

*lammps_sep.in*
        performs a separation optimization on bulk Cu.


Non-energy
__________

ARL_Transire offers a few features that do not require an energy calculation.
These options are useful for testing parameters before performing the energy
calculation or exporting the resulting structures to a non-supported program,

*non_energy_angle.in*
        generates structures for a stacked pair of graphene sheets under
        twist rotation in the range of 0 to 60 degrees with 10 degree
        step size.

*non_energy_insert.in*
        generates a stacked pair of graphene sheets with a rotated pait of
        Pt atoms inserted at the interface.
