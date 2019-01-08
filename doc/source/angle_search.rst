.. module:: source.angle_search

The Twist Angle Search
======================

Class for generating interface structures by rotating one crsytal around the
z-axis relative to the other crystal.  This class can be used in a combination
of three ways:

    1) Calculate the energy of each new configuration and return the
       configuration with the lowest energy per interface area

    2) Calculate the energy of each new configuration and print out the angle
       and energy data to an ootput file

    2) Perform the rotations and output the coordinates of the interface
       structure using the xyz format.

The angles can be specified in one of three ways:

    1) specify each angle in a list using angles_list.

    2) specify the number of angles (number_of_angles) and the step size
       between each angle (angles_stepsize)

    3) specify the output file from a previous angle search that used the
       reducing angle search option to read in the list of angles.

All angles are specified in degrees and all rotations occur relative to the
starting configuration

The 'reducing angle search' or RAS automates the process of finding the angles
in a range that result in interfaces with the fewest number of total atoms.
The initial range of angles and stepsize is specified by the user as well as
the number of smallest interfaces to report (ras_factor) and the number of
levels to include in the search (ras_depth).  Each level after the first
involves decreasing the stepsize by a factor of 10 and testing the 10 angle
steps on either side of each angle located in the previous layer.  The
'ras_energy' option determines whether the energy is calculated for each of the
ngles that produce the smallest interfaces at the end of the search.

The separation between the interfaces can be optimized after each rotation to
allow the interface to relax using a non-gradient based method that locates a
local minimum to within the specified tolerance.

____________

.. autoclass:: AngleSearch
