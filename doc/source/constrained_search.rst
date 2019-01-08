.. module:: source.constrained_search

The Markov Chain Random Walk
============================

Class for generating interface structures by twist rotations and/or
translations with random step sizes.  Each step is used as the starting
configuration for the following step and the lowest energy interface
structure is returned at the end.  In addition, the translations and
rotation necessary to reproduce any step in the Markov Chain are printed
to the output file.

Twist rotations occur as a rotation around the axis perpendicular to
the interface (Z-axis).  Each time a rotation occurs, the interface
structure is regenerated from the unit cells to ensure orthornomarlity
of the cell vectors.

Translations occur by shifting the atoms in the top (crystal a) and then
wrapping the atoms back into the interface using the crystaline periodic
boundary condition.

____________

.. autoclass:: ConstrainedSearch
