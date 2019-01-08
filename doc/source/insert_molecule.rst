.. module:: source.insert_molecule

Molecule or Crystal Insertion
=============================

Class for inserting the contents of a coordinate file in any format
that ASE supports in between the two crystals at the interface.
By default, a certical space is created at the interface is created
equal to the Z dimension of the inserted object and the inserted object
is centered so that the center of the inserted object is colinear with
the centers of the two sides of the interface structure.

Supported options allow for adjusting the spacing above and below the
inserted object as well as rotating the inserted object around some or
all of the cardinal axes before inserting.

.. caution::
   performing a different search after the insertion will remove the
   inserted object.

____________

.. autoclass:: InsertMolecule
