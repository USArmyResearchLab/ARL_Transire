.. module:: source.interface

Interface Structure Constructor
===============================

Class for generating and directly modifying the interface structures
from unit cells to the full interface structure.
The structure is defined by the two unit cells making up the two sides of
interface, the Miller Indices for each side, and the number of layers
to be replicated away from the interface for each side.

Building the interface structure involves the following:

1) Determine the linear combination of unit cell vectors that results
   in axis intersections and construct an orthonormal representation of
   each crystal.

#) Populate the orthonormal cells with atoms from periodic translations
   of the unit cells.  Like a cookie cutter from the bulk crystal.

#) Determine the smallest number of each orthornomal cells needed to stack
   the two with stress below the user defined tolerance.

#) the two crystals are stacked along the Z axis with an optional vacuum
   separation at the interface.

Changes to the unit_cells for generating interfaces made by other the
parts of the code should be done by setting self.cut_cell_a based on
ICS.unit_cell_a where ICS is the interface storage unit currently being used.

____________

.. autoclass:: InterfaceSupercell
