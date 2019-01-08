.. module:: source.energy_calculators

Energy Calculators
============================

A class for streamlining the interaction with external calculators.
Currently CP2K and LAMMPS are included through the python wrappers
pycp2k and lammpslib respectively instead of the ASE internal interface.
This is done to allow the user greater control.

The energy calculation is called by the class and all options are
made in the input file.  Each calculator sets up the necessary input
file for the external program, executes the program, and returns the
total energy of the interface structure.  The total energy is divided
by the area of the X and Y vectors of the interface structure cell.
This energy per area is returned as the output of the class.

If binding energy calculation is turned on, then an energy calculation
is carried out on each side of the interface as slabs bounded by vacuum.
The returned energy per area is then the difference between the total
energy and the sum of the slab energies divided by the X-Y area of the
interface structure.

For samples of the inputs required for each external program, look
in the examples directory.

____________

.. autoclass:: EnergyCalc
