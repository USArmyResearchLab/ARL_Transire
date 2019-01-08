# copyright Caleb Michael Carlin (2018)
# Released uder Lesser Gnu Public License (LGPL)
# See LICENSE file for details.

from .utilities import InterfaceConfigStorage, printx
import re
from ase import Atoms
import sys


class EnergyCalc(object):
    """
    Class for handling all energy calculations.
    
    Parameters:

    input: InputReader object
        object with read in keywords
    method: string
        lowercase string that determines which program to use
        when calculating the energy
    ics: InterfaceConfigStorage object
        object containing the interface structure of which
        the energy will be determined
    exetype: string ('run' or 'write')
        option to write CP2K input file without executing the
        program
    """

    def __init__(self, input, method, ics, exetype):
        self.method = method
        self.input = input
        self.config = ics
        self.exetype = exetype
        if (method == 'cp2k'):
            self.config.energy = self.energy_cp2k()
        if (method == 'lammps'):
            self.config.energy = self.energy_lammps()
        ics.energy = self.config.energy

    def energy_cp2k(self):
        """Quick routine to calculate total energy using CP2K
        """
        from pycp2k import CP2K
        calc = CP2K()
        calc.working_directory = self.input.dict['working_directory']
        # Determine if using less than max_processers would be better
        if (int(self.input.dict['atoms_per_process']) != 0):
            natoms = len(self.config.atom)
            processes = natoms / float(self.input.dict['atoms_per_process'])
            if (processes <= 1.0):
                processes = 1
                printx("atoms_per_process may be set too high for system size")
                printx("Number of Atoms = " + str(natoms))
            elif (processes > int(self.input.dict['max_mpi_processes'])):
                processes = int(self.input.dict['max_mpi_processes'])
            else:
                processes = int(processes)
        else:
            processes = int(self.input.dict['max_mpi_processes'])
        calc.mpi_n_processes = processes
        self.config.atom.set_calculator(calc)
        # load in cp2k options
        execfile(self.input.dict['cp2k_input'])
        calc.project_name = self.config.file_name
        calc.write_input_file()
        # If we only want the input file written, then we stop here
        if (self.exetype == 'write'):
            energy_per_area = 2.0e10
            return energy_per_area
        # if print_debug is true, we run energy calc without "try" to
        # allow for better error message passing
        if (self.input.dict['print_debug'] != 'False'):
            calc.run()
        else:
            try:
                calc.run()
            except Exception as err:
                printx('error occured in the energy calculation.')
                energy_per_area = 1.0e10
                return energy_per_area

        # Read in the output file to find the final energy.
        with open(calc.output_path, "r") as fin:
            regex = re.compile(
                " ENERGY\| Total FORCE_EVAL \( QS \)"
                " energy \(a\.u\.\):\s+(.+)\n")
            for line in fin:
                match = regex.match(line)
                if match:
                    printx("Final energy: {}".format(match.groups()[0]))
                    energy = match.groups()[0]
        area = self.config.atom.cell[0, 0] * self.config.atom.cell[1, 1]
        if (self.input.dict['print_debug'] != 'False'):
            printx('Area used for energy '+self.config.file_name+'is'
                    +str(area))
        energy_per_area = float(energy) / area

        return energy_per_area

    def energy_lammps(self):
        """
        Routine that sets up and runs a LAMMPS calculation through LAMMPSlib.
        A file that gives the LAMMPS commands as a list of strings is needed.
        """
        try:
            from lammpslib import LAMMPSlib
        except Exception as err:
            printx("Could not import lammpslib")
            sys.exit(1)

        parallel = True
        try:
            import mpi4py
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
        except Exception as err:
            printx("mpi4py not found, using serial")
            parallel = False

        # load the LAMMPS specific commands.  If file not found,
        # then we close things right away so we don't waste effort
        try:
            exec(compile(open(self.input.dict['lammps_input']).read(),
                 self.input.dict['lammps_input'], 'exec'))
        except Exception as err:
            printx("ERROR: couldn't find or load lammps input file")
            sys.exit(1)

        # Add command so that the box doesn't get collapsed by LAMMPS
        # due to the default behavior of lammpslib
        commands.append('change_box all boundary p p m')
        output_file = self.config.file_name + '.log'
        calc = LAMMPSlib(lmpcmds=commands,
                         log_file=output_file, keep_alive=False)

        self.config.atom.set_calculator(calc)
        if (parallel):
            me = MPI.COMM_WORLD.Get_rank()
            nprocs = MPI.COMM_WORLD.Get_size()
            mpi4py.rc.threaded = False
            printx("Proc %d out of %d procs has" % (me, nprocs))

        # perform the energy calculation.  Normally LAMMPS will kill
        # the python process if an error occurs, but there may be
        # some errors that can be caught without complete closure
        printx("Beginning LAMMPS energy calculation")
        try:
            total_energy = self.config.atom.get_potential_energy()
        except Exception as err:
            printx("ERROR in LAMMPS energy calculation")
            energy_per_area = 2.0e10
            return energy_per_area

        printx("Finished LAMMPS energy calculation")
        area = self.config.atom.cell[0, 0] * self.config.atom.cell[1, 1]
        energy_per_area = float(total_energy) / area

        return energy_per_area


def calculate_binding_energy(ics, input):
    """
    calculates the energy of each slab in the vacuum for determining the
    binding energy of the full interface
    """

    side_1 = ics.tag_copy(1)
    printx("side_1 =" + str(side_1.atom.cell))
    side_2 = ics.tag_copy(2)
    printx("side_2 =" + str(side_2.atom.cell))

    EnergyCalc(input, input.dict['energy_method'], side_1, 'run')
    EnergyCalc(input, input.dict['energy_method'], side_2, 'run')

    area = ics.atom.cell[0][0] * ics.atom.cell[1][1]

    printx("ICS ENERGY = "+str(ics.energy))
    printx("side 1 ENERGY = "+str(side_1.energy))
    printx("side 2 ENERGY = "+str(side_2.energy))

    return (ics.energy - side_1.energy - side_2.energy)
