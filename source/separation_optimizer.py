# copyright Caleb Michael Carlin (2018)
# Released uder Lesser Gnu Public License (LGPL)
# See LICENSE file for details.

from .utilities import InterfaceConfigStorage, printx, almost_zero
import numpy as np
from .energy_calculators import EnergyCalc


class SeparationOpt(object):
    """
    Class to determine the optimimum separation between the two sides
    of the interface.

    Parameters:

    ics: InterfaceConfigStorage object
        object containing the interface to be manipulated
    input: InputReader object
        object with read in keywords
    surf_a: tuple
        tuple of Miller Indices for the top crystal. Used for generating
        filenames.
    surf_b: tuple
        tuple of Miller Indices for the bottom crystal. Used for generating
        filenames.
    """

    def __init__(self, ics, input, surf_a, surf_b):
        self.input = input
        self.structure = ics.copy()
        self.surf_a = surf_a
        self.surf_b = surf_b
        self.direction = None
        self.step_size = float(self.input.dict['sep_initial_step'])
        self.initial_sep = float(self.input.dict['separation'])
        self.file = self.input.dict['output_file']
        dist = self.locate_minima()
        input.dict['sep_guess'] = dist
        ics = self.structure.copy(no_filename=True)
        #ics.enery = self.structure.energy

    def locate_minima(self):
        """
        Instead of generating the structure each time with a different
        separation we shift the atoms with the tag 2 and change the cell size.
        """
        temp_structure = self.structure.copy()
        # if initial guess and separation are the same then no optimization
        # will take place.  Here we force a difference in energy.
        energy_1 = temp_structure.energy
        delta_d = float(self.input.dict['sep_guess']) - self.initial_sep
        dist = float(self.input.dict['sep_guess'])
        self.direction = cmp(delta_d, 0.0)

        for i in range(int(self.input.dict['sep_max_steps'])):
            # translate the atoms on one side
            temp_structure.atom.positions[:, 2] = [
                x.position[2] + delta_d if x.tag == 2 else x.position[2]
                for x in temp_structure.atom]
            temp_structure.atom.cell[2, 2] += delta_d
            temp_structure.file_name = "sep_opt" + '.' + str(i)
            EnergyCalc(
                self.input, self.input.dict['energy_method'],
                temp_structure, 'run')

            # print out the details at each step to allow restarts
            printx("separation = " + str(dist), self.file)
            printx("energy = " + str(temp_structure.energy), self.file)

            # figure out the next step and check to see if we have converged
            delta_d = self.get_next_step(
                energy_1, temp_structure.energy, delta_d)
            dist += delta_d
            if abs(delta_d) < float(self.input.dict['sep_tolerance']):
                self.print_output(i, dist, delta_d)
                self.structure = temp_structure.copy()
                return dist
            else:
                energy_1 = temp_structure.energy
        printx("Warning, reached end of max_sep_steps without converging. "
               "Last step was " +
               str(delta_d), self.file)
        self.print_output(int(self.input.dict['sep_max_steps']),
                          dist, delta_d)
        self.structure = temp_structure.copy()

        return dist

    def get_next_step(self, energy_1, energy_2, delta_d):
        """
        determine the next step based on the change in energy from the
        previous step. If the energy decreases, then we keep going in
        the same direction. If the energy increases, then we change direction
        and decrease the stepsize.
        """
        delta_energy = energy_2 - energy_1

        # if the change in energy is almost zero, then we are done
        if almost_zero(delta_energy):
            delta_d = 0.0e00
            return delta_d
        elif delta_energy < 0.0:
            delta_d = self.step_size * self.direction
            return delta_d
        else:
            self.direction *= -1
            self.step_size *= 0.7
            delta_d = self.step_size * self.direction
            return delta_d

    def print_output(self, step, sep, change):
        """
        Print out the final results of the separtion optimization
        """
        printx("=====Separation Optimization Completed=====", self.file)
        printx("Number of steps : " + str(step), self.file)
        printx("Final separation : " + str(sep), self.file)
        printx("Final change in separation : " + str(change), self.file)
        printx("===========================================\n", self.file)
