# copyright Caleb Michael Carlin (2018)
# Released uder Lesser Gnu Public License (LGPL)
# See LICENSE file for details.

import random
from math import pi
import numpy as np
import re
import sys
from ase import Atoms
from ase.io import write
from .utilities import InterfaceConfigStorage as ICS
from .utilities import write_xyz, printx
from traceback import print_exc
from .energy_calculators import EnergyCalc, calculate_binding_energy


class ConstrainedSearch(object):
    """
    Class for performing a Markov Chain random walk search of twist
    rotations and/or translations parallel to the interface.


    Paramters:

    interface: InterfaceSupercell object
        object for generating and manipulating interfaces
    input: InputReader object
        object with read in keywords
    initial: InterfaceConfigStorage object
        object containing the interface to be manipulated
    surf_a: tuple
        tuple of Miller Indices for the top crystal. Used for generating
        filenames.
    surf_b: tuple
        tuple of Miller Indices for the bottom crystal. Used for generating
        filenames.
    """

    def __init__(self, interface, input, initial, surf_a, surf_b):
        self.interface = interface
        self.input = input
        self.max_rotation = 2.0 * pi
        self.max_translation = 2.0
        self.lowest_energy_conf = initial.copy()
        self.size_reject = 0
        self.energy_reject = 0
        self.energy_accept = 0
        self.initial = initial
        self.surf_a = surf_a
        self.surf_b = surf_b
        self.file = self.input.dict['output_file']

    def generate_step(self):
        """
        Generate all steps by iterating over all combinations of surfaces
        and a number of translations and rotations for each surface
        combination. The restart option can be used to continue the MC
        search or reload a step in a previous MC search by setting
        'number_of_steps' to 0.
        """
        temporary_store = self.lowest_energy_conf.copy()
        energy_old = temporary_store.energy
        translate_array = [0.0, 0.0, 0.0]
        random.seed()
        # Reload a specific configuration
        if (self.input.dict['mc_restart'] != 'False'):
            rotation = float(self.input.dict['mc_rotate']) * pi / 180.0
            temporary_store.unit_cell_a = self.interface.rotate_cell(
                temporary_store.unit_cell_a, rotation)
            translate = [float(self.input.dict['mc_translate_x']), float(
                self.input.dict['mc_translate_y']), 0.0]
            temporary_store.unit_cell_a = self.interface.translate_cell(
                temporary_store.unit_cell_a, translate)
            self.interface.cut_cell_a = temporary_store.unit_cell_a.copy()
            try:
                self.interface.generate_interface()
            except Exception as err:
                [printx(x, self.file) for x in err.args]
                if self.input.dict['print_debug'] != 'False':
                    print_exc()
                printx(
                    "Could not create structure at restart coordinates. "
                    "Closing out\n")
                sys.exit(1)
            temporary_store.atom = self.interface.interface.copy()
            self.lowest_energy_conf = temporary_store.copy()
            temporary_store.file_name = (temporary_store.file_name +
                                         '.mc_restart')
            write_xyz(temporary_store)

        # generate name used to track results for this interface
        file_name = (str(self.surf_a[0]) + '_' + str(self.surf_a[1]) +
                     '_' + str(self.surf_a[2]) + '_' +
                     str(self.surf_b[0]) + '_' + str(self.surf_b[1]) +
                     '_' + str(self.surf_b[2]))
        # Create a file to track the Markov steps and the simplest steps need
        # to reproduce each configuration
        file_2 = (self.initial.file_name + '.MC.log')
        printx('==========MC steps==========', file_2)
        printx(file_name, file_2)
        printx('============================', file_2)
        mc_track = [[0.0, 0.0, 0.0], [0.0]]
        # Initialize the storage for this set of indices
        for i in range(int(self.input.dict['number_of_steps'])):

            # TorR is translation or rotation
            if (self.input.dict['markov_type'] == '2'):
                TorR = int(random.random() * 2.0)
            elif (self.input.dict['markov_type'] == '1'):
                TorR = 1
            else:
                TorR = 0

            # Translation step
            if (TorR):
                translate_array[0] = random.random(
                ) * self.max_translation - 1.0
                translate_array[1] = random.random(
                ) * self.max_translation - 1.0
                self.interface.cut_cell_a = self.interface.translate_cell(
                    temporary_store.unit_cell_a, translate_array)
                mc_track[0][0] += translate_array[0]
                mc_track[0][1] += translate_array[1]
                printx('step = ' + str(i) + '  Translation : (' +
                       str(translate_array[0]) + ',' + str(translate_array[1])
                       + ')', file_2)
                self.print_mc_coordinates(mc_track, file_2)
            # Rotation step
            else:
                rotation = random.random() * self.max_rotation - pi
                self.interface.cut_cell_a = self.interface.rotate_cell(
                    temporary_store.unit_cell_a, rotation)
                mc_track[1][0] += rotation * 180.0 / pi
                printx('step = ' + str(i) + '  Rotation : (' +
                       str(rotation * 180.0 / pi) + ')', file_2)
                self.print_mc_coordinates(mc_track, file_2)

            # Build a new supercell/supercell interface and calculate energy
            try:
                self.interface.generate_interface()
            except Exception as err:
                [printx(x, self.file) for x in err.args]
                if self.input.dict['print_debug'] != 'False':
                    print_exc()
                self.size_reject += 1
                continue
            # if successfully build interface, then store it to temporary_store
            temporary_store.atom = self.interface.interface
            temporary_store.file_name = (self.input.dict['project_name']
                                         + '.' + str(i) + '.' + file_name)
            temporary_store.step = i
            temporary_store.unit_cell_a = self.interface.cut_cell_a
            EnergyCalc(
                self.input, self.input.dict['energy_method'],
                temporary_store, 'run')
            if self.input.dict['calculate_binding_energy'] != 'False':
                temporary_store.energy = calculate_binding_energy(
                    temporary_store, self.input)
            printx('new energy')
            printx(str(temporary_store.energy))

            printx("energy_old =" + str(energy_old))
            printx("energy new = " + str(temporary_store.energy))
        # If first step for this interface, populate configuration storage
            if (self.lowest_energy_conf.energy is None):
                printx("Error occured in constrained search")
        # If the new configuration reduces the energy, we keep it saved
            elif (temporary_store.energy < energy_old):
                energy_old = temporary_store.energy
                self.lowest_energy_conf = temporary_store.copy()
                self.energy_accept += 1
            else:
                self.energy_reject += 1
        self.print_counters()
        # output a restart coordinate file at the end if requested
        if (self.input.dict['mc_write_coord_file'] != 'False'):
            write_xyz(self.lowest_energy_conf)

        return self.lowest_energy_conf

    def print_counters(self):
        """print out statistics for the specific interface"""
        file = self.input.dict['output_file']
        printx('\n', file)
        printx('======Search Results======', file)
        printx('for interface = ' + str(self.surf_a) +
               '/' + str(self.surf_b), file)
        printx('Energy at first step = ' +
               str(self.initial.energy), file)
        printx('Number of steps rejected for size = ' +
               str(self.size_reject), file)
        printx('Number of steps with higher energy = ' +
               str(self.energy_reject), file)
        printx('Number of steps with lower energy = ' +
               str(self.energy_accept), file)
        printx('Energy of lowest energy config = ' +
               str(self.lowest_energy_conf.energy), file)
        printx('File name of lowest energy config = ' +
               str(self.lowest_energy_conf.file_name), file)
        printx('==========================' + '\n', file)

    def print_mc_coordinates(self, points, file):
        """
        print the translation and rotation needed to reproduce a given state
        from the Markov Chain
        """
        rot = points[1][0] % 360.0
        printx("Translate X : %5f  Y : %5f  Rotate %5f \n" %
               (points[0][0], points[0][1], rot), file)
        printx("---------------------\n", file)

        return
