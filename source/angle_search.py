# copyright Caleb Michael Carlin (2018)
# Released uder Lesser Gnu Public License (LGPL)
# See LICENSE file for details.

from ase import Atoms
from .energy_calculators import EnergyCalc, calculate_binding_energy
from ase.io import read, write, Trajectory
from .utilities import InterfaceConfigStorage as ICS
from .utilities import write_xyz, printx
import sys
import os
from math import pi
from .separation_optimizer import SeparationOpt
from traceback import print_exc


class AngleSearch(object):
    """
    Class for rotating one side of the interface around the axis parallel
    to the interface.

    Paramters:

    input: InputReader object
        object with read in keywords
    interface: InterfaceSupercell object
        object for generating and manipulating interfaces
    ics_initial: InterfaceConfigStorage object
        object containing the interface to be manipulated
    surf_a: tuple
        tuple of Miller Indices for the top crystal. Used for generating
        filenames.
    surf_b: tuple
        tuple of Miller Indices for the bottom crystal. Used for generating
        filenames.
    """

    def __init__(self, input, interface, ics_initial, surf_a, surf_b):
        self.input = input
        self.interface = interface
        self.initial = ics_initial
        self.surf_a = surf_a
        self.surf_b = surf_b
        self.file = self.input.dict['output_file']
        self.is_ras = (int(self.input.dict['ras_depth']) != 1)

    def angle_search(self):
        """Main driver for twist angle search."""
        temp_storage = ICS()
        lowest_energy = self.initial.copy()

        # backup input values
        stepsize_backup = self.input.dict['angles_stepsize']
        angle_list_backup = self.input.dict['angles_list']
        number_of_angles_backup = self.input.dict['number_of_angles']
        start_angle_backup = self.input.dict['starting_angle']
        # turn off energy calculations if using ras
        if self.is_ras:
            self.input.dict['angle_calculate_energy'] = 'False'
            self.input.dict['angle_write_energy_file'] = 'False'
            ras_array = [[] for y in range(2)]

        # check to read in angles from a previous RAD run
        if self.input.dict['read_in_ras_file'] is not None:
            try:
                self.input.dict['angles_list'] = self.read_in_angle_ras_file()
            except Exception as err:
                [printx(x, self.file) for x in err.args]
                if self.input.dict['print_debug'] != 'False':
                    print_exc()
                sys.exit(1)

        # set up main output with instant printing (no buffer)
        printx('Beginning twist angle search\n', self.file)
        printx('Surface Indices : ' + str(self.surf_a) +
               '/' + str(self.surf_b) + '\n', self.file)

        # If iterating over a number angles, we need to populate the
        # angles_to_gen list
        self.determine_angles()
        # iterate over angles
        lowest_energy, ras_array = self.loop_over_angles(
            temp_storage, lowest_energy)

        # set up initial arrays for ras
        if (self.input.dict['ras_all_angles'] != 'False'):
            self.input.dict['ras_factor'] = len(ras_array[0])

        # pull angles to use for RAS
        if (len(ras_array[0]) < int(self.input.dict['ras_factor'])
                or not(self.is_ras)):
            printx("Number of acceptable interfaces found is less than"
                   "'ras_factor'")
            printx("Suggest either increasing 'max_atoms' or reducing"
                   "'ras_factor'")
            self.input.dict['angles_stepsize'] = stepsize_backup
            self.input.dict['angles_list'] = angle_list_backup
            self.input.dict['number_of_angles'] = number_of_angles_backup
            self.input.dict['starting_angle'] = start_angle_backup
            return
        else:
            print('xgeox')
            print(ras_array)
            print(self.is_ras)
            while len(ras_array[0]) > int(self.input.dict['ras_factor']):
                largest = max(
                    xrange(len(ras_array[1])), key=ras_array[1].__getitem__)
                ras_array[1].pop(largest)
                ras_array[0].pop(largest)

        # If using RAS, begin the repeated searches for RAS level 2 through
        # 'ras_depth'
        for j in range(int(self.input.dict['ras_depth']) - 1):
            printx('\n================\nBeginning ras search ' +
                   str(j + 2) + '\n================\n', self.file)
            step = float(self.input.dict['angles_stepsize'])
            self.input.dict['angles_stepsize'] = step / 10.0
            self.input.dict['number_of_angles'] = 20
            ras_temp = [[] for y in range(2)]
            for k in range(len(ras_array[0])):
                self.input.dict['starting_angle'] = ras_array[0][k] - step
                # If iterating over a number angles, we need to populate the
                # angles_to_gen list
                self.input.dict['angles_list'] = None
                self.determine_angles()
                # iterate over angles
                replace_array = [[] for y in range(2)]
                lowest_energy, replace_array = self.loop_over_angles(
                    temp_storage, lowest_energy)
                try:
                    smallest = min(
                        xrange(len(replace_array[1])),
                        key=replace_array[1].__getitem__)
                    ras_temp[0].append(replace_array[0].pop(smallest))
                    ras_temp[1].append(replace_array[1].pop(smallest))
                except Exception as err:
                    [printx(x) for x in err.args]
                    printx("Error : angle search returned no successful "
                           "interfaces")
                    ras_temp[0].append(ras_array[0][k])
                    ras_temp[1].append(ras_array[1][k])
            ras_array = ras_temp

        # Setup energy vs angle file
        if (self.input.dict['angle_write_energy_file'] == 'True'):
            angle_energy = open(
                self.input.dict['project_name'] + '.a_e.log', 'a', 0)
            angle_energy.write('===================\n\n')
            angle_energy.close()
        # If using RAS, resets the input values and recovers the lowest energy
        # interface located to pass back.
        if (int(self.input.dict['ras_depth']) != 1):
            if (self.input.dict['ras_energy'] != 'False'):
                self.input.dict['angles_list'] = ras_array[0]
                self.input.dict['angle_calculate_energy'] = 'True'
                self.input.dict['angle_write_energy_file'] = 'True'
                lowest_energy, replace_array = self.loop_over_angles(
                    temp_storage, lowest_energy)
            printx('\n ====== RAS results ====== ', self.file)
            for l in range(int(self.input.dict['ras_factor'])):
                printx(
                    'angle = ' + str(ras_array[0][l]) + ' size = ' +
                    str(ras_array[1][l]), self.file)
            printx('\n', self.file)
        # reset values in case twist_search is called more than once
        self.input.dict['angles_stepsize'] = stepsize_backup
        self.input.dict['angles_list'] = angle_list_backup
        self.input.dict['number_of_angles'] = number_of_angles_backup
        self.input.dict['starting_angle'] = start_angle_backup

        if (self.input.dict['angle_return_initial'] != 'True'):
            return lowest_energy
        else:
            return self.initial

    def determine_angles(self):
        """
        if angles_list is empty, populate it using
        angles_to_iter and number_of_angles
        """
        if (self.input.dict['angles_list'] is None):
            if (self.input.dict['angles_stepsize'] is None or
                    self.input.dict['number_of_angles'] is None):
                printx('Error: additional information for angles needed in'
                       'the input file.')
                sys.exit(1)
            temp_array = []
            angle = float(self.input.dict['starting_angle'])
            for i in range(int(self.input.dict['number_of_angles'])):
                angle += float(self.input.dict['angles_stepsize'])
                temp_array.append(angle)
            self.input.dict['angles_list'] = temp_array

    def loop_over_angles(self, temp_storage, lowest_energy):

        backup_initial = lowest_energy.copy()
        temp_angles = [[] for y in range(2)]
        # setup energy output file header
        if (self.input.dict['angle_write_energy_file'] == 'True'):
            angle_energy = open(
                self.input.dict['project_name'] + '.a_e.log', 'a')
            angle_energy.write("===================\nSurface : " +
                               str(self.surf_a) + '/' +
                               str(self.surf_b) + '\n')
            angle_energy.write(
                "===================\nAngle   Energy   Atoms\n\n")
            angle_energy.close()

        # loop over angles
        for j in self.input.dict['angles_list']:
            printx('===========\n', self.file)
            printx('beginning rotation of ' +
                   str(j) + ' degrees:\n', self.file)
            Backup_atom = self.initial.unit_cell_a.copy()
            radians = pi * j / 180.0
            self.interface.cut_cell_a = self.interface.rotate_cell(
                self.initial.unit_cell_a, radians).copy()
            try:
                self.interface.generate_interface()
            except Exception as err:
                [printx(x, self.file) for x in err.args]
                if self.input.dict['print_debug'] != 'False':
                    print_exc()
                self.initial.unit_cell_a = Backup_atom.copy()
                continue
            printx('printing output: angle = ' + str(j) + ' atoms = ' +
                   str(len(self.interface.interface)) + '\n', self.file)
            temp_angles[0].append(j)
            temp_angles[1].append(len(self.interface.interface))
            # calculate energy and return lowest energy configuration
            if (self.input.dict['angle_calculate_energy'] != 'False'):
                temp_storage.atom = self.interface.interface.copy()
                surf_index_a = ''.join(str(x) for x in self.surf_a)
                surf_index_b = ''.join(str(y) for y in self.surf_b)
                temp_storage.file_name = (self.input.dict['project_name'] +
                                          '.' + str(
                    j) + '.' + surf_index_a + surf_index_b)
                temp_storage.step = j
                temp_storage.unit_cell_a = self.interface.cut_cell_a.copy()
                temp_storage.unit_cell_b = self.interface.cut_cell_b.copy()

                # Check to see if this combination of surfaces and angle has
                # already been calculated to save on work done
                Exist = os.path.isfile(temp_storage.file_name + '.out')
                if Exist:
                    printx("Duplicate angle for energy calculation detected."
                           "Skipping energy calculation")
                    temp_storage.energy = 1.0e10
                else:
                    # print out a traj file that can be used for ET calcs
                    # without calculating energy again
                    if (self.input.dict['angle_write_restart'] != 'False'):
                        temp_name = temp_storage.file_name + '.traj'
                        traj = Trajectory(filename=(temp_name),
                                          mode='w', atoms=temp_storage.atom)
                        traj.write()
                        traj.close()
                    try:
                        EnergyCalc(
                            self.input, self.input.dict['energy_method'],
                            temp_storage, 'run')
                    except Exception as err:
                        [printx(x) for x in err.args]
                        printx(
                            'An error occured when calculating the energy '
                            'for rotation ' + str(j) + ' degrees')
                        continue
                    if self.input.dict['angle_optimize_separation'] != 'False':
                        temp_storage = self.optimize_separation(temp_storage)
                    try:
                        if (self.input.dict['calculate_binding_energy']
                                != 'False'):
                            temp_storage.energy = calculate_binding_energy(
                                temp_storage, self.input)
                    except Exception as err:
                        [printx(x) for x in err.args]
                        printx(
                            'An error occured when calculating the energy for '
                            'rotation ' + str(j) + ' degrees')
                        continue

                # determine if this is lowest energy  structure
                if (temp_storage.energy < lowest_energy.energy):
                    lowest_energy = temp_storage.copy()
                if (self.input.dict['angle_write_energy_file'] == 'True'
                        and not Exist):
                    angle_energy = open(
                        self.input.dict['project_name'] + '.a_e.log', 'a', 0)
                    angle_energy.write(str(j) + '   ' +
                                       str(temp_storage.energy)
                                       + '    ' +
                                       str(len(self.interface.interface))
                                       + '\n')
                    angle_energy.close()
            # set up information in case we didn't run the energy calculation
            else:
                temp_storage.atom = self.interface.interface.copy()
                surf_index_a = ''.join(str(x) for x in self.surf_a)
                surf_index_b = ''.join(str(y) for y in self.surf_b)
                temp_storage.file_name = (self.input.dict['project_name'] +
                                          '.' + str(
                    j) + '.' + surf_index_a + surf_index_b)
                temp_storage.step = j
                temp_storage.unit_cell_a = self.interface.cut_cell_a.copy()
                temp_storage.unit_cell_b = self.interface.cut_cell_b.copy()
                lowest_energy = temp_storage.copy()
            # write the xyz file
            if (self.input.dict['angle_write_coord_file'] == 'True'):
                write_storage = ICS()
                if self.input.dict['angle_calculate_energy'] != 'False':
                    write_storage.atom = temp_storage.atom.copy()
                else:
                    write_storage.atom = self.interface.interface.copy()
                write_storage.file_name = (self.input.dict['project_name']
                                           + '.' + str(j))
                write_xyz(write_storage)
            self.initial.unit_cell_a = Backup_atom.copy()

        printx("\nLowest Energy configuration correspond to angle " +
               str(lowest_energy.step) + '\n', self.file)

        if (self.input.dict['angle_return_initial'] != 'True'):
            return lowest_energy, temp_angles
        else:
            return backup_initial, temp_angles

    def read_in_angle_ras_file(self):
        """
        read in an output file from a previously completed run to get
        the set of angles to use.
        The file read in is assumed to be the 'output_file' from a
        previous RAS run.
        """

        file_name = self.input.dict['read_in_ras_file']
        try:
            with open(file_name) as f:
                lines = f.readlines()
        except Exception as err:
            raise Exception("Error couldn't open old output file for reading")

        angles = []
        switch = False
        for i in lines:
            if switch and i != '\n':
                line = i.split()
                if line[0] == 'angle':
                    angles.append(float(line[2]))
                else:
                    continue
            if i == " ====== RAS results ====== \n":
                switch = True
        if not switch:
            raise Exception("Error, no RAS results in output file")

        return angles

    def optimize_separation(self, ics):
        """
        run the separation optimizer to allow the structure to relax at
        each new angle
        """

        backup_binding_energy = self.input.dict['calculate_binding_energy']
        self.input.dict['calculate_binding_energy'] = 'False'

        SeparationOpt(ics, self.input, self.surf_a, self.surf_b)
        self.input.dict['calculate_binding_energy'] = backup_binding_energy

        return ics.copy()
