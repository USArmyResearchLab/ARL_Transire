# copyright Caleb Michael Carlin (2018)
# Released uder Lesser Gnu Public License (LGPL)
# See LICENSE file for details.

import re
import sys
import numpy as np
from .utilities import printx, almost_zero


class ReadInput(object):

    def __init__(self):
        """ master list of all inputs with default values"""
        self.dict = {'crys_a_file': None,
                     'crys_a_surface': None,
                     'crys_a_layers': 3,
                     'crys_b_file': None,
                     'crys_b_surface': None,
                     'crys_b_layers': 3,
                     'cp2k_input': None,
                     'separation': 0.0,
                     'max_mpi_processes': 32,
                     'atoms_per_process': 13,
                     'working_directory': "./",
                     'project_name': "default_name",
                     'range_surface_h_a': None,
                     'range_surface_k_a': None,
                     'range_surface_l_a': None,
                     'range_surface_h_b': None,
                     'range_surface_k_b': None,
                     'range_surface_l_b': None,
                     'number_of_steps': 25,
                     'max_atoms': 75000,
                     'output_file': "transire.out",
                     'search_list': None,
                     'angles_list': None,
                     'angles_stepsize': None,
                     'number_of_angles': None,
                     'tolerance': 5.0e-2,
                     'number_of_layers_a': 2,
                     'number_of_layers_b': 2,
                     'perform_ET': 'False',
                     'restart_path': "./",
                     'ET_restart': 'False',
                     'calculate_initial_energy': 'True',
                     'angle_write_energy_file': 'True',
                     'angle_write_coord_file': 'True',
                     'angle_calculate_energy': 'True',
                     'angle_write_restart': 'True',
                     'exclude_coupling': 'False',
                     'energy_levels_ET': None,
                     'restart_file': 'ET_restart_coord',
                     'insert_file': None,
                     'calc_insert_energy': 'False',
                     'starting_angle': 0.0,
                     'sep_guess': 0.5,
                     'sep_max_steps': 25,
                     'sep_tolerance': 1e-5,
                     'sep_initial_step': 0.1,
                     'markov_type': 2,
                     'print_debug': 'False',
                     'remove_duplicates': 'False',
                     'ras_factor': 5,
                     'ras_depth': 1,
                     'ras_energy': 'True',
                     'ras_all_angles': 'False',
                     'lammps_input': None,
                     'energy_method': 'cp2k',
                     'surface_search': 'False',
                     'ignore_initial_structure': 'False',
                     'mc_translate_x': 0.0,
                     'mc_translate_y': 0.0,
                     'mc_rotate': 0.0,
                     'mc_restart': 'False',
                     'mc_write_coord_file': 'False',
                     'angle_return_initial': 'False',
                     'insert_vacuum_below': 0.5,
                     'insert_vacuum_above': 0.0,
                     'read_in_ras_file': None,
                     'calculate_binding_energy': 'True',
                     'orthonormal_overlap': 'False',
                     'rotate_insert_x': None,
                     'rotate_insert_y': None,
                     'rotate_insert_z': None,
                     'angle_optimize_separation': 'True',
                     'full_periodicity': 'False',
                     'z_axis_vacuum' : '0.0',
                     'et_method': 'ASE'}

        # inputs to convert to tuples
        self.tuples = ['crys_a_surface',
                       'crys_b_surface',
                       'range_surface_h_a',
                       'range_surface_k_a',
                       'range_surface_l_a',
                       'range_surface_h_b',
                       'range_surface_k_b',
                       'range_surface_l_b',
                       ]

        # inputs to convert to a list of floats
        self.lists = ['angles_list',
                      'search_list',
                      'energy_levels_ET'
                      ]

        # inputs that are required
        self.required = ['crys_a_file',
                         'crys_a_surface',
                         'crys_b_file',
                         'crys_b_surface']
        # keywords with either True or False strings
        self.logicals = ['calculate_initial_energy',
                         'perform_ET',
                         'ET_restart',
                         'angle_write_energy_file',
                         'angle_write_coord_file',
                         'angle_calculate_energy',
                         'exclude_coupling',
                         'calc_insert_energy',
                         'print_debug',
                         'remove_duplicates',
                         'ras_energy',
                         'ras_all_angles',
                         'surface_search',
                         'angle_write_restart',
                         'ignore_initial_structure',
                         'mc_restart',
                         'mc_write_coord_file',
                         'angle_return_initial',
                         'calculate_binding_energy',
                         'orthonormal_overlap',
                         'full_periodicity',
                         'angle_optimize_separation']

    def read_input(self, file_name):
        """
        Scans through the dict keywords and searches for each in
        the input file.  The line is split into pieces, taking the
        value to be anything after the '='

        Parameter:

        file_name: string
            file name of Transire input file with path if not in
            current directory.
        """
        in_file = open(file_name).read()
        for i in self.dict:
            found = re.search('[#]*' + i + '[A-Za-z0-9\t_ .=\-~/#]*', in_file)
            if (found is not None):
                line = found.group()
                if (line.startswith('#')):
                    continue
                value = re.search('(?<==)[A-Za-z0-9\t_ .\-~/]*', line)
                if (value is not None):
                    value = value.group()
                    value = value.strip()
                    self.dict[i] = value

        # convert to tuples and lists
        for j in self.tuples:
            if (self.dict[j] is not None):
                self.dict[j] = self.string_to_tuple(self.dict[j])
        for k in self.lists:
            if (self.dict[k] is not None):
                self.dict[k] = self.string_to_list(self.dict[k])

        # check certain values to see if they are in the acceptable range
        self.error_check()

        if (self.dict['print_debug'] != 'False'):
            for i in self.dict:
                printx(i + " = " + str(self.dict[i]))

        return

    def string_to_tuple(self, string_in):
        """converts a string of space separated integers into a tuple"""
        tuple_out = []
        split_string = string_in.split()
        for line in split_string:
            tuple_out.append(int(line))
        tuple_out = tuple(tuple_out)

        return tuple_out

    def string_to_list(self, string_in):
        """converts a string to a list of floats"""
        string_out = string_in.split()
        for i in range(len(string_out)):
            string_out[i] = float(string_out[i])

        return string_out

    def error_check(self):
        if (int(self.dict['crys_a_layers']) <= 0):
            printx('Error: crys_a_layers must be 1 or greater')
            sys.ext(1)

        if (int(self.dict['crys_b_layers']) <= 0):
            printx('Error: crys_b_layers must be 1 or greater')
            sys.exit(1)

        if (float(self.dict['separation']) < 0.0):
            printx('Error: separation can not be negative')
            sys.exit(1)

        if (int(self.dict['max_atoms']) <= 0):
            printx('Error: max_atoms must be greater than 0')
            sys.exit(1)

        if (int(self.dict['max_mpi_processes']) < 1):
            printx('Error: max_mpi_processes must be 1 or greater')
            sys.exit(1)

        if (int(self.dict['atoms_per_process']) < 1):
            printx('Error: atoms_per_process must be 1 or greater')
            sys.exit(1)

        if (np.mod(int(self.dict['number_of_layers_a']), 2) != 0):
            printx('Error: number_of_layers_a must be a multiple of 2')
            sys.exit(1)

        if (np.mod(int(self.dict['number_of_layers_b']), 2) != 0):
            printx('Error: number_of_layers_b must be a multiple of 2')
            sys.exit(1)

        if (self.dict['energy_levels_ET'] is not None):
            if (len(self.dict['energy_levels_ET']) != 3):
                printx('Error: energy_levels_ET must have all 3'
                       ' values specified')
                sys.exit(1)

        if almost_zero(float(self.dict['separation']) -
                       float(self.dict['sep_guess'])):
            printx(
                "Error: separation and sep_guess can't be the same."
                "  Nudging sep_guess")
            self.dict['sep_guess'] = float(self.dict['sep_guess']) + 0.1

        energy_options = ['lammps', 'cp2k']
        match = 0
        for k in energy_options:
            if (self.dict['energy_method'] == k):
                match = 1
        if (match != 1):
            printx("Error: energy_method is not a valid option.")
            sys.exit(1)

        for i in self.logicals:
            if (self.dict[i] != 'True' and self.dict[i] != 'False'):
                printx('Error: {0} must be either True or False'.format(i))
                sys.exit(1)
        return
