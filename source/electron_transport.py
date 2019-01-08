# copyright Caleb Michael Carlin (2018)
# Released uder Lesser Gnu Public License (LGPL)
# See LICENSE file for details.

import re
import ase
from ase import Atoms
from ase.transport.calculators import TransportCalculator as TCalc
import sys
import numpy as np
from ase.io import write, read
from ase.io import Trajectory
import warnings
import gc
from .utilities import printx
from .interface import InterfaceSupercell

class ElectronTransport(object):
    """
    Class for setting up and performing an electron transport calculation using the
    Non-Equilibrium Green's Function method in ASE.

    Parameters :

    input: InputReader object
        object with read in keywords
    interface: InterfaceSupercell object
        object for generating and manipulating interfaces
    restart: Boolean
        switch for whether this is the continuation of a previous
        electron transport calculation
    """

    def __init__(self, input, interface, restart):
        self.basis_sizes = {}
        self.number_atoms = {}
        self.array_dict = {}
        self.input = input
        self.interface = interface
        self.restart = restart
        # g_o = 2e^2/h or the conversion from transmission to conductance
        # in siemens
        self.g_o = 7.7480917310e-5
        self.file = self.input.dict['output_file']
	if (self.input.dict['et_method'] == 'slymer'):
	#	import pyTightWrapper as ptw
		import pTW as ptw


    def calculate_ET(self, unit_cell_a, unit_cell_b):
        """
        Set up ET calculator and calculate transmission of
        interface_config_structure
        """
        # either write or read the restart file depending on if we performed
        # the initial energy or interface search
        if (self.restart is not False):
            traj = Trajectory(
                self.input.dict['restart_path'] + '/' +
                self.input.dict['restart_file'] + '.traj')
            self.interface.atom = traj[0]
            self.interface.file_name = self.input.dict['restart_file']
            traj.close()
        else:
            temp = Trajectory(filename='ET_restart_coord.traj',
                              mode='w', atoms=self.interface.atom)
            temp.write()
            temp.close()
        # read the output file to get the basis set sizes
        scatter_basis, scatter_number = self.get_basis_sizes(
            self.input.dict['restart_path'] +
            self.interface.file_name + '.out')
        # determine the important dimensions
        llead_dimens, llead_apl, llead_number = self.determine_dimensions(
            unit_cell_a, 1)
        rlead_dimens, rlead_apl, rlead_number = self.determine_dimensions(
            unit_cell_b, 2)
        # gather what we need from the last output of the whole system energy
        # calculation and the initial unit cells
        basis_a = 0
        basis_b = 0
        for i in scatter_basis:
            basis_a += scatter_basis[i] * llead_number[i]
            basis_b += scatter_basis[i] * rlead_number[i]

        Total_basis = 0
        for i in scatter_basis:
            Total_basis += scatter_basis[i] * scatter_number[i]
        N1 = abs(int((basis_a / llead_dimens[2, 2]) *
                 int(self.input.dict['number_of_layers_a'])))
        N2 = abs(int((basis_b / rlead_dimens[2, 2]) *
                 int(self.input.dict['number_of_layers_b'])))
        N = int(Total_basis - N1 - N2)

        # take a moment to print out the filename being used
        printx('\n' + self.interface.file_name + '\n', self.file)

        collected = gc.collect()
        printx("Garbage collector: collected %d objects." % (collected))

        if (self.input.dict['et_method'] is not 'ASE'):
            Hamil, Overlap = self.get_whole_matrix(
                Total_basis, self.input.dict['restart_path'] +
                self.interface.file_name + '-hamiltonian-1_0.Log')
        else:    
            scatter_hamil = self.get_matrices(
                N, N1, N2, self.input.dict['restart_path'] +
                self.interface.file_name + '-hamiltonian-1_0.Log')

        if (self.input.dict['et_method'] is not 'ASE'):
            points = self.input.dict['energy_levels_ET']
            starting = points[0]
            step_size = points[2]
            steps = int((points[1]-points[0])/points[2])
            Transmission = ptw.call_slymer(Hamil,Overlap, N1, N2, steps,
                starting, step_size)
            printx("Transmission values \n", self.file)
            et_energy = starting
            for output in Transmission:
                printx(str(et_energy)+'   '+str(output),self.file)
                et_energy += step_size
        else:
            ase_et(scatter_hamil)

        return


    def ase_et(self, scatter_hamil):
        
        # Set up options if we want to call ASE with special cases
        if self.input.dict['orthonormal_overlap'] != 'False':
            scatter_hamil['s'] = None
            scatter_hamil['s1'] = None
            scatter_hamil['s2'] = None
        # initialize the calculator
        ETran = TCalc(h=scatter_hamil['h'], h1=scatter_hamil['h1'],
                      h2=scatter_hamil['h2'],
                      s=scatter_hamil['s'], s1=scatter_hamil['s1'],
                      s2=scatter_hamil['s2']
                      )

        # Disable reading in the off diagonal matrices
        if (self.input.dict['exclude_coupling'] == 'False'):
            ETran.set(hc1=scatter_hamil['hc1'],
                      hc2=scatter_hamil['hc2'],
                      sc1=scatter_hamil['sc1'], sc2=scatter_hamil['sc2'])
        # either calculate the electron transport at the Fermi level or a
        # range of values relative to the Fermi level
        if (self.input.dict['energy_levels_ET'] is None):
            temp_form = [(0.0)]
            ETran.set(energies=[(0.0)])
        else:
            temp_form = self.input.dict['energy_levels_ET']
            ETran.set(energies=np.arange(
                temp_form[0], temp_form[1], temp_form[2]))

        # ase has some warnings that we can't do anything about,
        # so suppress the warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            T = ETran.get_transmission()
        # print results
        if (self.input.dict['energy_levels_ET'] is not None):
            j = 0
            for i in np.arange(temp_form[0], temp_form[1], temp_form[2]):
                printx("Energy = {0:.3f} Transmission = {1}".format(
                    i, T[j]), self.file)
                j += 1
            self.print_conductance(T, temp_form[2], temp_form[0],
                                   temp_form[1])
        else:
            printx("Energy = 0.0  Transmission = {0}".format(
                T[0]), self.file)
        printx("\n", self.file)

        return

    def get_basis_sizes(self, file_name):
        """
        read output file to get number of basis functions for each element type
        """
        number_atoms = {}
        basis_sizes = {}
        with open(file_name) as output:
            data = output.read()
            found = re.findall('Atomic kind:' + '[A-Za-z0-9\t :]*', data)
            found2 = re.findall(
                'Number of spherical basis functions:' + '[0-9 ]*', data)
            j = 0
            for i in found:
                line = i.split()
                line2 = found2[j].split()
                number_atoms[line[2]] = int(line[6])
                basis_sizes[line[2]] = int(line2[5])
                j += 1

        return basis_sizes, number_atoms

    def get_whole_matrix(self,dimen, file):
        """
        read in the Hamiltonian and Overlap matrices.
        Because Slymer divides things internally, we
        don't do that here.
        """
        temp = []
        section = 0
        H = np.ndarray(shape=(dimen, dimen), dtype=np.float64)
        S = np.ndarray(shape=(dimen, dimen), dtype=np.float64)

        with open(file) as f:
            for line in f:
                temp = line.split()
                length = len(temp)
                if length == 4:
                    y += 4
                    x = 0
                elif length > 4:
                    if section == 2:
                        for val in range(length - 4):
                            H[x,y+val] = float(temp[val+4])
                    if section == 1:
                        for val in range(length - 4):
                            S[x,y+val] = float(temp[val+4])
                    x += 1
                elif length == 3:
                    section = 2  # HAMILTONIAN section
                    y = -4
                elif length == 2:
                    section = 1  # OVERLAP section
                    y = -4

        return H, S

    def get_matrices(self, N, N1, N2, file):
        """
        Read matrices from previous output into
        left, right, and scatter regions
        """
        # should be rewritten to streamline
        temp = []
        section = 0
        NA = N + N1
        NB = NA + N2
        NC = N1 // 2
        ND = N2 // 2
        # take a moment to print out useful information for debugging
        printx('===Extracting the relavant matrices===', self.file)
        printx('Important dimensions (number of basis functions):',
               self.file)
        printx('scatter region  = ' + str(N), self.file)
        printx('left lead = ' + str(N1), self.file)
        printx('right lead = ' + str(N2) + '\n', self.file)
        h = np.ndarray(shape=(N, N), dtype=complex)
        h1 = np.ndarray(shape=(N1, N1), dtype=complex)
        h2 = np.ndarray(shape=(N2, N2), dtype=complex)
        hc1 = np.ndarray(shape=(NC, N), dtype=complex)
        hc2 = np.ndarray(shape=(ND, N), dtype=complex)
        s = np.ndarray(shape=(N, N), dtype=complex)
        s1 = np.ndarray(shape=(N1, N1), dtype=complex)
        s2 = np.ndarray(shape=(N2, N2), dtype=complex)
        sc1 = np.ndarray(shape=(NC, N), dtype=complex)
        sc2 = np.ndarray(shape=(ND, N), dtype=complex)

        with open(file) as f:
            for line in f:
                temp = line.split()
                length = len(temp)
                if length == 4:
                    y += 4
                elif length > 4:
                    if section == 1:
                        self.overlap(s, s1, s2, sc1, sc2, temp,
                                     y, N1, NA, NB, N, N2)
                    elif section == 2:
                        self.hamiltonian(h, h1, h2, hc1, hc2,
                                         temp, y, N1, NA, NB, N, N2)
                    else:
                        printx("length = {0}  section = {1}".format(
                            length, section))
                        printx("There is an error reading matrices")
                elif length == 3:
                    section = 2  # HAMILTONIAN section
                    y = -4
                elif length == 2:
                    section = 1  # OVERLAP section
                    y = -4
        array_dict = {'h': h, 'h1': h1, 'h2': h2, 'hc1': hc1, 'hc2': hc2,
                      's': s, 's1': s1, 's2': s2, 'sc1': sc1, 'sc2': sc2}

        return array_dict

    def hamiltonian(self, h, h1, h2, hc1, hc2,
                    split_line, y, N1, NA, NB, N, N2):
        """
        Read hamiltonian and splits it into pieces to match ASE ET calc
        """
        silly = N1 // 2 + N1
        silly2 = N2 // 2 + NA
        x = int(split_line[0]) - 1
        line_length = len(split_line)
        # RL
        if (x < N1 and y < N1):
            yval = y
            for counter in range(4, line_length):
                h1[x, yval] = float(split_line[counter]) + 0.j
                yval += 1
                if (yval == N1):
                    break
        # Scatter
        elif (x >= N1 and x < NA and y >= N1 and y < NA):
            a = x - N1
            b = y - N1
            yval = b
            for counter in range(4, line_length):
                h[a, yval] = float(split_line[counter]) + 0.j
                yval += 1
                if (yval == N):
                    break
        # LL
        elif (x >= NA and x < NB and y >= NA and y < NB):
            a = x - NA
            b = y - NA
            yval = b
            for counter in range(4, line_length):
                h2[a, yval] = float(split_line[counter]) + 0.j
                yval += 1
                if (yval == N2):
                    break
        # RL/Scatter coupling
        elif (x >= N1 and x < silly and y < N):
            a = x - N1
            yval = y
            for counter in range(4, line_length):
                hc1[a, yval] = float(split_line[counter]) + 0.j
                yval += 1
                if (yval == N1):
                    break
        # LL/Scatter coupling
        elif (x >= NA and x < silly2 and y >= N1 and y < NA):
            a = x - NA
            b = y - N1
            yval = b
            for counter in range(4, line_length):
                hc2[a, yval] = float(split_line[counter]) + 0.j
                yval += 1
                if (yval == N):
                    break
        return

    def overlap(self, s, s1, s2, sc1, sc2, split_line, y, N1, NA, NB, N, N2):
        """
        Read overlap and splits it into pieces to match ASE ET calc
        """
        silly = N1 // 2 + N1
        silly2 = N2 // 2 + NA
        x = int(split_line[0]) - 1
        line_length = len(split_line)
        # RL
        if (x < N1 and y < N1):
            yval = y
            for counter in range(4, line_length):
                s1[x, yval] = float(split_line[counter]) + 0.j
                yval += 1
                if (yval == N1):
                    break
        # Scatter
        elif (x >= N1 and x < NA and y >= N1 and y < NA):
            a = x - N1
            b = y - N1
            yval = b
            for counter in range(4, line_length):
                s[a, yval] = float(split_line[counter]) + 0.j
                yval += 1
                if (yval == N):
                    break
        # LL
        elif (x >= NA and x < NB and y >= NA and y < NB):
            a = x - NA
            b = y - NA
            yval = b
            for counter in range(4, line_length):
                s2[a, yval] = float(split_line[counter]) + 0.j
                yval += 1
                if (yval == N2):
                    break
        # RL/Scatter coupling
        elif (x >= N1 and x < silly and y < N1):
            a = x - N1
            yval = y
            for counter in range(4, line_length):
                sc1[a, yval] = float(split_line[counter]) + 0.j
                yval += 1
                if (yval == N1):
                    break
        # LL/Scatter coupling
        elif (x >= NA and x < silly2 and y >= N1 and y < NA):
            a = x - NA
            b = y - N1
            yval = b
            for counter in range(4, line_length):
                sc2[a, yval] = float(split_line[counter]) + 0.j
                yval += 1
                if (yval == N):
                    break
        return

    def determine_dimensions(self, unit_cell, tag):
        """
        use interface_cell /unit_cell to get x and y dimensions
        use number of atoms per layer to calculate number of layers.
        """
        interface_cell = self.interface.atom.cell
        # determine the orthorhombic representation calling
        # protect_periodicity from interface.py
        side_cell = self.call_protect_periodicity(unit_cell)
        lead_num = {}

        x_dimens = interface_cell[0, 0] / side_cell[0, 0]
        x_dimens = np.around(x_dimens)
        y_dimens = interface_cell[1, 1] / side_cell[1, 1]
        y_dimens = np.around(y_dimens)

        atom_count = 0
        assign_array = self.interface.atom.get_tags()
        element_array = self.interface.atom.get_chemical_symbols()
        temp = set(element_array)
        for j in temp:
            lead_num[j] = 0
        for i in range(len(assign_array)):
            if (assign_array[i] == tag):
                atom_count += 1
                lead_num[element_array[i]] = lead_num[element_array[i]] + 1
        atoms_per_layer = x_dimens * y_dimens * len(unit_cell)
        layers = atom_count / atoms_per_layer
        P = np.eye(3)
        P[0, 0] = x_dimens
        P[1, 1] = y_dimens
        P[2, 2] = layers

        return P, atoms_per_layer, lead_num

    def call_protect_periodicity(self, unit_cell):
        """We need to call an instance of the intercace builder
        so that we can use protect_periodicity to find the
        orthorhombic representation of the unit cell"""

        temp_IS = InterfaceSupercell(unit_cell,unit_cell,self.input)
        side_cell, throw_away = temp_IS.protect_periodicity(unit_cell)

        return side_cell

    def print_conductance(self, transmission, step_size,
                          low_range, high_range):

        total_conductance = 0
        for i in transmission:
            total_conductance += i * step_size
        total_conductance *= self.g_o

        printx('=============Conductance=============\n', self.file)
        printx("In range from {0} eV to {1} eV relative to Fermi level\n"
               .format(low_range, high_range), self.file)
        printx("Conductance = {0} siemens\n".format(
            total_conductance), self.file)
        printx('=====================================\n', self.file)

        return
