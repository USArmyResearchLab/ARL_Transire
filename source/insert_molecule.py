# copyright Caleb Michael Carlin (2018)
# Released uder Lesser Gnu Public License (LGPL)
# See LICENSE file for details.

from ase import Atoms
from ase.build import stack
from ase.io import write, read
from .utilities import InterfaceConfigStorage as ICS
from .utilities import printx, write_xyz
from ase.geometry import get_layers
import numpy as np
from .energy_calculators import EnergyCalc, calculate_binding_energy
import sys
from math import pi


class InsertMolecule(object):
    """
    Insert molecule or object in between surfaces of interface.

    Parameters:

    input: InputReader object
        object with read in keywords
    ics: InterfaceConfigStorage object
        object containing the interface to be manipulated
    surf_a: tuple
        tuple of Miller Indices for the top crystal. Used for generating
        filenames.
    surf_b: tuple
        tuple of Miller Indices for the bottom crystal. Used for generating
        filenames. 
    """

    def __init__(self, input, ics, surf_a, surf_b):
        self.input = input
        self.structure = ics
        self.surf_a = surf_a
        self.surf_b = surf_b
        if (input.dict['insert_file'] is None):
            printx("Error: a coordinate file of molecule to insert must"
                   " be specified.")
            sys.exit(1)

    def mol_insert(self):
        """
        Extract the atoms from one side of interface using tags,
        add the molecule to be inserted and re-stack the removed
        side on top of the inserted molecule.
        """
        temp = ICS()
        right_lead = Atoms()
        left_lead = Atoms()
        tags_list = self.structure.atom.get_tags()
        left_lead.cell[0][0] = self.structure.atom.cell[0][0]
        left_lead.cell[1][1] = self.structure.atom.cell[1][1]

        # extract the atoms from each layer of the interface
        for i in range(len(self.structure.atom)):
            # since pop deletes the atom in the interface,
            # we always take atom 0
            if (tags_list[i] == 2):
                right_lead.append(self.structure.atom.pop(0))
            else:
                left_lead.append(self.structure.atom.pop(0))

        lead_height = max(left_lead.positions[:, 2])

        file = self.input.dict['insert_file']
        insert = read(file)
        if self.input.dict['rotate_insert_x'] is not None:
            rot = float(self.input.dict['rotate_insert_x']) * pi / 180.0
            insert.rotate('x', rot, rotate_cell=True)
        if self.input.dict['rotate_insert_y'] is not None:
            rot = float(self.input.dict['rotate_insert_y']) * pi / 180.0
            insert.rotate('y', rot, rotate_cell=True)
        if self.input.dict['rotate_insert_z'] is not None:
            rot = float(self.input.dict['rotate_insert_z']) * pi / 180.0
            insert.rotate('z', rot, rotate_cell=True)
        self.shift_insert(insert)
#        for i in range(3):
#            insert.cell[i][i] = max(insert.positions[:, i])
        shift = (insert.cell[2][2] +
                 float(self.input.dict['insert_vacuum_below']) +
                 float(self.input.dict['insert_vacuum_above']))

        insert.positions[:, 2] += lead_height + \
            float(self.input.dict['insert_vacuum_below'])

        base = self.center_stack(left_lead, insert)
        for i in range(len(right_lead)):
            right_lead.positions[0, 2] += shift
            base.append(right_lead.pop(0))
        base.cell[2, 2] = self.structure.atom.cell[2, 2] + shift
        base.cell[0, 0] = self.structure.atom.cell[0, 0]
        base.cell[1, 1] = self.structure.atom.cell[1, 1]

        temp.atom = base.copy()
        temp.step = 'inserted'
        temp.file_name = self.input.dict['project_name'] + '.inserted'
        temp.energy = 1.0e10
        temp.unit_cell_a = self.structure.unit_cell_a
        temp.unit_cell_b = self.structure.unit_cell_b
        # write an xyz file to make it easier to view the result
        write_xyz(temp)
        if (self.input.dict['calc_insert_energy'] == 'True'):
            EnergyCalc(
                self.input, self.input.dict['energy_method'], temp, 'run')
            if self.input.dict['calculate_binding_energy'] != 'False':
                temp.energy = calculate_binding_energy(temp, self.input)
            printx("\n======Insertion Energy======\nEnergy = "
                    +str(temp.energy)+"\n",self.input.dict['output_file'])
        else:
            EnergyCalc(
                self.input, self.input.dict['energy_method'], temp, 'write')

        return temp

    def center_stack(self, atom_1, atom_2):
        """
        stacks atom_2 on atom_1 so that the center of the two cells are in line
        """
        center_1x = (atom_1.cell[0][0] + atom_1.cell[1][0]) / 2.0
        center_1y = (atom_1.cell[0][1] + atom_1.cell[1][1]) / 2.0
        center_2x = (atom_2.cell[0][0] + atom_2.cell[1][0]) / 2.0
        center_2y = (atom_2.cell[0][1] + atom_2.cell[1][1]) / 2.0

        del_x = center_2x - center_1x
        del_y = center_2y - center_1y
        atom_1.cell[2][2] += atom_2.cell[2][2]

        for i in range(len(atom_2)):
            atom_2.positions[0][0] -= del_x
            atom_2.positions[0][1] -= del_y
            atom_1.append(atom_2.pop(0))

        return atom_1

    def shift_insert(self, molecule):
        """
        check to see if we need to shift the inserted molecule entirely into
        the positive portion so that cp2k doesn't give errors.
        """
        x_min = min(molecule.positions[:, 0])
        y_min = min(molecule.positions[:, 1])
        z_min = min(molecule.positions[:, 2])

        if x_min < 0.0:
            molecule.positions[:, 0] -= x_min
        if y_min < 0.0:
            molecule.positions[:, 1] -= y_min
        if z_min < 0.0:
            molecule.positions[:, 2] -= z_min

        return
