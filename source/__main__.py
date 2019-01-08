# copyright Caleb Michael Carlin (2018)
# Released uder Lesser Gnu Public License (LGPL)
# See LICENSE file for details.

import ase
import re
import io
import sys
import numpy as np
import argparse
# -------------------
from ase import Atoms
from ase.io import read
from .interface import InterfaceSupercell
from .input_reader import ReadInput
from ase.transport.calculators import TransportCalculator as TCalc
from .constrained_search import ConstrainedSearch
from .electron_transport import ElectronTransport
from .utilities import InterfaceConfigStorage as ICS
from .utilities import write_xyz, printx
from .utilities import generate_surfaces_list
from .energy_calculators import EnergyCalc, calculate_binding_energy
from .angle_search import AngleSearch
from .insert_molecule import InsertMolecule
from .separation_optimizer import SeparationOpt


def main():
    """
    Main driver for Transire.  Divided into 4 parts.
    1. read input
    2. build initial interface structure
    3. perform searches
    4. perfrom ET calculation.
    """
# input reading
# =============
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, nargs=1, dest='input_file',
                        help="input file with general keywords")
    parser.add_argument('-o', type=str, nargs=1, dest='output_file',
                        help="output file for interface search summary.")
    parser.add_argument('-c', type=str, nargs=1, dest='cp2k_file',
                        help="file with cp2k specific input options")
    args = parser.parse_args()

    input = ReadInput()
    input.read_input(args.input_file[0])
    if (input.dict['cp2k_input'] is None):
        try:
            input.dict['cp2k_input'] = args.cp2k_file[0]
        except Exception as err:
            pass
    if (input.dict['output_file'] is None):
        try:
            input.dict['output_file'] = args.output_file[0]
        except Exception as err:
            pass

# building initial structure
# ==========================
    Crystal1 = read(input.dict['crys_a_file'])
    Crystal2 = read(input.dict['crys_b_file'])

    BB = InterfaceSupercell(Crystal1, Crystal2, input)
    BB.cut_surface(BB.surface_a, BB.surface_b)
    try:
        BB.generate_interface()
    except Exception as err:
        [printx(x, input.dict['output_file']) for x in err.args]
        if (input.dict['ignore_initial_structure'] != 'False'):
            input.dict['calculate_initial_energy'] = 'False'
        else:
            sys.exit(1)

    lowest_energy_config = ICS()
    lowest_energy_config.atom = BB.interface.copy()
    lowest_energy_config.file_name = input.dict['project_name'] + '.initial'
    lowest_energy_config.unit_cell_a = BB.cut_cell_a
    lowest_energy_config.unit_cell_b = BB.cut_cell_b
    lowest_energy_config.step = 'initial'

    if (input.dict['calculate_initial_energy'] != 'True'):
        write_xyz(lowest_energy_config)
        lowest_energy_config.energy = 1.0e8
    else:
        #    if (input.dict['calculate_initial_energy'] != 'False'):
        EnergyCalc(
            input, input.dict['energy_method'], lowest_energy_config, 'run')
        printx('initial interface energy is= ' +
               str(lowest_energy_config.energy))
        file = input.dict['output_file']
        printx('=========Initial Energy Completed=========', file)
        printx('Interface Energy = ' + str(lowest_energy_config.energy), file)
        if (input.dict['calculate_binding_energy'] != 'False'):
            binding_energy = calculate_binding_energy(
                lowest_energy_config, input)
            printx("Binding_energy = " + str(binding_energy), file)

# Searches
# ========
    # generate list of unique surfaces or copy over initial surface
    if (input.dict['surface_search'] != 'False'):
        surface_list_a, surface_list_b = generate_surfaces_list(input)
    else:
        surface_list_a = [input.dict['crys_a_surface']]
        surface_list_b = [input.dict['crys_b_surface']]
    # loop over surfaces and searches
    for x in surface_list_a:
        for y in surface_list_b:
            BB.cut_surface(x, y)
            lowest_energy_config.unit_cell_a = BB.cut_cell_a.copy()
            lowest_energy_config.unit_cell_b = BB.cut_cell_b.copy()
            if (input.dict['search_list'] is not None):
                for i in input.dict['search_list']:
                    if (i == 0.0):
                        MCS = ConstrainedSearch(
                            BB, input, lowest_energy_config, x, y)
                        lowest_energy_config = MCS.generate_step()
                    elif (i == 1.0):
                        ANS = AngleSearch(
                            input, BB, lowest_energy_config, x, y)
                        lowest_energy_config = ANS.angle_search()
                    elif (i == 2.0):
                        IM = InsertMolecule(
                            input, lowest_energy_config, x, y)
                        lowest_energy_config = IM.mol_insert()
                    elif (i == 3.0):
                        SeparationOpt(lowest_energy_config, input, x, y)
                    else:
                        printx('Invalid number given for search_list')
# ET calculation
# ==============
    if (input.dict['perform_ET'] == 'True'):
        #   make sure that the restart path isn't used if not restarting
        if(input.dict['ET_restart'] != 'True'):
            input.dict['restart_path'] = './'
            restart = False
        else:
            restart = True
        ET = ElectronTransport(input, lowest_energy_config, restart)
        ET.calculate_ET(lowest_energy_config.unit_cell_a,
                        lowest_energy_config.unit_cell_b)
