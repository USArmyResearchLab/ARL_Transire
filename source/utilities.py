# copyright Caleb Michael Carlin (2018)
# Released uder Lesser Gnu Public License (LGPL)
# See LICENSE file for details.

import io
import numpy as np
import fractions
from ase import Atoms
import itertools
import sys
from numpy.linalg import norm, solve
from ase.utils import gcd
from ase.build import bulk


def lcm(a, b):
    """
    determines the least common multiple of two integers (a and b)
    """
    if (a and b):
        j = abs(a * b) / fractions.gcd(a, b)
        return j
    else:
        return 0


def write_xyz(atom):
    """
    prints out the cell dimensions and coordinates of an atom object
    in extended xyz format.
    """
    filename = atom.file_name + '.xyz'

    printx(str(len(atom.atom)),filename)
    celldim = ''
    for i in range(3):
        for j in range(3):
            celldim += str(atom.atom.cell[i][j]) + ' '
    printx('Lattice = \"' + celldim + '\"',filename)

    symbols = atom.atom.get_chemical_symbols()
    coords = atom.atom.get_positions()
    for k in range(len(atom.atom)):
        printx(symbols[k] + ' ' + str(coords[k][0]) + ' ' +
                     str(coords[k][1]) + ' ' +
                     str(coords[k][2]),filename)



class InterfaceConfigStorage(object):
    """
    object to hold an atom object, the energy_per_area, output name,
    step number if applicable, and the unstacked left supercell
    (for insertion purposes).

    atom: ASE Atoms object
        The interface structure
    step: integer
        The step number or value used for printing unique filenames
    energy: float
        Either the total energy or binding energy per unit area of interface.
        Units are determined by method used to calculate.
    file_name: string
        File name used for printing output and testing for duplicate
        energy calculations.
    unit_cell_a: ASE Atoms object
        Top unit cell structure used to construct interface structure
    unit_cell_b: ASE Atoms object
        Bottom unit cell structure used to construct interface structure
    ortho_cell_a: ASE Atoms object
        The orthorhombic representation of unit_cell_a
    ortho_cell_b: ASE Atoms object
        The orthorhombic representation of unit_cell_b
    """

    def __init__(self):
        self.atom = None
        self.step = None
        self.energy = None
        self.file_name = None
        self.unit_cell_a = None
        self.unit_cell_b = None

    def copy(self, no_filename=False):
        """
        creates a new copy of each of the objects in the ICS and returns
        it as a new ICS.  If no_filename is set to True, then the
        filename object is skipped.
        """
        temp_ics = InterfaceConfigStorage()
        temp_ics.atom = self.atom.copy()
        if (self.step is not None):
            temp_ics.step = self.step
        if (self.energy is not None):
            temp_ics.energy = self.energy
        if (self.file_name != 0 and not no_filename):
            temp_ics.file_name = self.file_name
        if (self.unit_cell_a is not None):
            temp_ics.unit_cell_a = self.unit_cell_a.copy()
        if (self.unit_cell_b is not None):
            temp_ics.unit_cell_b = self.unit_cell_b.copy()

        return temp_ics

    def tag_copy(self, tag):
        """
        copies the ics but returns only the atoms that match tag.
        """

        temp_ics = InterfaceConfigStorage()
        temp_atoms = [x for x in self.atom if (x.tag == tag)]
        temp_ics.atom = Atoms()
        temp_ics.atom.pbc = self.atom.pbc
        [temp_ics.atom.extend(x) for x in temp_atoms]
        temp_ics.atom.cell = self.atom.cell
        if (self.file_name != 0):
            temp_ics.file_name = self.file_name + '_' + str(tag)

        return temp_ics

    def len(self, tag):
        """
        returns the number of atoms with 'tag'
        """

        tag_list = self.atom.get_tags()
        side_tags = []
        [side_tags.append(x) for x in tag_list if x == tag]

        return len(side_tags)


def generate_surfaces_list(input):
    """
    Read in the user generated list of possible values for the miller planes
    for each surface.  Generate a list of tuples with all possible
    combinations. remove duplicates and any combinations involving (0,0,0)
    """
    x_range_a = input.dict['range_surface_h_a']
    y_range_a = input.dict['range_surface_k_a']
    z_range_a = input.dict['range_surface_l_a']
    x_range_b = input.dict['range_surface_h_b']
    y_range_b = input.dict['range_surface_k_b']
    z_range_b = input.dict['range_surface_l_b']
    surface_list_a = []
    surface_list_b = []

    for i, j, k in itertools.product(x_range_a, y_range_a, z_range_a):
        sum = abs(i) + abs(j) + abs(k)
        if (sum != 0):
            surface_list_a.append((i, j, k))

    for i, j, k in itertools.product(x_range_b, y_range_b, z_range_b):
        sum = abs(i) + abs(j) + abs(k)
        if (sum != 0):
            surface_list_b.append((i, j, k))

    surface_list_a = set(surface_list_a)
    surface_list_b = set(surface_list_b)

    return surface_list_a, surface_list_b

def printx(text, file=None):
    """
    Function for handling printing to standard out and files.
    If no file is given, then output is printed to standard out.
    If mpi4py is being used, then only the head node performs the printing.
    """

    parallel = True
    try:
        from mpi4py import MPI
    except Exception as err:
        parallel = False

    if (parallel):
        if (file is None):
            if (MPI.COMM_WORLD.Get_rank() == 0):
                print(text)
        else:
            if (MPI.COMM_WORLD.Get_rank() == 0):
                with open(file, 'a') as outfile:
                    outfile.write(text + '\n')
    else:
        if (file is None):
            print(text)
        else:
            with open(file, 'a') as outfile:
                outfile.write(text + '\n')


def angle_between(v1, v2):
    """
    Function for finding the angle in radians between two vectors.
    """
    unit_v1 = unit_vector(v1)
    unit_v2 = unit_vector(v2)
    return np.arccos(np.clip(np.dot(unit_v1, unit_v2), -1.0, 1.0))


def unit_vector(vector):
    """
    function for finding the unit vector of an np.array.
    """
    return vector / np.linalg.norm(vector)


def surface_from_ase(lattice, indices, layers, vacuum=None, tol=1e-10):
    """Create surface from a given lattice and Miller indices.

    lattice: Atoms object or str
        Bulk lattice structure of alloy or pure metal.  Note that the
        unit-cell must be the conventional cell - not the primitive cell.
        One can also give the chemical symbol as a string, in which case the
        correct bulk lattice will be generated automatically.
    indices: sequence of three int
        Surface normal in Miller indices (h,k,l).
    layers: int
        Number of equivalent layers of the slab.
    vacuum: float
        Amount of vacuum added on both sides of the slab.
    
    Originally from ASE v. 3.11
    """

    indices = np.asarray(indices)

    if indices.shape != (3,) or not indices.any() or indices.dtype != int:
        raise ValueError('%s is an invalid surface type' % indices)

    if isinstance(lattice, str):
        lattice = bulk(lattice, cubic=True)

    h, k, l = indices
    h0, k0, l0 = (indices == 0)
    if h0 and k0 or h0 and l0 or k0 and l0:  # if two indices are zero
        if not h0:
            c1, c2, c3 = [(0, 1, 0), (0, 0, 1), (1, 0, 0)]
        if not k0:
            c1, c2, c3 = [(0, 0, 1), (1, 0, 0), (0, 1, 0)]
        if not l0:
            c1, c2, c3 = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    else:
        p, q = ext_gcd(k, l)
        a1, a2, a3 = lattice.cell

        # constants describing the dot product of basis c1 and c2:
        # dot(c1,c2) = k1+i*k2, i in Z
        k1 = np.dot(p * (k * a1 - h * a2) + q * (l * a1 - h * a3),
                    l * a2 - k * a3)
        k2 = np.dot(l * (k * a1 - h * a2) - k * (l * a1 - h * a3),
                    l * a2 - k * a3)

        if abs(k2) > tol:
            i = -int(round(k1 / k2))  # i corresponding to the optimal basis
            p, q = p + i * l, q - i * k

        a, b = ext_gcd(p * k + q * l, h)

        c1 = (p * k + q * l, -p * h, -q * h)
        c2 = np.array((0, l, -k)) // abs(gcd(l, k))
        c3 = (b, a * p, a * q)

    surf = build_from_ase(lattice, np.array([c1, c2, c3]), layers, tol)
    if vacuum is not None:
        surf.center(vacuum=vacuum, axis=2)
    return surf


def build_from_ase(lattice, basis, layers, tol):
    """
    older version of build from ase to counteract the update that broke
    how arl_transire works.

    Originally from ASE v. 3.11
    """
    surf = lattice.copy()
    scaled = solve(basis.T, surf.get_scaled_positions().T).T
    scaled -= np.floor(scaled + tol)
    surf.set_scaled_positions(scaled)
    surf.set_cell(np.dot(basis, surf.cell), scale_atoms=True)
    surf *= (1, 1, layers)

    a1, a2, a3 = surf.cell
    surf.set_cell([a1, a2,
                   np.cross(a1, a2) * np.dot(a3, np.cross(a1, a2)) /
                   norm(np.cross(a1, a2))**2])

    # Change unit cell to have the x-axis parallel with a surface vector
    # and z perpendicular to the surface:
    a1, a2, a3 = surf.cell
    surf.set_cell([(norm(a1), 0, 0),
                   (np.dot(a1, a2) / norm(a1),
                    np.sqrt(norm(a2)**2 - (np.dot(a1, a2) / norm(a1))**2), 0),
                   (0, 0, norm(a3))],
                  scale_atoms=True)

    surf.pbc = (True, True, False)

    # Move atoms into the unit cell:
    scaled = surf.get_scaled_positions()
    scaled[:, :2] %= 1
    surf.set_scaled_positions(scaled)

    return surf


def ext_gcd(a, b):
    """
    calculate the greatest common denominator or return integers if
    b is 0 or a is an integer multiple of b.
    """
    if b == 0:
        return 1, 0
    elif a % b == 0:
        return 0, 1
    else:
        x, y = ext_gcd(b, a % b)
        return y, x - y * (a // b)


def almost_zero(number):
    """
    quick function to tell if a number should be zero but is off
    due to computer error
    """

    return (abs(number) < 1e-10)
