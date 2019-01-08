# copyright Caleb Michael Carlin (2018)
# Released uder Lesser Gnu Public License (LGPL)
# See LICENSE file for details.

import ase
from ase import Atoms, Atom
import numpy as np
from numpy.linalg import norm
import itertools
import fractions
from math import pi, floor
from ase.build import cut, make_supercell
from ase.build import stack, surface
from .utilities import lcm, printx, angle_between, surface_from_ase
from .utilities import almost_zero
from ase.io import write as ase_write
#import sympy
from traceback import print_exc
import re
import sys
from .spiral import Hyperspiral


class InterfaceSupercell(object):
    """
    Class for generating interface structures from two unit cells.
   
    Parameters:

    unit_cell_a: ASE atoms object
        unit cell atoms object for the top side of the interface
    unit_cell_b: ASE atoms object
        unit cell atoms object for the bottom side of the interface
    input: InputReader object
        object with read in keywords
    """

    def __init__(self, unit_cell_a, unit_cell_b, input):
        self.raw_cell_a = unit_cell_a
        self.raw_cell_b = unit_cell_b
        self.cut_cell_a = unit_cell_a
        self.cut_cell_b = unit_cell_b
        self.layers_a = int(input.dict['crys_a_layers'])
        self.layers_b = int(input.dict['crys_b_layers'])
        self.surface_a = input.dict['crys_a_surface']
        self.surface_b = input.dict['crys_b_surface']
        self.interface = None
        self.distance = float(input.dict['separation'])
        self.super_cell1 = None
        self.super_cell2 = None
        self.input = input
        self.duplicates = []
        self.ase_version = self.determine_version()

    def cut_surface(self, surface_a, surface_b):
        """
        cut the raw unit cell to produce the cut unit cell
        """
        self.cut_cell_a = surface_from_ase(
            self.raw_cell_a, surface_a, layers=1)
        self.cut_cell_b = surface_from_ase(
            self.raw_cell_b, surface_b, layers=1)

    def generate_interface(self):
        """
        main function to generate interface from two cut_cells.
        Returns .True. if an error occured, otherwise returns .False.
        """

        # set up the output file so we can put the error messages
        # in the output file
        file = self.input.dict['output_file']

        # copy over the atom units so that cut_cells are unchanged
        unit_cell_a = self.cut_cell_a.copy()
        unit_cell_b = self.cut_cell_b.copy()

        # remove small numbers from computer error
#        unit_cell_a.cell = self.check_zero_diag(unit_cell_a.cell)
#        unit_cell_a = self.check_cell_rotation(unit_cell_a, unit_cell_b)

        # =====Debug=====
        if (self.input.dict['print_debug'] != 'False'):
            printx("========Starting Cell 1========")
            printx(str(unit_cell_a.cell))
            printx("atoms = " + str(len(unit_cell_a)))
            printx("========Starting Cell 2========")
            printx(str(unit_cell_b.cell))
            printx("atoms = " + str(len(unit_cell_b)))

        # replicate the unit cells so that periodicity is always preserved
        periodic_cell_a, max_coeff_a = self.protect_periodicity(unit_cell_a)
        periodic_cell_b, max_coeff_b = self.protect_periodicity(unit_cell_b)

        # populate the new cell using cookie cutter method on generated lattice
        try:
            self.super_cell1 = self.populate_new_cell(
                unit_cell_a, periodic_cell_a, max_coeff_a)
            self.super_cell2 = self.populate_new_cell(
                unit_cell_b, periodic_cell_b, max_coeff_b)
        except Exception as err:
            [printx(x) for x in err.args]
            raise Exception("Too many atoms, skipping to next step")

        # =====Debug=====
        if (self.input.dict['print_debug'] != 'False'):
            printx("========Ortho Cell 1========")
            printx(str(self.super_cell1.cell))
            printx("atoms = " + str(len(self.super_cell1)))
            printx("========Ortho Cell 2========")
            printx(str(self.super_cell2.cell))
            printx("atoms = " + str(len(self.super_cell2)))

        # calculate the smallest supercells needed to minimize
        # stress in the interface
        P_list, R_list = self.generate_interface_transform(
            self.super_cell1, self.super_cell2)
        P_tuple = tuple(P_list + [int(self.input.dict['crys_a_layers'])])
        R_tuple = tuple(R_list + [int(self.input.dict['crys_b_layers'])])
        # generate new supercells
        try:
            self.super_cell1 *= P_tuple
        except Exception as err:
            raise Exception(
                "Error in generating supercell_a in interface step")
        try:
            self.super_cell2 *= R_tuple
        except Exception as err:
            raise Exception(
                "Error in generating supercell_b in interface step")

        # =====Debug=====
        if (self.input.dict['print_debug'] != 'False'):
            printx("Replication A = " + str(P_tuple))
            printx("Replication B = " + str(R_tuple))

        # check that total size isn't too big before we continue
        total = len(self.super_cell1) + len(self.super_cell2)
        if (total >= int(self.input.dict['max_atoms'])):
            raise Exception("Error: interface is too large: " + str(total))

        # tag the two supercells so that they can be separated later
        self.super_cell1.set_tags(1)
        self.super_cell2.set_tags(2)


        # add a vacuum between the layers.
        if (self.distance is not None):
            self.super_cell1.cell[2, 2] += self.distance

        # =====Debug=====
        if (self.input.dict['print_debug'] != 'False'):
            printx("========Super Cell 1========")
            printx(str(self.super_cell1.cell))
            printx("atoms = " + str(len(self.super_cell1)))
            printx("========Super Cell 2========")
            printx(str(self.super_cell2.cell))
            printx("atoms = " + str(len(self.super_cell2)))

        # stack the supercells on top of each other and set pbc to xy-slab
        try:
            self.interface, self.super_cell1, self.supercell2 = stack(
                self.super_cell1, self.super_cell2,
                output_strained=True, maxstrain=None)
        except Exception as err:
            raise Exception(
                "Error in generating interface during the stack step")

        # set pbc to infinite slab or fully periodic setting
        if (self.input.dict['full_periodicity'] != 'False'):
            self.interface.pbc = [1, 1, 1]
        else:
            self.interface.pbc = [1, 1, 0]

        #add explicit vacuum above and below 
        if (self.input.dict['z_axis_vacuum'] != '0.0'):
            self.interface = self.z_insert_vacuum(self.interface)

        # use merge sort to identify and remove duplicate atoms.
        if (self.input.dict['remove_duplicates'] != 'False'):
            try:
                self.remove_duplicates(self.interface)
            except Exception as err:
                raise Exception("Error in checking for atom overlaps")

        return

    def match_sizes(self, cell_side_a, cell_side_b):
        """
        the unit cells must be replicated an integer number of times.
        Using a back and forth method, the two integers are determined that
        reduces the difference between the two values is less
        than the tolerance given in the input.
        """
        a = 1.0
        b = 1.0
        convergence = cell_side_a / cell_side_b
        upper_bound = 1.0 + float(self.input.dict['tolerance'])
        lower_bound = 1.0 - float(self.input.dict['tolerance'])
        while ((convergence < lower_bound) or (convergence > upper_bound)):
            if (cell_side_a * a) < (cell_side_b * b):
                a += 1.0
                convergence = (cell_side_a * a) / (cell_side_b * b)
            else:
                b += 1.0
                convergence = (cell_side_a * a) / (cell_side_b * b)

        return a, b

    def generate_interface_transform(self, unit_cell_a, unit_cell_b):
        """
        A pair of lists for replicating the two cells so that they match
        up are generated.
        """
        P_list = []
        R_list = []

        for j in range(2):
            side_a = unit_cell_a.cell[j][j]
            side_b = unit_cell_b.cell[j][j]
            x, y = self.match_sizes(side_a, side_b)

            P_list.append(abs(int(x)))
            R_list.append(abs(int(y)))

        return P_list, R_list

    def protect_periodicity(self, unit_cell):
        """
        determines the number of copies of unit cell along each axis are
        needed to ensure that any further replication of the supercell will
        accurately match up the cells.
        """
        new_cell = unit_cell.cell.copy()

        max_coeff = [0, 0, 0]
        for i, j in [[0, 1], [0, 2], [1, 2]]:
            new_cell, max_coeff = self.twod_matching(i, j, new_cell, max_coeff)

        # if the cell has been rotated so that the vector lies on the negative
        # axis, then we need to for it to have a non-zero max_coeff
        #max_coeff = [x+1 for x in max_coeff if x == 0 else x]
        max_coeff = [x+1 if x == 0 else x for x in max_coeff]

        return new_cell, max_coeff

    def twod_matching(self, axis1, axis2, matrix, max_coeff):
        """
        Take a two by two matrix that is the projection of two lattice
        vectors on a plane and determine the retangular representation of
        the matrix.
        """

        a_b = np.array([matrix[axis1, axis1], matrix[axis1, axis2]])
        b_a = np.array([matrix[axis2, axis1], matrix[axis2, axis2]])
        coeff_storage = [[0, 0], [0, 0]]

        intersect = [0, 0]
        # if the two numbers we are looking at are not zero,
        # we find the integer multiples needed to get a zero intercept
        if not(almost_zero(a_b[1]) and almost_zero(b_a[0])):
            for i in range(2):
                if not almost_zero(a_b[i] * b_a[i]):
                    c, d = self.match_sizes(abs(a_b[i]), abs(b_a[i]))
                    if np.sign(a_b[i]) == np.sign(b_a[i]):
                        c *= -1
                else:
                    if almost_zero(a_b[i]):
                        c, d = 1, 0
                    else:
                        c, d = 0, 1
                coeff_storage[i] = [abs(int(c)), abs(int(d))]
                b = (i + 1) % 2
                intersect[b] = c * a_b[b] + d * b_a[b]
            # store the values of the zero intercept
            matrix[axis1, axis1] = abs(intersect[0])
            matrix[axis2, axis2] = abs(intersect[1])
            matrix[axis1, axis2], matrix[axis2, axis1] = 0.0, 0.0
            # store the max coefficient for later when we populate the
            # the new orthonormal cell
            max_1_coeff = coeff_storage[0][0] + coeff_storage[1][0]
            max_2_coeff = coeff_storage[0][1] + coeff_storage[1][1]
            max_coeff[axis1] = max(abs(max_1_coeff), max_coeff[axis1])
            max_coeff[axis2] = max(abs(max_2_coeff), max_coeff[axis2])

        return abs(matrix), max_coeff

    def check_zero_diag(self, cell_matrix):
        """
        check if any diagonal elements are zero.  If there are any, swap
        the rows around until the check is passed.
        """
        while True:
            for i in range(3):
                if almost_zero(cell_matrix[i, i]):
                    next_val = (i + 1) % 3
                    cell_matrix = self.swap_rows(cell_matrix, i, next_val)
                    break
            # only way we get here is if all the diagonals are non-zero
            return cell_matrix

    def swap_rows(self, cell, row1, row2):
        """
        swap two rows in an array.
        """
        cell[row1, :], cell[row2, :] = cell[row2, :].copy(), cell[row1,
                                                                  :].copy()

        return cell

    def translate_cell(self, cut_cell, translation):
        """
        translate the atoms in the cell and then wrap back into the cell.
        """
        cut_cell.translate(translation)
        cut_cell.wrap(pbc=(1, 1, 0))

        return cut_cell

    def rotate_cell(self, cut_cell, rotation):
        """
        rotate the atoms and the cell vectors.
        """
        cut_cell.rotate('z', rotation, rotate_cell=True)

        return cut_cell

#    def check_cell_rotation(self, atom_a, atom_b):
#        """
#        Translate the top unit cell so that each corner of the unit_cell
#        is moved to the origin and then the area of overlap with the lower
#        unit_cell.  Returns the unit cell with the largest overlap area.
#        """
#        large_area = atom_a.copy()
#        largest_area = 0.0
#        Y_one = (atom_b.cell[0, 0], atom_b.cell[0, 1])
#        Y_two = (atom_b.cell[1, 0], atom_b.cell[1, 1])
#        # -1 corresponds with an inversion around X and Y axes
#        chg = [(1, 1), (-1, 1), (1, -1), (-1, 1)]
#        position = (1, 1)
#
#        for value in (chg):
#            # modify the cell dimensions and wrap atoms back into cell
#            atom_a.cell[0, 0] *= value[0]
#            atom_a.cell[0, 1] *= value[0]
#            atom_a.cell[1, 0] *= value[1]
#            atom_a.cell[1, 1] *= value[1]
#            atom_a.wrap(pbc=(1, 1, 0))
#            X_one = (atom_a.cell[0, 0], atom_a.cell[0, 1])
#            X_two = (atom_a.cell[1, 0], atom_a.cell[1, 1])
#            area = self.convex_intersect(X_one, X_two, Y_one, Y_two)
#            if (area > largest_area):
#                large_area = atom_a.copy()
#                largest_area = area
#                position = value
#            # to ensure right-handedness is maintained,
#            # swap first and second row if only one switch was used
#        if (position == (-1, 1) or position == (1, -1)):
#            self.swap_rows(large_area.cell, 0, 1)
#
#        return large_area
#
#    def convex_intersect(self, X_one, X_two, Y_one, Y_two):
#        """
#        use sympy to calculate the area of the overlap of the
#        parallelograms in the XY-plane of the two unit cells.
#        """
#        X_far = (X_one[0] + X_two[0], X_one[1] + X_two[1])
#        Y_far = (Y_one[0] + Y_two[0], Y_one[1] + Y_two[1])
#        Origin = (0.0, 0.0)
#        vertices = []
#        # Creates the polygons in sympy.  The vertices must be given in
#        # clockwise or counter-clockwise order.  We do some gymnastics
#        # to convert our ASE atom objects into sympy objects.
#        X_vert = [Origin, X_one, X_far, X_two]
#        Y_vert = [Origin, Y_one, Y_far, Y_two]
#        poly_a = Polygon(*X_vert)
#        poly_b = Polygon(*Y_vert)
#        # returns the points and/or line segments of intersection
#        intersect = poly_a.intersect(poly_b)
#        # But because the result isn't always an iterable variable, we
#        # have to use regex to read out all the points or end points of
#        # the line segments
#        rex = re.compile('Point2D\([0-9/0-9 ]*,[0-9/0-9 ]*\)')
#        found = re.findall(rex, str(intersect))
#        for i in found:
#            temp = eval(i)
#            vertices.append(temp)
#        # Loop over vertices of polygons to see if any are enclosed in the
#        # other polygon.
#        for j in poly_a.vertices:
#            if poly_b.encloses_point(j):
#                vertices.append(j)
#        for k in poly_b.vertices:
#            if poly_a.encloses_point(k):
#                vertices.append(j)
#        # Remove duplicate points in our list of overlap vertices
#        vertices = list(set(vertices))
#        # Sort the vertices of the overlap vertices so that
#        # they are counter-clockwise.
#        # Create the overlap polygon and return the area
#        vertices = self.sort_points(vertices)
#        P_over = Polygon(*vertices)
#        try:
#            area = P_over.area
#        except Exception as err:
#            area = 0.0
#
#        return area

    def sort_points(self, vertices):
        """
        Sort a set of points so that they are in counter-clockwise order,
        starting at (0,0).
        """
        ordered = []
        dict_store = {'0': Point2D(0.0, 0.0)}
        # calculate the 2D determinant between each point (and origin)
        # and count up how many are positive.
        # exclude any i=j and if i or j is the origin.
        for i in range(len(vertices)):
            counter = 1
            if (vertices[i] == Point2D(0.0, 0.0)):
                continue
            for j in range(len(vertices)):
                if (vertices[j] == Point2D(0.0, 0.0)):
                    continue
                if (i == j):
                    continue
                det = (vertices[i][0] * vertices[j][1] -
                       vertices[i][1] * vertices[j][0])
                # if two points are colinear,
                # set the count so that ordering remains counter-clockwise
                if (det == 0.0):
                    if (vertices[i][0] == vertices[j][0]):
                        if (vertices[i][1] < vertices[j][1]):
                            counter += 1
                            continue
                    elif (vertices[i][1] == vertices[j][1]):
                        if (vertices[i][0] > vertices[j][0]):
                            counter += 1
                            continue
                    else:
                        a_val = abs(vertices[i][0] + vertices[i][1])
                        b_val = abs(vertices[j][0] + vertices[j][1])
                        if (a_val > b_val):
                            counter += 1
                            continue
                elif (det < 0):
                    counter += 1
        # Create a dictionary where the number of positive determinants
        # is the index and the point is the value
            dict_store[str(counter)] = vertices[i]
        for j in range(len(vertices)):
            ordered.append(dict_store[str(j)])

        return ordered

    def remove_duplicates(self, atom):
        """
        remove duplicates to 0.01 accuracy by turning the coordinates into
        a single string and then removing duplicates
        """
        positions = atom.get_positions()
        reduced_coord = []
        for i in positions:
            coord = (str(i[0].round(2)) +
                     str(i[1].round(2)) + str(i[2].round(2)))
            reduced_coord.append(coord)
        dupes = self.find_duplicates(reduced_coord)
        del atom[dupes]

    def find_duplicates(self, seq, idfun=None):
        """
        method for finding duplicates efficiently.
        """
        dupes = []
        if idfun is None:
            def idfun(x): return x
        seen = {}
        result = []
        for item in range(len(seq)):
            marker = idfun(seq[item])
            if marker in seen:
                dupes.append(item)
                continue
            seen[marker] = 1
            result.append(seq[item])
        return dupes

    def determine_version(self):
        """
        determine the version of ase being used since the update after 3.12
        the way of building structures changed.
        """

        asever = (ase.__version__).split('.')

        return int(asever[1])

    def in_new_cell(self, atom, cell, error):
        """
        quick function to see an atom is inside a cell with the given error.
        """
        if (atom[0] < 0.0 - error) or (atom[0] > (cell[0, 0] - error)):
            return False
        if (atom[1] < 0.0 - error) or (atom[1] > (cell[1, 1] - error)):
            return False
        if (atom[2] < 0.0 - error) or (atom[2] > (cell[2, 2] - error)):
            return False
        return True

    def populate_new_cell(self, unit_cell, new_cell, max_coeff):
        """
        Fill up an orthorhombic cell wiht the atoms from a unit cell.
        Each atom is translated by a multiple of the old lattice vectors,
        and accepted atoms are added to the new object until the atom density
        matches that of the unit cell.
        """

        super_cell = Atoms()
        super_cell.set_cell(new_cell)
        # setup storage for rejected atoms in case we need them
        rejects = Atoms()
        volume = unit_cell.get_volume()
        new_volume = super_cell.get_volume()
        atoms = int(round(float(len(unit_cell)) * new_volume / volume))
        # quick check to see if the new cell will have too many atoms
        if (atoms > int(self.input.dict['max_atoms'])):
            raise Exception("too many atoms in supercell")
        vectors = np.asarray(unit_cell.cell)
        spiral = Hyperspiral(max_coeff)
        atom_positions = unit_cell.get_positions()
        # have to zero out infinitesimal values in atom_positions

        # =====Debug=====
        if self.input.dict['print_debug'] != "False":
            printx("old cell = " + str(unit_cell.cell))
            printx("new cell = " + str(new_cell))
            printx("max_coeff = " + str(max_coeff))
            printx("atoms = " + str(atoms))

        # move through the representations of the initial unit cell along a
        # spiral pattern on a grid of integer values.  first the spiral on
        # a plane is completed and then the spiral is shifted down and then up
        # along the third coordinate.
        while True:
            shift = np.matmul(spiral.position, vectors)
            for i in range(len(unit_cell)):
                atom_prime = np.add(shift, atom_positions[i])
                if self.in_new_cell(atom_prime, new_cell, 1e-7):
                    new_atom = unit_cell[i]
                    new_atom.position = atom_prime
                    super_cell.append(new_atom)
                    atoms -= 1
                    # satisfying this condition means success
                    if atoms == 0:
                        return super_cell
                else:
                    new_atom = unit_cell[i]
                    new_atom.position = atom_prime
                    rejects.append(new_atom)

            # if we get to the end of the spirals then we check
            # the edges for barely rejected atoms to add in
            try:
                spiral.tick()
            except Exception as err:
                [printx(x) for x in err.args]
                if self.input.dict['print_debug'] != 'False':
                    print_exc()
                try:
                    super_cell = self.check_edges(
                        rejects, new_cell, super_cell, atoms)
                except Exception as err:
                    raise Exception(err.args[0])
                return super_cell

        return super_cell

    def check_edges(self, rejects, new_cell, super_cell, atoms):
        """
        go through the rejected atoms to find one that is close enough to our
        boundries that we can add it in for edge cases.
        """

        for i in range(len(rejects)):
            if self.in_new_cell(rejects[i].position, new_cell, 1e-3):
                super_cell.append(rejects[i])
                atoms -= 1
                if atoms == 0:
                    return super_cell

        # if we get here, then we have failed to make the super_cell
        raise Exception("Error: failed to populate the cell")

        return super_cell

    def z_insert_vacuum(self, interface):
        """
        Add a vacuum above and below the crystal slabs by editing cell
        and shifting atoms.
        """

        vacuum = float(self.input.dict['z_axis_vacuum'])/2.0

        interface.cell[2][2] += vacuum*2

        for x in range(len(interface)):
            interface.positions[x,2] += vacuum

        return interface
