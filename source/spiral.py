# copyright Caleb Michael Carlin (2018)
# Released uder Lesser Gnu Public License (LGPL)
# See LICENSE file for details.

import numpy as np
try:
	from itertools import izip as zip
except ImportError: # will be 3.x series
	pass
from itertools import cycle, count
from operator import add
import sys


class Hyperspiral(object):
    """
    class for generating the positions of a 3D spiral so that we can
    move through a crystal lattice for the purpose of filling the
    orthorhombic representation of the crystal.  Each call of tick returns
    the next set of coordiantes on the grid, spiraling through the values of
    x and y bounded by max_x and max_y before repeating for different values
    of z up to max_z.
    """

    def __init__(self, max_values):
        self.max_x = max_values[0]
        self.max_y = max_values[1]
        self.max_z = max_values[2]
        self.z_range = range(-self.max_z, self.max_z + 1)
        self.z_spot = 0
        self.position = [0, 0, 0]
        self.movement = self.spiral_movements()

    def spiral_distances(self):
        """
        track the distance that is traveled each loop.
        """
        for distance in count(1):
            for _ in (0, 1):
                yield distance

    def clockwise_directions(self):
        """
        defines the movements that the step can take.
        """
        left = (-1, 0)
        right = (1, 0)
        up = (0, -1)
        down = (0, 1)
        return cycle((right, down, left, up))

    def spiral_movements(self):
        """
        factory for moving.
        """
        for distance, direction in zip(self.spiral_distances(),
                                        self.clockwise_directions()):
            for _ in range(distance):
                yield direction

    def tick(self):
        """
        manages advancing the spiral by one step.
        """
        # first check to see if we need to start a new layer or stop
        if self.check_end_of_layer():
            if len(self.z_range) != 1:
                # if we have completed a spiral in a plane, then we
                # need a new Z value and restarting the spiral
                self.movement = self.spiral_movements()
                if self.z_range[self.z_spot] != 0.0:
                    self.position = [0, 0, self.z_range[self.z_spot]]
                    self.z_spot += 1
                else:
                    # skip the z = 0 spot since we already did that
                    self.position = [0, 0, self.z_range[self.z_spot + 1]]
                    self.z_spot += 2
            else:
                raise Exception('Completed spiral without filling cell')
        # get the next displacement
        dx, dy = next(self.movement)
        self.position = list(map(add, self.position, [dx, dy, 0]))
        while self.out_of_bounds():
            dx, dy = next(self.movement)
            self.position = list(map(add, self.position, [dx, dy, 0]))
            if ((abs(self.position[0]) > self.max_x) and
               (abs(self.position[1]) > self.max_y)):
                break
        return self.position

    def out_of_bounds(self):
        """
        Test to see if the next step is outside our defined box.
        """
        if abs(self.position[0]) > self.max_x:
            return True
        if abs(self.position[1]) > self.max_y:
            return True
        if abs(self.position[2]) > self.max_z:
            return True
        return False

    def check_end_of_layer(self):
        """
        Checks to see if the spiral has completed all positions on a single
        plane.
        """
        if self.position[0] <= self.max_x:
            return False
        if self.position[1] <= self.max_y:
            return False
        return True
