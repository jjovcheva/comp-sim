'''
Uses periodic boundary conditions to:
1. Return the image of a point inside a cube of box size
set by mdutilities
2. Return the image of the point closest to the origin
'''

import mdutilities
import numpy as np

def pbc(vector, box_size):
    '''
    Returns the image of a point inside a cube of size
    set by mdutilities.

    :param vector: [3] float array
    :param box_size: box size set by mdutilities

    :return pos: position of the particle in the box
    '''
    pos = np.mod(vector, box_size)
    return pos

def mic(vector, box_size):
    '''
    Returns the image of the point closest to the origin.

    :param vector: [3] float array
    :param box_size: box size set by mdutilities

    :return: image of point closest to origin
    '''
    return np.subtract((np.mod(vector+box_size/2, box_size)),box_size/2)