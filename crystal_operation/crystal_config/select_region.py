__all__ = ['select_sphere_region']

import numpy as np
from .lattice import Lattice, Duplicate

def select_sphere_region(points, center, radius):
    """
    Select points that are within a sphere defined by center and radius.
    
    Parameters:
    points: np.array[n,3], atoms position x, y, z
    center: [x,y,z], the center of sphere
    radius: the radius of sphere
    
    Returns:
    region_index: np.array[n], binary array where 1 indicates point is inside sphere
    """
    # Calculate squared distances from each point to the center
    squared_distances = np.sum((points - center)**2, axis=1)
    
    # Compare with squared radius to avoid sqrt operation
    region_index = (squared_distances <= radius**2).astype(int)
    
    return region_index