from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier
import numpy as np
from scipy.spatial.transform import Rotation as R
import pandas as pd
# Extracts quaternion orientations from an OVITO-compatible dump file
def dump_ovito_quaternion(dump_name = "min-1200.data"):
    pipeline = import_file(dump_name)
    ptm = PolyhedralTemplateMatchingModifier()
    ptm.output_orientation = True
    pipeline.modifiers.append(ptm)
    data = pipeline.compute()
    orientations = data.particles['Orientation']
    return orientations

# Converts quaternions to [hkl] directions based on crystal orientation
def quaternions_to_hkls(quaternions):
    orientations = [
        np.array([0, 0, 0, 1]) if np.array_equal(x, np.array([0, 0, 0, 0])) else x
        for x in quaternions
    ]
    # Vectorized rotation calculation
    rots = R.from_quat(orientations)
    R_matrices = rots.as_matrix()
    hkls_crystal = np.linalg.inv(R_matrices) @ np.array([0, 0, 1])
    return hkls_crystal

# Normalize [hkl] directions to reduced integer form
def normalize_hkls(hkls):
    result = []
    for hkl in hkls:
        # Convert to NumPy array for processing
        hkl_arr = np.array(hkl, dtype=float)
        x, y, z = np.abs(hkl_arr)

        error_h, error_l = 100, 0.05
        for a, b, i, j in [(x, y, 0, 1), (x, z, 0, 2), (y, z, 1, 2)]:
            ratio = a / b if b != 0 else float('inf')
            if abs(ratio) > error_h or abs(ratio) < error_l:
                hkl_arr[j if a > b else i] = 0

        zero_count = np.count_nonzero(hkl_arr == 0)
        if zero_count == 2:
            # If two elements are zero, set the non-zero one to 1
            hkl_arr[hkl_arr != 0] = 1
        else:
            # Normalize to minimum non-zero component
            non_zero = hkl_arr[np.abs(hkl_arr) > 1e-6]
            if len(non_zero) > 0:
                norm_val = np.min(np.abs(non_zero))
                hkl_arr = hkl_arr / norm_val
                # Reduce to lowest integer ratio
                # int_hkl = np.round(hkl_arr).astype(int)
                # gcd_val = np.gcd.reduce(int_hkl)
                # if gcd_val != 0:
                #     hkl_arr = int_hkl / gcd_val

        result.append(hkl_arr.tolist())
    return result

# Stereographic projection of multiple crystal directions
def stereographic_projection(hkls, cartesian=True):
    """
    Calculate stereographic projection of multiple crystal planes [h, k, l] (using (001) as projection plane)
    
    Parameters:
        hkls: NumPy array with shape (n, 3), representing n sets of Miller indices
        cartesian: Boolean, if True returns Cartesian coordinates (x, y), if False returns polar coordinates (phi, r)
    
    Returns:
        If cartesian == True: 
            xy: NumPy array with shape (n, 2), projection coordinates (x, y) for each plane
        If cartesian == False:
            polar: NumPy array with shape (n, 2), polar coordinates (phi, r) for each plane
    """
    hkls = np.asarray(hkls)
    if hkls.ndim == 1:
        hkls = hkls.reshape(1, -1)
    
    h, k, l = hkls[:, 0], hkls[:, 1], hkls[:, 2]
    
    norms = np.linalg.norm(hkls, axis=1, ord=2)
    theta = np.arccos(l / norms)
    
    phi = np.arctan2(k, h)
    r = np.tan(theta / 2)
    
    if cartesian:
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        return np.column_stack((x, y))
    else:
        return np.column_stack((phi, r))

# Reduce orientations to the stereographic fundamental zone using cubic symmetry
def reduce_fundamental_zone(hkls):
    """
    Reduce Miller indices to the fundamental zone using cubic symmetry operations
    
    Parameters:
        hkls: NumPy array with shape (n, 3), representing n sets of Miller indices
    
    Returns:
        reduced_hkls: List of tuples, each containing the original hkl and its reduced form
    """
    symmetry_ops = R.create_group('O').as_matrix()
    reduced_hkls = []
    
    # Reference points defining the fundamental zone boundary
    ref_points = np.array([
        [0.0, 0.0],
        [0.3660, 0.3660],
        [0.4142, 0.0000]
    ])
    
    for hkl in hkls:
        # Generate all symmetrically equivalent directions
        all_equiv_rots = np.array([s @ hkl for s in symmetry_ops])
        
        # Project all directions to the stereographic plane
        all_projection = stereographic_projection(all_equiv_rots)
        
        # Filter valid projections within the fundamental zone
        for i, (x, y) in enumerate(all_projection):
            if x >= 0 and y >= 0:
                distances = np.sqrt(np.sum((np.array([x, y]) - ref_points)**2, axis=1))
                print(distances)
                if sum(distances) < 0.97:
                    print(sum(distances))
                    if np.all(distances <= 0.518):
                        reduced_hkls.append(all_equiv_rots[i])
                        break
    return reduced_hkls

# Example usage
if __name__ == "__main__":
    orientations = dump_ovito_quaternion("test.dump")
    print(len(orientations))
    # orientations = orientations[20000:20005]
    # hkls_crystal = quaternions_to_hkls(orientations)
    # print(hkls_crystal)
    # hkls_norm = normalize_hkls(hkls_crystal)
    # print(hkls_norm)
    # hkls_reduce = reduce_fundamental_zone(hkls_norm)
    # print(hkls_reduce)
