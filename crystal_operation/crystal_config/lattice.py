import numpy as np

##### 基础函数 ######
# 定义晶体结构
class Lattice:
    def __init__(self, name, point, axis, size=None, point_types=None, point_properties=None):
        self.name = name
        self.point = np.array(point)
        self.axis = np.array(axis)
        self.size = np.array(size) if size is not None else np.array([1, 1, 1])
        if point_types is None:
            self.point_types = np.zeros(len(point), dtype=int)
        else:
            self.point_types = np.array(point_types)
        self.point_properties = point_properties
        
def setting_type(Lattice, index, value=None):
    """
    Set point types in a Lattice object based on the given index array.
    
    Parameters:
    Lattice: Lattice object containing point_types to modify
    index: np.array[n] - binary array indicating which points to modify
    value: int - value to set for selected points (default: None)
    
    Returns:
    Lattice: Modified Lattice object
    
    Raises:
    ValueError: If length of index array doesn't match number of points
    """
    if len(Lattice.point_types) != len(index):
        raise ValueError("The number of atoms does not match atoms type array length")
    
    if value is None:
        value = 1
    # Set point_types where index is non-zero to the specified value  
    Lattice.point_types[index != 0] = value
    
    return Lattice
        
# 立方
def SC():  # Simple Cubic
    return Lattice('SC', [[0, 0, 0]], np.eye(3))

def BCC():  # Body-Centered Cubic
    return Lattice('BCC', [[0, 0, 0], [0.5, 0.5, 0.5]], np.eye(3))

def FCC():  # Face-Centered Cubic
    return Lattice('FCC', [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]], np.eye(3))

# 六方 / 三方
def HCP(c_a_ratio=np.sqrt(8/3)):  # Hexagonal Close-Packed
    axis = [[1, 0, 0], 
            [-0.5, np.sqrt(3)/2, 0], 
            [0, 0, c_a_ratio]]
    return Lattice('HCP', [[0, 0, 0], [2/3, 1/3, 0.5]], axis)

# 四方
def ST(z = 1.5):  # Simple Tetragonal
    axis = [[1, 0, 0], [0, 1, 0], [0, 0, z]]
    return Lattice('ST', [[0, 0, 0]], axis)

def BCT(z = 1.5):  # Body-Centered Tetragonal
    axis = [[1, 0, 0], [0, 1, 0], [0, 0, z]]
    return Lattice('BCT', [[0, 0, 0], [0.5, 0.5, 0.5]], axis)

# 正交
def SO(x = 1, y = 2, z = 3):  # Simple Orthorhombic
    axis = [[x, 0, 0], [0, y, 0], [0, 0, z]]
    return Lattice('SO', [[0, 0, 0]], axis)

def ABO(x = 1, y = 2, z = 3):  # Base-Centered Orthorhombic
    axis = [[x, 0, 0], [0, y, 0], [0, 0, z]]
    return Lattice('BCO', [[0, 0, 0], [0.5, 0.5, 0]], axis)

def FCO(x = 1, y = 2, z = 3):  # Face-Centered Orthorhombic
    axis = [[x, 0, 0], [0, y, 0], [0, 0, z]]
    return Lattice('FCO', [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]], axis)

def BCO(ratio_yx = 2, ratio_zx = 3):  # Body-Centered Orthorhombic
    axis = [[1, 0, 0], [0, ratio_yx, 0], [0, 0, ratio_zx]]
    return Lattice('ICO', [[0, 0, 0], [0.5, 0.5, 0.5]], axis)

# 单斜
def SM(x = 1, y = 2, z = 3, alpha_deg=100):  # Simple Monoclinic
    alpha = np.radians(alpha_deg)
    axis = [[x, 0, 0], [0, y, 0], [np.cos(alpha) * z, 0, np.sin(alpha) * z]]
    return Lattice('SM', [[0, 0, 0]], axis)

def BCM(x = 1, y = 2, z = 3, alpha_deg=100):  # Base-Centered Monoclinic
    alpha = np.radians(alpha_deg)
    axis = [[x, 0, 0], [0, y, 0], [np.cos(alpha) * z, 0, np.sin(alpha) * z]]
    return Lattice('BCM', [[0, 0, 0], [0.5, 0.5, 0]], axis)

# 三斜
def TRI(length = [1, 2, 3], angle = [80, 85, 75]):  # Triclinic
    alpha = np.radians(angle[0])
    beta = np.radians(angle[1])
    gamma = np.radians(angle[2])
    a = length[0]
    b = length[1]
    c = length[2]
    ax = [a, 0, 0]
    ay = [b * np.cos(gamma), b * np.sin(gamma), 0]
    cz_x = c * np.cos(beta)
    cz_y = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz_z = np.sqrt(c**2 - cz_x**2 - cz_y**2)
    az = [cz_x, cz_y, cz_z]
    return Lattice('TRI', [[0, 0, 0]], [ax, ay, az])

# 菱方
def RHL(alpha_deg=75):  # Rhombohedral (Trigonal)
    alpha = np.radians(alpha_deg)
    a1 = [1, 0, 0]
    a2 = [np.cos(alpha), np.sin(alpha), 0]
    a3 = [
        np.cos(alpha),
        (np.cos(alpha) - np.cos(alpha)**2)/(np.sin(alpha)),
        np.sqrt(1 - 2 * np.cos(alpha)**2 + np.cos(alpha)**3)/(np.sin(alpha))
    ]
    return Lattice('RHL', [[0, 0, 0]], [a1, a2, a3])

# 定义晶体结构复制
def Duplicate(lattice, nx=1, ny=1, nz=1):
    """
    Replicate a unit cell in the x, y, z directions.

    Parameters:
    - lattice: Lattice object with attributes .point (N,3), .axis (3,3), and .size (3,)
    - nx, ny, nz: Positive integers specifying the number of replications along each axis.

    Returns:
    - A new Lattice object with replicated atomic positions and updated size.
    """
    if not all(isinstance(n, int) and n > 0 for n in [nx, ny, nz]):
        raise ValueError("nx, ny, nz must be positive integers.")

    # Generate translation vectors in lattice units
    shifts = np.array([[i, j, k] for i in range(nx) for j in range(ny) for k in range(nz)])
    shifts = shifts * lattice.size  # convert to Cartesian translation in terms of unit cell

    # Add translations to each atom in the original lattice
    all_points = shifts[:, None, :] + lattice.point[None, :, :]
    all_points = all_points.reshape(-1, 3)

    # Update size
    new_size = lattice.size * np.array([nx, ny, nz])

    return Lattice(lattice.name, all_points, lattice.axis, new_size)

import numpy as np
import matplotlib.pyplot as plt

class Lattice_plot:
    def __init__(self, outfile='lattice', type_fig='svg'):
        """
        Initialize the LatticePlot.
        Args:
            outfile (str): Base filename for output files.
            type_fig (str): Vector format ('svg', 'pdf', 'eps', …).
        """
        self.outfile = outfile
        self.type_fig = type_fig

    def save_fig(self, fig, name):
        fig.savefig(
            f"{self.outfile}_{name}.{self.type_fig}",
            format=self.type_fig,
            bbox_inches='tight',     
            transparent=True,        
            pad_inches=0             
        )

    def _init_ax(self):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_axis_off()
        ax.view_init(elev=18, azim=39)
        return fig, ax

    def plot_body_center(self, a=1, b=1, c=1, save_path=None):
        # Define corner atoms
        corners = np.array([
            [0, 0, 0],
            [a, 0, 0],
            [a, b, 0],
            [0, b, 0],
            [0, 0, c],
            [a, 0, c],
            [a, b, c],
            [0, b, c]
        ])

        # Define body center atom
        body_center = np.array([[0.5 * a, 0.5 * b, 0.5 * c]])

        fig, ax = self._init_ax()

        # Plot corner atoms
        ax.scatter(corners[:,0], corners[:,1], corners[:,2], color='red', s=100)
        # Plot body center atom
        ax.scatter(body_center[:,0], body_center[:,1], body_center[:,2], color='red', s=100, alpha=0.5)

        # Draw edges of the cube
        for start in corners:
            for end in corners:
                dif = np.abs(start - end)
                if np.count_nonzero(dif) == 1 and (np.sum(dif) == a or np.sum(dif) == b or np.sum(dif) == c):
                    ax.plot(*zip(start, end), color='black', linewidth=1)

        # Draw dashed lines from corners to center
        for corner in corners:
            ax.plot(*zip(corner, body_center[0]), color='black', linestyle='--', linewidth=0.8)

        # Set axes properties
        ax.set_box_aspect([a,b,c])

        if save_path:
            self.save_fig(fig, save_path)
        plt.show()

    def plot_face_center(self, a=1, b=1, c=1, save_path=None):
        # Define corner atoms
        corners = np.array([
            [0, 0, 0],
            [a, 0, 0],
            [a, b, 0],
            [0, b, 0],
            [0, 0, c],
            [a, 0, c],
            [a, b, c],
            [0, b, c]
        ])

        # Define face-centered atoms
        faces = np.array([
            [0.5 * a, 0.5 * b, 0],
            [0.5 * a, 0.5 * b, c],
            [0.5 * a, 0, 0.5 * c],
            [0.5 * a, b, 0.5 * c],
            [0, 0.5 * b, 0.5 * c],
            [a, 0.5 * b, 0.5 * c]
        ])

        fig, ax = self._init_ax()

        # Plot corner atoms
        ax.scatter(corners[:,0], corners[:,1], corners[:,2], color='red', s=100)
        # Plot face-centered atoms
        ax.scatter(faces[:,0], faces[:,1], faces[:,2], color='red', s=100, alpha=0.5)

        # Draw edges of the cube
        for start in corners:
            for end in corners:
                diff = np.abs(start - end)
                if np.count_nonzero(diff) == 1 and (np.sum(diff) in [a, b, c]):
                    ax.plot(*zip(start, end), color='black', linewidth=1)

        # Define corners of each face
        face_corners = [
            [[0,0,0], [a,0,0], [a,b,0], [0,b,0]],      # bottom face z=0
            [[0,0,c], [a,0,c], [a,b,c], [0,b,c]],      # top face z=c
            [[0,0,0], [a,0,0], [a,0,c], [0,0,c]],      # front face y=0
            [[0,b,0], [a,b,0], [a,b,c], [0,b,c]],      # back face y=b
            [[0,0,0], [0,b,0], [0,b,c], [0,0,c]],      # left face x=0
            [[a,0,0], [a,b,0], [a,b,c], [a,0,c]]       # right face x=a
        ]

        # Draw dashed lines from face center atoms to corners of their face
        for fc_atom, fc_corners in zip(faces, face_corners):
            for corner in fc_corners:
                ax.plot(*zip(fc_atom, corner), color='black', linestyle='--', linewidth=0.8)

        # Set aspect ratio and hide axes
        ax.set_box_aspect([a,b,c])
        if save_path:
            self.save_fig(fig, save_path)
        plt.show()

    def plot_angle_arc(self, ax, u, v, angle, radius, label, center=np.zeros(3), error = [0.1, 0.1,0.1]):
        """
        Draws an arc between vectors u and v starting at `center`.

        Parameters:
        - ax: matplotlib 3D axis
        - u, v: vectors defining the angle plane
        - angle: angle between u and v (in radians)
        - radius: radius of the arc
        - label: text label (e.g., r'$\alpha$')
        - center: 3D coordinates of the arc center (default: [0, 0, 0])
        """
        u_hat = u / np.linalg.norm(u)
        v_hat = v / np.linalg.norm(v)
        n = np.cross(u_hat, v_hat)
        n /= np.linalg.norm(n)
        e1 = u_hat
        e2 = np.cross(n, u_hat)

        theta = np.linspace(0, angle, 100)
        arc = (radius * (np.outer(np.cos(theta), e1) +
                        np.outer(np.sin(theta), e2))) + center

        ax.plot(arc[:, 0], arc[:, 1], arc[:, 2], 'b-', linewidth=1.5)

        mid_theta = angle / 2
        mid_point = radius * (np.cos(mid_theta) * e1 + np.sin(mid_theta) * e2) + center + error
        ax.text(*mid_point, label, color='blue', fontsize=12, fontweight='bold')

    def plot_primitive(self, a=1, b=1, c=1, alpha=90, beta=90, gamma=90, save_path=None):
        alpha = np.radians(alpha)
        beta = np.radians(beta)
        gamma = np.radians(gamma)
        basis_vectors = np.array([
            [a, 0, 0],
            [b * np.cos(gamma), b * np.sin(gamma), 0],
            [c * np.cos(beta), c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
            c * np.sqrt(1 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma) -
                        np.cos(alpha) ** 2 - np.cos(beta) ** 2 - np.cos(gamma) ** 2) / np.sin(gamma)]
        ])

        corners = []
        for i in [0, 1]:
            for j in [0, 1]:
                for k in [0, 1]:
                    corner = i * basis_vectors[0] + j * basis_vectors[1] + k * basis_vectors[2]
                    corners.append(corner)
        corners = np.array(corners)

        fig, ax = self._init_ax()

        ax.scatter(corners[:, 0], corners[:, 1], corners[:, 2], color='red', s=100)

        edge_vectors = basis_vectors
        for i in range(len(corners)):
            for j in range(len(corners)):
                diff = corners[j] - corners[i]
                for vec in edge_vectors:
                    if np.allclose(diff, vec) or np.allclose(diff, -vec):
                        ax.plot(*zip(corners[i], corners[j]), color='black', linewidth=1)

        ax.text(basis_vectors[0, 0]/2, 0, 0.05, 'a', fontsize=14)
        mid_b = basis_vectors[1]/2
        ax.text(mid_b[0], mid_b[1], 0.05, 'b', fontsize=14)
        mid_c = basis_vectors[2]/2
        ax.text(mid_c[0]+0.1, mid_c[1], mid_c[2], 'c', fontsize=14)

        radius = min(a, b, c) * 0.2
        self.plot_angle_arc(ax, basis_vectors[1], basis_vectors[2], alpha, radius, r'$\alpha$', basis_vectors[0])
        self.plot_angle_arc(ax, basis_vectors[0], basis_vectors[2], beta, radius, r'$\beta$', basis_vectors[1], [0.2, 0.1, 0.1])
        self.plot_angle_arc(ax, basis_vectors[0], basis_vectors[1], gamma, radius, r'$\gamma$')

        # Determine the full box extent in x, y, z 
        x_range = corners[:, 0].max() - corners[:, 0].min()
        y_range = corners[:, 1].max() - corners[:, 1].min()
        z_range = corners[:, 2].max() - corners[:, 2].min()
        ax.set_box_aspect([x_range, y_range, z_range])
        if save_path:
            self.save_fig(fig, save_path)
        plt.show()

    def plot_base_center(self, a=1, b=1, c=1, alpha=90, beta=90, gamma=90, save_path=None):
        alpha = np.radians(alpha)
        beta = np.radians(beta)
        gamma = np.radians(gamma)

        # Lattice basis vectors
        basis_vectors = np.array([
            [a, 0, 0],
            [b * np.cos(gamma), b * np.sin(gamma), 0],
            [
                c * np.cos(beta),
                c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
                c * np.sqrt(1 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma)
                            - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2) / np.sin(gamma)
            ]
        ])

        # Corner atoms
        corners = []
        for i in [0, 1]:
            for j in [0, 1]:
                for k in [0, 1]:
                    corner = i * basis_vectors[0] + j * basis_vectors[1] + k * basis_vectors[2]
                    corners.append(corner)
        corners = np.array(corners)

        # Base-centered atoms: center of bottom and top faces
        base_atoms = []
        for k in [0, 1]:
            center = 0.5 * (basis_vectors[0] + basis_vectors[1]) + k * basis_vectors[2]
            base_atoms.append(center)
        base_atoms = np.array(base_atoms)

        fig, ax = self._init_ax()

        # Plot corner atoms
        ax.scatter(corners[:,0], corners[:,1], corners[:,2], color='red', s=100)

        # Plot base-centered atoms
        ax.scatter(base_atoms[:,0], base_atoms[:,1], base_atoms[:,2], color='red', s=100, alpha=0.5)

        # Draw edges
        for i in range(len(corners)):
            for j in range(i+1, len(corners)):
                diff = corners[j] - corners[i]
                for vec in basis_vectors:
                    if np.allclose(diff, vec) or np.allclose(diff, -vec):
                        ax.plot(*zip(corners[i], corners[j]), color='black', linewidth=1)

        # Draw dashed lines from base-centered atoms to their respective face corners
        bottom_face_corners = [
            corners[0],
            corners[6],
            corners[4],
            corners[2],
        ]
        top_face_corners = [
            corners[3],
            corners[5],
            corners[7],
            corners[1]
        ]

        for fc_atom, fc_corners in zip(base_atoms, [bottom_face_corners, top_face_corners]):
            for corner in fc_corners:
                ax.plot(*zip(fc_atom, corner), color='black', linestyle='--', linewidth=0.8)

        # Set aspect ratio and view
        x_range = corners[:,0].ptp()
        y_range = corners[:,1].ptp()
        z_range = corners[:,2].ptp()
        ax.set_box_aspect([x_range, y_range, z_range])
        ax.set_axis_off()
        if save_path:
            self.save_fig(fig, save_path)
        plt.show()

    def plot_Hexagonal_Close_Packed(self, c = 1.633, save_path=None):
        # Lattice constants
        a = 1.0
        c = c * a  # ideal HCP

        # Top and bottom layer hexagon positions
        hex_xy = np.array([
            [0, 0],
            [a, 0],
            [0.5*a, np.sqrt(3)/2*a],
            [-0.5*a, np.sqrt(3)/2*a],
            [-a, 0],
            [-0.5*a, -np.sqrt(3)/2*a],
            [0.5*a, -np.sqrt(3)/2*a]
        ])

        # Bottom layer (z=0)
        bottom = np.column_stack((hex_xy, np.zeros(7)))
        # Top layer (z=c)
        top = np.column_stack((hex_xy, c*np.ones(7)))

        # Middle-layer atoms at fractional positions inside triangles
        center_vector = np.array([2/3, 1/3, 0.5])
        center_vector = center_vector[0] * np.array([1, 0, 0]) + center_vector[1] * np.array([-0.5, np.sqrt(3)/2, 0]) + center_vector[2] * np.array([0, 0, c])
        middle = np.array([
            bottom[0] + center_vector,
            bottom[5] + center_vector,
            bottom[4] + center_vector,
        ])

        fig, ax = self._init_ax()

        # Plot bottom and top layers
        ax.scatter(bottom[:,0], bottom[:,1], bottom[:,2], color='red', s=100)
        ax.scatter(top[:,0], top[:,1], top[:,2], color='red', s=100)

        # Plot middle-layer atoms
        ax.scatter(middle[:,0], middle[:,1], middle[:,2], color='red', s=100, alpha=0.5)

        # Connect vertical edges
        for i in range(7):
            ax.plot(*zip(bottom[i], top[i]), color='black', linewidth=1)

        # Connect hexagon edges (bottom and top)
        bottom_c = bottom[0]
        top_c = top[0]
        for i in range(6):
            ax.plot(*zip(bottom_c, bottom[(i+1)]), color='black', linewidth=1)
            ax.plot(*zip(top_c, top[(i+1)]), color='black', linewidth=1)

        bottom = bottom[1:]
        top = top[1:]
        for i in range(6):
            ax.plot(*zip(bottom[i], bottom[(i+1)%6]), color='black', linewidth=1)
            ax.plot(*zip(top[i], top[(i+1)%6]), color='black', linewidth=1)

        # Connect middle edges
        for i in range(3):
            ax.plot(*zip(middle[i], middle[(i+1)%3]), color='black', linewidth=1)

        # Axes and view
        # Aspect ratio based on data range
        all_points = np.vstack((bottom, top, middle))
        ranges = all_points.max(axis=0) - all_points.min(axis=0)
        ax.set_box_aspect(ranges)
        if save_path:
            self.save_fig(fig, save_path)
        plt.show()

    def BCC(self): # Body-Centered Cubic
        self.plot_body_center(a=1, b=1, c=1, save_path='BCC')

    def BCT(self, c = 1.5):  # Body-Centered Tetragonal
        self.plot_body_center(c = c, save_path= "BCT")

    def BCO(self, a = 1.1, b =1.2, c= 1.3): # Body-Centered Orthorhombic
        self.plot_body_center(a, b, c, save_path= 'BCO')

    def FCC(self):  # Face-Centered Cubic
        self.plot_face_center(a=1, b=1, c=1, save_path= "FCC")

    def FCO(self, a=1.1, b=1.2, c=1.3):
        # Validate lattice type
        equal_count = sum([a == b, a == c, b == c])
        if equal_count == 1:
            raise ValueError(
                "Not a face-centered lattice: must be either FCC (a=b=c) or FCO (a≠b≠c)"
            )
        if equal_count == 3:
            self.FCC()
        else:
            self.plot_face_center(a, b, c, save_path= "FCO")

    def SC(self):  # Simple Cubic
        self.plot_primitive(a=1, b=1, c=1, save_path="SC")

    def ST(self, c=1.5):  # Simple Tetragonal
        if c == 1:
            print("It's actually Simple Cubic")
            self.SC()
        else:
            self.plot_primitive(a=1, b=1, c=c, save_path=f"ST")

    def SO(self, a=1.1, b=1.2, c=1.3):  # Simple Orthorhombic
        equal_count = sum([a == b, a == c, b == c])
        if equal_count >= 2:
            print("It's actually Simple Cubic")
            self.SC()
        elif equal_count == 1:
            print("It's actually Simple Tetragonal")
            name = f"ST"
            self.plot_primitive(a=a, b=b, c=c, save_path=name)
        else:
            name = f"SO"
            self.plot_primitive(a=a, b=b, c=c, save_path=name)

    def SM(self, a=1.1, b=1.2, c=1.3, alpha_deg=100):  # Simple Monoclinic
        self.plot_primitive(a=a, b=b, c=c, alpha=alpha_deg, beta=90, gamma=90, save_path=f"SM")

    def TRI(self, lengths=[1, 2, 3], angles=[80, 85, 75]):  # Triclinic
        a, b, c = lengths
        alpha, beta, gamma = angles
        self.plot_primitive(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma,save_path=f"TRI")

    def RHL(self, alpha_deg=75):  # Rhombohedral (Trigonal)
        self.plot_primitive(a=1, b=1, c=1, alpha=alpha_deg, beta=alpha_deg, gamma=alpha_deg, save_path=f"RHL")

    def ABO(self, a=1.1, b=1.2, c=1.3):  # Base-Centered Orthorhombic
        self.plot_base_center(a=a, b=b, c=c, save_path=f"ABO")

    def BCM(self, a=1.1, b=1.2, c=1.3, alpha_deg=100):  # Base-Centered Monoclinic
        # Monoclinic: only alpha != 90°
        self.plot_base_center(a=a, b=b, c=c, alpha=alpha_deg, beta=90, gamma=90, save_path=f"BCM")

    def HCP(self, c=1.633):  # Hexagonal Close Packed
        self.plot_Hexagonal_Close_Packed(c=c, save_path=f"HCP")

if __name__ == "__main__":
    lp = Lattice_plot(outfile="lattice_plot", type_fig="svg")

    print("Testing Body-Centered Cubic...")
    lp.BCC()

    print("Testing Body-Centered Tetragonal...")
    lp.BCT(c=1.5)

    print("Testing Body-Centered Orthorhombic...")
    lp.BCO(a=1.1, b=1.2, c=1.3)

    print("Testing Face-Centered Cubic...")
    lp.FCC()

    print("Testing Face-Centered Orthorhombic...")
    lp.FCO(a=1.1, b=1.2, c=1.3)

    print("Testing Simple Cubic...")
    lp.SC()

    print("Testing Simple Tetragonal...")
    lp.ST(c=1.5)

    print("Testing Simple Orthorhombic...")
    lp.SO(a=1.1, b=1.2, c=1.3)

    print("Testing Simple Monoclinic...")
    lp.SM(a=1.1, b=1.2, c=1.3, alpha_deg=100)

    print("Testing Triclinic...")
    lp.TRI(lengths=[1, 1.2, 1.4], angles=[80, 85, 75])

    print("Testing Rhombohedral...")
    lp.RHL(alpha_deg=75)

    print("Testing Base-Centered Orthorhombic...")
    lp.ABO(a=1.1, b=1.2, c=1.3)

    print("Testing Base-Centered Monoclinic...")
    lp.BCM(a=1.1, b=1.2, c=1.3, alpha_deg=100)

    print("Testing Hexagonal Close Packed...")
    lp.HCP(c=1.633)

    print("All lattice types tested.")

