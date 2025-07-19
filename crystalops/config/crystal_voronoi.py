import numpy as np
from scipy.spatial import Voronoi

def compute_voronoi_face_centers(points, box_vectors = np.eye(3), center_index=0, layers=1, center_surface = True):
    real_points = points @ box_vectors
    shifts = np.arange(-layers, layers+1)
    translations = np.array([[i, j, k] for i in shifts for j in shifts for k in shifts])
    real_translations = translations @ box_vectors
    all_points = (real_points[:, None, :] + real_translations[None, :, :]).reshape(-1, 3)
    vor = Voronoi(all_points)
    
    num_images = len(translations)
    zero_image = np.where((translations == np.array([0, 0, 0])).all(axis=1))[0][0]
    center_idx = center_index * num_images + zero_image
    
    # Identify the region of the center point
    region_index = vor.point_region[center_idx]
    region = vor.regions[region_index]
    
    if -1 in region:
        raise ValueError("Voronoi region is unbounded; consider increasing layers.")
    if center_surface:
        # List to store the face centers
        face_centers = []
        
        # Iterate over the ridge points to find the face centers
        for ridge_points, ridge_vertices in zip(vor.ridge_points, vor.ridge_vertices):
            if center_idx in ridge_points and all(v >= 0 for v in ridge_vertices):
                p1 = all_points[ridge_points[0]]
                p2 = all_points[ridge_points[1]]
                center = (p1 + p2) / 2
                face_centers.append(center)
        # print(real_points[center_index])
        face_centers = np.array(face_centers) - real_points[center_index]
        inv_box = np.linalg.inv(box_vectors) 
        return face_centers, face_centers @ inv_box
    else:
        vertices = vor.vertices[region] - real_points[center_index]
        
        # 收集所有和中心点有关的 ridge（边界面）对应的顶点索引
        faces = []
        for (p1, p2), ridge_vertices in zip(vor.ridge_points, vor.ridge_vertices):
            if center_idx in (p1, p2) and -1 not in ridge_vertices:
                # 找到这个面对应于中心点的顶点索引在 region 内的相对位置
                try:
                    face = [region.index(v) for v in ridge_vertices]
                    if len(face) >= 3:  # 至少是一个面
                        faces.append(face)
                except ValueError:
                    continue  # 某个 ridge vertex 不在 region 中，跳过
        return vertices, faces
    
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

def plot_voronoi_polyhedron(vertices, faces, color='skyblue', edge_color='k', alpha=0.5, show_vertices=True, save_path=None):
    """
    Visualize a Voronoi polyhedron with better 3D effect.

    Parameters:
        vertices (ndarray): (N, 3) vertex coordinates relative to the center atom
        faces (List[List[int]]): list of face vertex indices
        color (str): face base color
        edge_color (str): edge color
        alpha (float): face transparency
        show_vertices (bool): whether to show vertex spheres
    """
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=20, azim=-16)
    # Build each face
    poly3d = [[vertices[idx] for idx in face] for face in faces]
    collection = Poly3DCollection(poly3d, facecolors=color, edgecolors=edge_color, alpha=alpha)
    ax.add_collection3d(collection)

    # Voronoi vertex spheres
    if show_vertices:
        ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2],
                   c='red', s=100, edgecolors='black', linewidths=0.5, depthshade=True)

    vertices_all = np.vstack([vertices, [[0, 0, 0]]])
    max_range = (vertices_all.max(axis=0) - vertices_all.min(axis=0)).max() / 2.0
    mid = vertices_all.mean(axis=0)
    for axis, m in zip([ax.set_xlim, ax.set_ylim, ax.set_zlim], mid):
        axis([m - max_range, m + max_range])

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect([1, 1, 1])
    ax.set_axis_off()

    if save_path:
        fig.savefig(
            f"{save_path}_voronoi.svg",
            format='svg',
            bbox_inches='tight',     
            transparent=True,        
            pad_inches=0             
        )
    plt.show()