import numpy as np

def aniso_alpha(n: np.ndarray) -> float:
    """
    Calculate anisotropy alpha from an averaged 3D vector n.

    Args:
        n (np.ndarray): Shape (3,), averaged displacement vector.

    Returns:
        float: Vibrational anisotropy alpha.
    """
    # Handle zero vector case
    if np.allclose(n, np.zeros(3)):
        return 0.0
    
    # Fabric tensor from mean vector
    F = np.outer(n, n)
    F = F.astype(np.float64)
    F /= np.trace(F)  # normalize to trace=1

    eigvals = np.linalg.eigvalsh(F)
    print(eigvals)

    alpha = 3 / np.sqrt(6) * np.sqrt(np.sum((eigvals - 1/3)**2))
    return alpha

def aniso_lattice(n, a0=1.0):
    """
    Compute the lattice deviation from bcc.

    Args:
        n (array-like): [a, b, c] lattice constants of the distorted structure
        a0 (float, optional): reference bcc lattice constant. Default is 1.0

    Returns:
        float: Δ_lattice
    """
    n = np.array(n, dtype=float)
    delta = 100
    for i in range(3):
        n1 = n / n[i]      # normalize so that axis i = 1
        ratios = n1 / a0
        norm_dev = np.linalg.norm(ratios - 1.0)
        print(f"{norm_dev:.6f}")
        delta = min(delta, norm_dev)
    return delta


if __name__ == "__main__":
    n_zero = np.array([1, 0, 0])
    print(f"Alpha: {aniso_alpha(n_zero)}")
    
    n = [1, 1.1, 1.2]
    print(f"Δ_lattice = {aniso_lattice(n):.6f}")
