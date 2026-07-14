import numpy as np
Plot = False


def fit_plane_normal(points, reference_direction=None):
    """
    Least-squares best-fit plane normal via SVD.

    Parameters
    ----------
    points : np.ndarray of shape (N, 3)
        Coordinates of N >= 3 points lying near a common plane.
    reference_direction : np.ndarray of shape (3,) or None
        Optional reference direction (e.g., detected symmetry axis or ROI
        direction). If provided, the normal is oriented to point in the
        same half-space (dot product > 0). If None, falls back to a
        deterministic sign rule (+Z, then +Y, then +X).

    Returns
    -------
    normal : np.ndarray of shape (3,)
        Unit normal vector of the best-fit plane (vh[-1] from SVD).
    rmsd : float
        Root-mean-square deviation of points from the fitted plane (Angstrom).
    singular_values : np.ndarray of shape (3,)
        Singular values s[0] >= s[1] >= s[2] from the SVD of centered points.
        s[2] relates to RMSD; s[1]/s[0] indicates in-plane point spread
        (ratio &lt;&lt; 1 implies near-collinear points, making the normal unreliable).

    Notes
    -----
    Unlike ``definePlane()``, this returns a **normalized** normal vector.
    For N=3 exactly coplanar points, rmsd will be 0.0.
    """
    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) < 3:
        raise ValueError(
            f"Expected an array of shape (N, 3) with N >= 3, got {points.shape}."
        )
    center = points.mean(axis=0)
    _, s, vh = np.linalg.svd(points - center)
    normal = vh[-1]
    # Orient normal: prefer reference direction, then +Z, +Y, +X
    if reference_direction is not None:
        reference_direction = np.asarray(reference_direction, dtype=np.float64)
        ref_norm = np.linalg.norm(reference_direction)
        if ref_norm < 1e-12:
            raise ValueError("reference_direction must be nonzero.")
        if np.dot(normal, reference_direction / ref_norm) < 0:
            normal = -normal
    else:
        for axis in (2, 1, 0):
            if abs(normal[axis]) > 1e-12:
                if normal[axis] < 0:
                    normal = -normal
                break
    residuals = (points - center) @ normal
    rmsd = np.sqrt(np.mean(residuals ** 2))
    return normal, rmsd, s


def select_5fold_ring(copies, pointA):
    """
    Select the five distance-ranked copies forming a 5-fold local ring.

    Parameters
    ----------
    copies : np.ndarray of shape (M, 3) with M >= 5
        Coordinates of all symmetry-related copies of the reference atom.
    pointA : np.ndarray of shape (3,)
        Coordinates of the reference copy.

    Returns
    -------
    ring : np.ndarray of shape (5, 3)
        The five copies sorted by increasing distance from pointA.

    Raises
    ------
    ValueError
        If fewer than 5 copies are provided.
    """
    if len(copies) < 5:
        raise ValueError(
            f"Need at least 5 copies for 5-fold ring selection, got {len(copies)}."
        )
    copies = np.asarray(copies, dtype=np.float64)
    pointA = np.asarray(pointA, dtype=np.float64)
    dists = np.linalg.norm(copies - pointA, axis=1)
    order = np.argsort(dists)
    return copies[order[:5]]


def compute_pentagon_diagnostics(ring, global_axis=None):
    """
    Lightweight pentagon quality check for a 5-fold ring.

    Computes plane RMSD, adjacent-edge distance variation, and the
    diagonal-to-side ratio which should approach the golden ratio
    φ ≈ 1.618 for a regular pentagon.

    Parameters
    ----------
    ring : np.ndarray of shape (5, 3)
        Five points forming a candidate 5-fold ring.
    global_axis : np.ndarray of shape (3,) or None
        Optional reference axis for alignment comparison.

    Returns
    -------
    dict
        Keys: plane_rmsd, non_collinearity, adj_dist_mean, adj_dist_std,
        diag_mean, diag_adj_ratio, phi, phi_deviation.
        Includes alignment_angle_deg when *global_axis* is provided.
    """
    ring = np.asarray(ring, dtype=np.float64)
    if ring.shape != (5, 3):
        raise ValueError(
            f"Expected ring of shape (5, 3), got {ring.shape}."
        )
    # SVD plane fit — same logic as fit_plane_normal() but we also need vh
    center = ring.mean(axis=0)
    _, s, vh = np.linalg.svd(ring - center)
    normal = vh[-1]
    rmsd = np.sqrt(np.mean(((ring - center) @ normal) ** 2))
    ncl = s[1] / s[0] if s[0] > 1e-12 else 0.0

    # Sort points by angular position in the best-fit plane
    in_plane = (ring - center) @ vh[:2].T  # (5, 2)
    angles = np.arctan2(in_plane[:, 1], in_plane[:, 0])
    order = np.argsort(angles)
    ring_ordered = ring[order]

    # Adjacent edge distances (cyclic)
    adj_dists = np.linalg.norm(
        ring_ordered - np.roll(ring_ordered, -1, axis=0), axis=1,
    )
    adj_mean = float(np.mean(adj_dists))
    adj_std = float(np.std(adj_dists))

    # Second-neighbour / diagonal distances (skip one, cyclic)
    diag_dists = np.linalg.norm(
        ring_ordered - np.roll(ring_ordered, -2, axis=0), axis=1,
    )
    diag_mean = float(np.mean(diag_dists))

    phi = (1 + np.sqrt(5)) / 2
    diag_adj_ratio = diag_mean / adj_mean if adj_mean > 1e-12 else 0.0
    phi_deviation = abs(diag_adj_ratio - phi) / phi

    result = {
        'plane_rmsd': float(rmsd),
        'non_collinearity': float(ncl),
        'adj_dist_mean': adj_mean,
        'adj_dist_std': adj_std,
        'diag_mean': diag_mean,
        'diag_adj_ratio': diag_adj_ratio,
        'phi': phi,
        'phi_deviation': phi_deviation,
    }

    if global_axis is not None:
        gax = np.asarray(global_axis, dtype=np.float64)
        gax = gax / np.linalg.norm(gax)
        dot = np.clip(np.abs(np.dot(normal, gax)), -1.0, 1.0)
        result['alignment_angle_deg'] = float(np.rad2deg(np.arccos(dot)))

    return result


def definePlane(x,y,z):
    coords = np.vstack((x,y,z)).T
    com = (np.mean(x), np.mean(y), np.mean(z))
    coords = coords - com

    x=coords[:,0]
    y=coords[:,1]
    z=coords[:,2]

    if Plot:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        ax.scatter(x, y, z, color='red')
        ax.scatter(0,0,0, color='goldenrod')

    ux, uy, uz = u = [x[1], y[1], z[1]]
    vx, vy, vz = v = [x[2], y[2], z[2]]
    # u_cross_v = [uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx]
    u_cross_v = np.cross(u,v)

    point = np.zeros(3)
    normal = np.array(u_cross_v)
    d = -point.dot(normal)

    if Plot:
        if abs(normal[2]) > 1e-12:
            xx, yy = np.meshgrid((-10,10), (-10,10))
            z = (-normal[0]*xx - normal[1]*yy - d) * 1 / normal[2]
            ax.plot_surface(xx,yy,z, alpha=0.6)
        else:
            print("[definePlane] Skipping plane surface plot because normal[2] is near zero.")
    # ax.scatter(normal[0], normal[1], normal[2], color='indigo')
    normalVectorLine = np.vstack((point,normal))

    normalVector = np.asarray((point - normal), dtype=float)

    if Plot:
        ax.plot(normalVectorLine[:,0], normalVectorLine[:,1], normalVectorLine[:,2], color='indigo')
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")

        lim=80
        ax.set_xlim(-1*lim,lim)
        ax.set_ylim(-1*lim,lim)
        ax.set_zlim(-1*lim,lim)
        plt.show()

    return normal
