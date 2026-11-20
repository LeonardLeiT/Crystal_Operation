"""Periodic connectivity of the four DSC transition regions."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def _periodic_components(mask: NDArray[np.bool_]) -> tuple[np.ndarray, np.ndarray]:
    """Label face-connected voxels, identifying opposite cell boundaries."""
    from scipy.ndimage import label

    components, count = label(mask)
    parents = np.arange(count + 1)

    def root(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = int(parents[index])
        return index

    for axis in range(3):
        lower = components.take(0, axis=axis)
        upper = components.take(-1, axis=axis)
        connected = (lower != 0) & (upper != 0) & (lower != upper)
        pairs = np.unique(np.column_stack((lower[connected], upper[connected])), axis=0)
        for first, second in pairs:
            left, right = root(int(first)), root(int(second))
            parents[max(left, right)] = min(left, right)

    periodic_ids = np.fromiter((root(index) for index in range(count + 1)), dtype=np.intp)
    components = periodic_ids[components]
    sizes = np.bincount(components.ravel(), minlength=count + 1)
    sizes[0] = 0
    return components, sizes


def reassign_dsc_islands(
    owners: np.ndarray,
    eligible: NDArray[np.bool_],
    protected: NDArray[np.bool_],
) -> tuple[np.ndarray, int]:
    """Attach detached twin voxels to an eligible neighboring main component.

    Labels 1 and 2 are the central regions, labels 3 through 6 are the
    four transition regions, and 0 denotes unfilled space. Eligibility
    channel 0 corresponds to owner 3, and so on. Eligibility describes each
    physical twin slab. A detached component can change owners only if
    all its voxels lie inside the destination slab and none are protected.
    Choose the destination with the most shared voxel faces against its
    original largest component, breaking ties by smaller owner label.
    """
    owners = np.asarray(owners)
    eligible = np.asarray(eligible, dtype=bool)
    protected = np.asarray(protected, dtype=bool)
    if (
        owners.ndim != 3
        or not all(owners.shape)
        or not np.issubdtype(owners.dtype, np.integer)
        or np.any((owners < 0) | (owners > 6))
    ):
        raise ValueError("owners must be a nonempty 3D integer grid with labels 0 through 6")
    if eligible.shape != (*owners.shape, 4) or protected.shape != owners.shape:
        raise ValueError("eligible must have four channels and protected must match the ownership grid")

    output = np.array(owners, copy=True, order="C")
    main_owners = np.zeros(owners.shape, dtype=np.int8)
    islands = []
    for owner in range(3, 7):
        components, sizes = _periodic_components(owners == owner)
        if not np.any(sizes):
            continue
        largest = int(sizes.argmax())
        main_owners[components == largest] = owner
        by_component = np.argsort(components.ravel(), kind="stable")
        stops = np.cumsum(np.bincount(components.ravel(), minlength=len(sizes)))
        for component in np.flatnonzero(sizes):
            if component != largest:
                islands.append(by_component[stops[component - 1]:stops[component]])
    if not islands:
        return output, 0

    contact_owners = [
        np.roll(main_owners, shift, axis=axis).ravel()
        for axis in range(3)
        for shift in (-1, 1)
    ]
    eligibility = eligible.reshape(-1, 4)
    protected_flat = protected.ravel()
    output_flat = output.reshape(-1)
    changed = 0
    for indices in islands:
        if protected_flat[indices].any():
            continue
        allowed = eligibility[indices].all(axis=0)
        if not allowed.any():
            continue
        neighbors = np.concatenate([contacts[indices] for contacts in contact_owners])
        contacts = np.bincount(neighbors, minlength=7)[3:]
        contacts[~allowed] = 0
        if contacts.any():
            output_flat[indices] = 3 + int(contacts.argmax())
            changed += len(indices)
    return output, changed
