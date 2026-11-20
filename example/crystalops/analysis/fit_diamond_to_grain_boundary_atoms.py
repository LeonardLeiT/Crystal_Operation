"""Fit an exact one-cell Diamond template to selected grain-boundary atoms."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from crystalops.analysis import identify_grain_boundary_atoms, mean_curvature_from_derivatives
from crystalops.analysis.diamond_fit import diamond_distances, fit_diamond_phase, sample_diamond_template
from crystalops.analysis.diamond_reference import exact_diamond_field
from crystalops.io import read_structure
from example.crystalops.analysis.analyze_grain_boundary_curvature import load_grain_labels_csv


def _distance_statistics(distances, weights=None):
    if not len(distances):
        return None
    if weights is None:
        weights = np.ones(len(distances))
    order = np.argsort(distances)
    fractions = np.cumsum(weights[order]) / np.sum(weights)
    median, p90 = distances[order][np.searchsorted(fractions, [0.5, 0.9])]
    return {
        "median_A": float(median), "p90_A": float(p90),
        "fraction_within_2p5_A": float(np.average(distances <= 2.5, weights=weights)),
        "fraction_within_5_A": float(np.average(distances <= 5, weights=weights)),
    }


def _fit_subset(template, positions, atoms, block_ids, holdout_blocks, pair_labels):
    training = ~np.isin(block_ids[atoms], holdout_blocks)
    if training.sum() < 20 or training.all():
        raise ValueError("the selected subset needs at least 20 training atoms and one held-out atom")
    phase, fit = fit_diamond_phase(positions[atoms[training]], template)
    d, closest, stationary = diamond_distances(positions[atoms], template, phase)
    centers = (template.centers + phase * template.lengths) % template.lengths
    reverse, _ = cKDTree(positions[atoms], boxsize=template.lengths).query(centers)

    # Random phases use the same sampled-distance estimator as the fitted
    # comparison, so the control is not favored by closest-point refinement.
    holdout_positions = positions[atoms[~training]]
    sampled_d, _ = template.tree.query(holdout_positions - phase * template.lengths)
    random_phases = np.random.default_rng(31).uniform(0, 0.5, size=(32, 3))
    random_support = [float(np.mean(template.tree.query(holdout_positions - shift * template.lengths)[0] <= 2.5))
                      for shift in random_phases]
    pairs = {}
    for pair in np.unique(pair_labels[atoms], axis=0):
        included = np.all(pair_labels[atoms] == pair, axis=1)
        key = f"{pair[0]}-{pair[1]}"
        pairs[key] = {
            "atom_count": int(included.sum()),
            "all_atoms": _distance_statistics(d[included]),
            "holdout_atoms": _distance_statistics(d[included & ~training]),
        }

    report = {
        "phase_fractional": phase.tolist(), "translation_relative_to_box_origin_A": (phase * template.lengths).tolist(),
        "optimizer": fit,
        "training_count": int(training.sum()), "holdout_count": int((~training).sum()),
        "all_atoms_to_template": _distance_statistics(d),
        "training_atoms_to_template": _distance_statistics(d[training]),
        "holdout_atoms_to_template": _distance_statistics(d[~training]),
        "template_to_atoms_area_weighted": _distance_statistics(reverse, template.areas),
        "closest_point_stationary_candidate_fraction": float(np.mean(stationary)),
        "holdout_sampled_distance_phase_control_within_2p5_A": {
            "fit": float(np.mean(sampled_d <= 2.5)),
            "random_median": float(np.median(random_support)),
            "random_maximum": max(random_support), "random_phase_count": len(random_support),
        },
        "grain_pairs": pairs,
    }
    return report, d, closest, ~training


def evaluate_structure(filename, labels_file, output_dir, *, cutoff=3.2, central_grain_id=None, resolution=64):
    structure = read_structure(filename)
    lengths = np.diag(structure.cell)
    if np.any(lengths <= 0) or not np.allclose(structure.cell, np.diag(lengths)) or not all(structure.pbc):
        raise ValueError("this aligned template comparison requires an orthogonal xyz-periodic Cartesian box")
    labels = load_grain_labels_csv(structure, labels_file)
    selection = identify_grain_boundary_atoms(structure, neighbor_cutoff=cutoff, grain_ids=labels)
    positions = (structure.positions - structure.origin) % lengths
    pair_labels = np.zeros((len(labels), 2), dtype=int)
    for pair, indices in selection.pair_atom_indices.items():
        pair_labels[indices[~selection.junction_mask[indices]]] = pair
    eligible = selection.boundary_mask & ~selection.junction_mask
    if not np.all(pair_labels[eligible] > 0):
        raise ValueError("expected positive grain pairs for all two-grain boundary atoms")
    blocks = np.minimum(3, np.floor(4 * positions / lengths).astype(int))
    block_ids = np.ravel_multi_index(blocks.T, (4, 4, 4))
    holdout_blocks = np.random.default_rng(14).choice(64, 16, replace=False)

    template = sample_diamond_template(lengths, resolution=resolution)
    _, gradient, hessian = exact_diamond_field(template.centers, lengths, order=2)
    h = mean_curvature_from_derivatives(gradient, hessian)
    quadrature = {
        "area_A2": float(template.areas.sum()),
        "mean_absolute_H_per_A": float(np.average(np.abs(h), weights=template.areas)),
        "rms_H_per_A": float(np.sqrt(np.average(h*h, weights=template.areas))),
    }
    subsets = {"all_two_grain_contacts": np.flatnonzero(eligible)}
    if central_grain_id is not None:
        if central_grain_id not in labels:
            raise ValueError("central_grain_id is not present in the provided labels")
        subsets["outer_grain_contacts"] = np.flatnonzero(eligible & ~(pair_labels == central_grain_id).any(axis=1))

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "input": Path(filename).name, "input_sha256": hashlib.sha256(Path(filename).read_bytes()).hexdigest(),
        "labels_sha256": hashlib.sha256(Path(labels_file).read_bytes()).hexdigest(),
        "box_lengths_A": lengths.tolist(), "box_origin_A": structure.origin.tolist(),
        "neighbor_cutoff_A": cutoff, "total_atoms": len(labels), "excluded_junction_count": int(selection.junction_mask.sum()),
        "central_grain_id": central_grain_id,
        "template": "one full exact Schwarz Diamond periodic cell, scaled separately along the Cartesian box axes",
        "fit_parameters": "three translations only, fit to selected atoms outside 16 held-out spatial blocks",
        "template_sampling_resolution": resolution, "holdout_blocks": sorted(holdout_blocks.tolist()),
        "conditional_template_curvature": quadrature, "subsets": {},
    }
    ids = structure.properties.get("atom_ids", np.arange(len(labels)) + 1)
    for name, atoms in subsets.items():
        fit, distances, closest, holdout = _fit_subset(template, positions, atoms, block_ids, holdout_blocks, pair_labels)
        report["subsets"][name] = fit
        np.savez_compressed(
            output_dir / f"{name}.npz", atom_id=ids[atoms], atom_xyz=structure.positions[atoms],
            closest_xyz=closest + structure.origin, distance=distances, holdout=holdout, grain_pair=pair_labels[atoms],
            cell=structure.cell, origin=structure.origin, phase_fractional=np.asarray(fit["phase_fractional"]),
        )
    path = output_dir / "diamond_fit.json"
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--grain-labels", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--neighbor-cutoff", type=float, default=3.2)
    parser.add_argument("--central-grain-id", type=int)
    parser.add_argument("--resolution", type=int, default=64)
    args = parser.parse_args()
    report = evaluate_structure(args.input, args.grain_labels, args.output_dir, cutoff=args.neighbor_cutoff,
                                central_grain_id=args.central_grain_id, resolution=args.resolution)
    for name, result in report["subsets"].items():
        print(name, json.dumps({"atom_to_surface": result["holdout_atoms_to_template"],
                                "surface_to_atom": result["template_to_atoms_area_weighted"]}))


if __name__ == "__main__":
    main()
