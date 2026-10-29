"""Select reproducible SigmaZ2 twin sets from the periodic face graph."""

from __future__ import annotations

import itertools
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
from scipy.spatial import Voronoi
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parent
SIGMA_SOURCE = ROOT.parent / "Sigma"
if str(SIGMA_SOURCE) not in sys.path:
    sys.path.insert(0, str(SIGMA_SOURCE))

from build_sigma import (
    C_OVER_A,
    EQUAL_VOLUME_WEIGHTS,
    GROUP_IDS,
    SIGMA_CENTERS,
    SITE_CLASSES,
    TWIN_FACE_PAIRS,
)


DEFAULT_CONFIG = ROOT / "twin_designs.json"
Z2_CENTERS = np.concatenate(
    (
        SIGMA_CENTERS * np.asarray([1.0, 1.0, 0.5]),
        SIGMA_CENTERS * np.asarray([1.0, 1.0, 0.5])
        + np.asarray([0.0, 0.0, 0.5]),
    )
)
Z2_CLASSES = SITE_CLASSES + SITE_CLASSES
Z2_WEIGHTS = EQUAL_VOLUME_WEIGHTS[np.tile(GROUP_IDS, 2)]


def _laguerre_face_margin(first: int, second: int, vector: np.ndarray) -> float:
    """Return the largest normalized interior margin on a weighted common face."""
    cell = np.diag([1.0, 1.0, 2.0 * C_OVER_A])
    shifts = np.asarray(
        list(itertools.product((-1.0, 0.0, 1.0), repeat=3)), dtype=np.float64
    )
    selected = vector @ cell
    selected_limit = float(
        selected @ selected + Z2_WEIGHTS[first] - Z2_WEIGHTS[second]
    )
    rows: list[list[float]] = []
    limits: list[float] = []
    for neighbor_index, neighbor in enumerate(Z2_CENTERS):
        for shift in shifts:
            if neighbor_index == first and np.all(shift == 0.0):
                continue
            relative = (neighbor + shift - Z2_CENTERS[first]) @ cell
            if neighbor_index == second and np.linalg.norm(relative - selected) < 1e-9:
                continue
            normal = 2.0 * relative
            norm = float(np.linalg.norm(normal))
            limit = float(
                relative @ relative
                + Z2_WEIGHTS[first]
                - Z2_WEIGHTS[neighbor_index]
            )
            rows.append([*(normal / norm), 1.0])
            limits.append(limit / norm)
    result = linprog(
        [0.0, 0.0, 0.0, -1.0],
        A_ub=np.asarray(rows),
        b_ub=np.asarray(limits),
        A_eq=[[*(2.0 * selected), 0.0]],
        b_eq=[selected_limit],
        bounds=[(None, None)] * 4,
        method="highs",
    )
    return float(result.x[3]) if result.success else float("-inf")


def _periodic_face_edges() -> list[dict[str, Any]]:
    """Return canonical true-face edges excluding every original Sigma pair."""
    cell = np.diag([1.0, 1.0, 2.0 * C_OVER_A])
    shifts = np.asarray(
        list(itertools.product((-1.0, 0.0, 1.0), repeat=3)), dtype=np.float64
    )
    points: list[np.ndarray] = []
    labels: list[tuple[int, tuple[int, int, int]]] = []
    for shift in shifts:
        for index, center in enumerate(Z2_CENTERS):
            points.append((center + shift) @ cell)
            labels.append((index, tuple(int(value) for value in shift)))
    tessellation = Voronoi(np.asarray(points))
    canonical: dict[tuple[int, ...], tuple[int, int, np.ndarray]] = {}
    for (point_a, point_b), vertices in zip(
        tessellation.ridge_points, tessellation.ridge_vertices, strict=True
    ):
        if -1 in vertices:
            continue
        for source, neighbor in ((point_a, point_b), (point_b, point_a)):
            first, first_shift = labels[source]
            second, second_shift = labels[neighbor]
            if first_shift != (0, 0, 0):
                continue
            shift = np.asarray(second_shift, dtype=np.int64)
            forward = (first, second, *shift.tolist())
            reverse = (second, first, *(-shift).tolist())
            if reverse < forward:
                key = reverse
                canonical[key] = (
                    second,
                    first,
                    Z2_CENTERS[first] - shift - Z2_CENTERS[second],
                )
            else:
                key = forward
                canonical[key] = (
                    first,
                    second,
                    Z2_CENTERS[second] + shift - Z2_CENTERS[first],
                )

    inherited = {frozenset((first, second)) for first, second, _ in TWIN_FACE_PAIRS}
    result: list[dict[str, Any]] = []
    for key in sorted(canonical):
        first, second, vector = canonical[key]
        if frozenset((first % 30, second % 30)) in inherited:
            continue
        if first < 30 and second < 30:
            placement = "lower"
        elif first >= 30 and second >= 30:
            placement = "upper"
        else:
            placement = "cross-layer"
        pair_class = "-".join(sorted((Z2_CLASSES[first], Z2_CLASSES[second])))
        normal = vector @ cell
        normal /= np.linalg.norm(normal)
        direction_bin = tuple(np.rint(np.abs(normal) * 4.0).astype(int).tolist())
        result.append(
            {
                "first": first,
                "second": second,
                "vector": vector,
                "placement": placement,
                "pair_class": pair_class,
                "direction_bin": direction_bin,
            }
        )
    return result


def _class_matches(edge: dict[str, Any], class_filter: str) -> bool:
    if class_filter == "any":
        return True
    if class_filter == "contains-Z12":
        return "Z12" in edge["pair_class"]
    if class_filter == "contains-Z15":
        return "Z15" in edge["pair_class"]
    return edge["pair_class"] == class_filter


def select_design(
    specification: dict[str, Any],
    *,
    selection_trials: int,
    edges: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Select a disjoint, diverse matching satisfying one design specification."""
    candidates = [
        edge
        for edge in (_periodic_face_edges() if edges is None else edges)
        if _class_matches(edge, specification.get("class_filter", "any"))
    ]
    quotas = specification.get("placements")
    target = int(specification["target_pairs"])
    rng = np.random.default_rng(int(specification["selection_seed"]))
    best: list[dict[str, Any]] | None = None
    best_score: tuple[int, int, int, int] | None = None
    for _ in range(selection_trials):
        used: set[int] = set()
        selected: list[dict[str, Any]] = []
        counts = {"lower": 0, "upper": 0, "cross-layer": 0}
        for candidate_index in rng.permutation(len(candidates)):
            edge = candidates[int(candidate_index)]
            if edge["first"] in used or edge["second"] in used:
                continue
            if quotas and counts[edge["placement"]] >= int(
                quotas.get(edge["placement"], 0)
            ):
                continue
            selected.append(edge)
            used.update((edge["first"], edge["second"]))
            counts[edge["placement"]] += 1
            if len(selected) == target:
                break
        exact = len(selected) == target and (
            not quotas
            or all(counts[placement] == int(value) for placement, value in quotas.items())
        )
        diversity = len({edge["pair_class"] for edge in selected}) + len(
            {edge["direction_bin"] for edge in selected}
        )
        score = (int(exact), len(selected), min(counts.values()), diversity)
        if best_score is None or score > best_score:
            best_score = score
            best = selected
    if best is None or best_score is None or best_score[0] != 1:
        raise RuntimeError(
            f"cannot satisfy SigmaZ2 design {specification['id']}: {best_score}"
        )
    pairs = [
        (
            int(edge["first"]),
            int(edge["second"]),
            tuple(float(value) for value in edge["vector"]),
            str(edge["placement"]),
        )
        for edge in best
    ]
    placement_counts = {
        placement: sum(edge["placement"] == placement for edge in best)
        for placement in ("lower", "upper", "cross-layer")
    }
    class_counts: dict[str, int] = {}
    for edge in best:
        class_counts[edge["pair_class"]] = class_counts.get(edge["pair_class"], 0) + 1
    face_margins = [
        _laguerre_face_margin(edge["first"], edge["second"], edge["vector"])
        for edge in best
    ]
    if min(face_margins) <= 1e-6:
        raise RuntimeError(
            f"selected pair in {specification['id']} is not a finite weighted Laguerre face"
        )
    return {
        **specification,
        "pairs": pairs,
        "placement_counts": placement_counts,
        "pair_class_counts": class_counts,
        "ordinary_grain_count": 60 - 2 * len(pairs),
        "inherits_sigma_pair_count": 0,
        "minimum_laguerre_face_margin": min(face_margins),
    }


def load_design_library(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1:
        raise ValueError("unsupported SigmaZ2 design schema")
    edges = _periodic_face_edges()
    trials = int(payload["selection_trials"])
    designs = [
        select_design(specification, selection_trials=trials, edges=edges)
        for specification in payload["designs"]
    ]
    return {**payload, "selected_designs": designs, "available_edge_count": len(edges)}


def load_design(design_id: str, path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    library = load_design_library(path)
    for design in library["selected_designs"]:
        if design["id"] == design_id:
            return design
    raise ValueError(f"unknown SigmaZ2 design {design_id!r}")


__all__ = ["DEFAULT_CONFIG", "load_design", "load_design_library"]
