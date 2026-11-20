"""A signed network average must not be mistaken for local minimality."""

from types import SimpleNamespace

import numpy as np

from crystalops.analysis.grain_boundary import _statistics


def _surface(values, areas):
    values = np.asarray(values, dtype=float)
    count = len(values)
    return SimpleNamespace(
        face_mean_curvature=values, face_area=np.asarray(areas, dtype=float),
        face_valid=np.ones(count, dtype=bool), face_junction=np.zeros(count, dtype=bool),
        face_projection_distance=np.zeros(count), faces=np.zeros((count, 3), dtype=int),
    )


def test_zero_signed_average_does_not_imply_zero_local_curvature():
    result = _statistics([_surface([-0.1, 0.1], [2.0, 2.0])])
    assert result["mean_signed_H"] == 0
    assert np.isclose(result["mean_absolute_H"], 0.1)
    assert np.isclose(result["rms_H"], 0.1)
    assert result["signed_cancellation_fraction"] == 1
    assert result["exact_zero_H_valid_face_count"] == 0
    assert result["positive_H_area_fraction"] == 0.5
    assert result["negative_H_area_fraction"] == 0.5


def test_true_zero_and_invalid_samples_remain_distinct():
    surface = _surface([0.0, np.nan], [1.0, 3.0])
    surface.face_valid[1] = False
    result = _statistics([surface])
    assert result["mean_absolute_H"] == 0
    assert result["rms_H"] == 0
    assert result["exact_zero_H_valid_face_count"] == 1
    assert result["signed_cancellation_fraction"] is None
    assert result["valid_area_fraction"] == 0.25


def test_flipping_one_interface_normal_changes_signed_mean_only():
    first = _statistics([_surface([0.2], [1]), _surface([-0.1], [2])])
    second = _statistics([_surface([-0.2], [1]), _surface([-0.1], [2])])
    assert first["mean_signed_H"] != second["mean_signed_H"]
    assert first["mean_absolute_H"] == second["mean_absolute_H"]
    assert first["rms_H"] == second["rms_H"]
