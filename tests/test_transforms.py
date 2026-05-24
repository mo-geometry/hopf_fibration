"""Tests for geometric transforms — rotations, look-at, coordinate conversion."""

from __future__ import annotations

import numpy as np
import pytest

from hopf_fibration.core.transforms import (
    axis_angle_to_rotation,
    euler_to_rotation,
    look_at,
)


class TestEulerToRotation:
    """Euler-angle to rotation matrix conversion."""

    def test_identity_for_zero_angles(self) -> None:
        r = euler_to_rotation(0.0, 0.0, 0.0)
        np.testing.assert_array_almost_equal(r, np.eye(3))

    def test_output_is_orthogonal(self) -> None:
        r = euler_to_rotation(0.3, 0.5, 0.7)
        np.testing.assert_array_almost_equal(r @ r.T, np.eye(3), decimal=12)

    def test_determinant_is_one(self) -> None:
        r = euler_to_rotation(1.0, -0.5, 2.0)
        assert np.linalg.det(r) == pytest.approx(1.0, abs=1e-12)

    def test_90_degree_yaw(self) -> None:
        r = euler_to_rotation(0.0, 0.0, np.pi / 2)
        rotated = r @ np.array([1.0, 0.0, 0.0])
        np.testing.assert_array_almost_equal(rotated, [0.0, 1.0, 0.0])


class TestAxisAngleToRotation:
    """Axis-angle (Rodrigues) to rotation matrix conversion."""

    def test_zero_angle_gives_identity(self) -> None:
        axis = np.array([0.0, 0.0, 1.0])
        r = axis_angle_to_rotation(axis, 0.0)
        np.testing.assert_array_almost_equal(r, np.eye(3))

    def test_180_about_z(self) -> None:
        axis = np.array([0.0, 0.0, 1.0])
        r = axis_angle_to_rotation(axis, np.pi)
        np.testing.assert_array_almost_equal(r @ np.array([1, 0, 0]), [-1, 0, 0])

    def test_orthogonality(self) -> None:
        axis = np.array([1.0, 1.0, 1.0])
        r = axis_angle_to_rotation(axis, 1.23)
        np.testing.assert_array_almost_equal(r @ r.T, np.eye(3), decimal=12)


class TestLookAt:
    """Look-at view matrix construction."""

    def test_output_shape(self) -> None:
        eye = np.array([0.0, 0.0, 5.0])
        target = np.array([0.0, 0.0, 0.0])
        mat = look_at(eye, target)
        assert mat.shape == (4, 4)

    def test_identity_when_aligned(self) -> None:
        eye = np.array([0.0, 0.0, 0.0])
        target = np.array([0.0, 0.0, -1.0])
        mat = look_at(eye, target)
        np.testing.assert_array_almost_equal(mat[:3, :3], np.eye(3))

    def test_translation_component(self) -> None:
        eye = np.array([3.0, 4.0, 5.0])
        target = np.array([0.0, 0.0, 0.0])
        mat = look_at(eye, target)
        eye_h = np.array([3.0, 4.0, 5.0, 1.0])
        transformed = mat @ eye_h
        np.testing.assert_array_almost_equal(transformed[:3], [0, 0, 0], decimal=10)
