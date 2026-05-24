"""Geometric transforms — rotation representations, coordinate conversions."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt


def euler_to_rotation(roll: float, pitch: float, yaw: float) -> npt.NDArray[np.float64]:
    """Convert Euler angles (ZYX convention, radians) to a 3x3 rotation matrix.

    Args:
        roll: Rotation about the x-axis.
        pitch: Rotation about the y-axis.
        yaw: Rotation about the z-axis.

    Returns:
        A 3x3 rotation matrix R = Rz(yaw) @ Ry(pitch) @ Rx(roll).
    """
    cr, sr = np.cos(roll), np.sin(roll)
    cp, sp = np.cos(pitch), np.sin(pitch)
    cy, sy = np.cos(yaw), np.sin(yaw)

    return np.array(
        [
            [cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr],
            [sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr],
            [-sp, cp * sr, cp * cr],
        ],
        dtype=np.float64,
    )


def axis_angle_to_rotation(axis: npt.NDArray[np.float64], angle: float) -> npt.NDArray[np.float64]:
    """Convert an axis-angle representation to a 3x3 rotation matrix.

    Uses Rodrigues' rotation formula.

    Args:
        axis: Unit vector representing the rotation axis (3,).
        angle: Rotation angle in radians.

    Returns:
        A 3x3 rotation matrix.
    """
    axis = axis / np.linalg.norm(axis)
    k = np.array(
        [
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0],
        ],
        dtype=np.float64,
    )
    result: npt.NDArray[np.float64] = (
        np.eye(3, dtype=np.float64) + np.sin(angle) * k + (1 - np.cos(angle)) * (k @ k)
    )
    return result


def look_at(
    eye: npt.NDArray[np.float64],
    target: npt.NDArray[np.float64],
    up: npt.NDArray[np.float64] | None = None,
) -> npt.NDArray[np.float64]:
    """Compute a 4x4 view matrix looking from *eye* toward *target*.

    Args:
        eye: Camera position in world space (3,).
        target: Point the camera looks at (3,).
        up: World-space up vector; defaults to [0, 1, 0].

    Returns:
        A 4x4 view matrix.
    """
    if up is None:
        up = np.array([0.0, 1.0, 0.0], dtype=np.float64)

    forward = target - eye
    forward = forward / np.linalg.norm(forward)

    right = np.cross(forward, up)
    right = right / np.linalg.norm(right)

    true_up = np.cross(right, forward)

    mat = np.eye(4, dtype=np.float64)
    mat[0, :3] = right
    mat[1, :3] = true_up
    mat[2, :3] = -forward
    mat[:3, 3] = -mat[:3, :3] @ eye
    return mat
