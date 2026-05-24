"""Unit sphere mesh generation via UV parameterisation."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt


def create_uv_sphere(
    radius: float = 1.0,
    slices: int = 48,
    stacks: int = 24,
    color: tuple[float, float, float] = (0.4, 0.6, 0.85),
) -> tuple[
    npt.NDArray[np.float32],
    npt.NDArray[np.float32],
    npt.NDArray[np.float32],
    npt.NDArray[np.uint32],
]:
    """Generate a UV sphere mesh with per-vertex normals and colours.

    On a unit sphere the normal at each vertex equals the position,
    which is handy for differential-geometry visualisations.

    Args:
        radius: Sphere radius.
        slices: Number of longitudinal divisions.
        stacks: Number of latitudinal divisions.
        color: RGB colour for the sphere surface.

    Returns:
        Tuple of (positions, normals, colors, indices) where:
        - positions: (N, 3) float32 vertex positions
        - normals: (N, 3) float32 per-vertex normals (= normalised position)
        - colors: (N, 3) float32 per-vertex RGB colours
        - indices: (M,) uint32 triangle indices
    """
    positions = []
    normals = []

    for i in range(stacks + 1):
        phi = np.pi * i / stacks  # 0 at north pole, pi at south
        for j in range(slices + 1):
            theta = 2.0 * np.pi * j / slices

            x = radius * np.sin(phi) * np.cos(theta)
            y = radius * np.cos(phi)
            z = radius * np.sin(phi) * np.sin(theta)

            positions.append([x, y, z])
            # Normal = position / radius (unit sphere: same as position)
            normals.append([x / radius, y / radius, z / radius])

    pos_array = np.array(positions, dtype=np.float32)
    norm_array = np.array(normals, dtype=np.float32)
    color_array = np.full_like(pos_array, color, dtype=np.float32)

    # Build triangle indices
    indices = []
    for i in range(stacks):
        for j in range(slices):
            first = i * (slices + 1) + j
            second = first + slices + 1

            indices.extend([first, second, first + 1])
            indices.extend([second, second + 1, first + 1])

    idx_array = np.array(indices, dtype=np.uint32)

    return pos_array, norm_array, color_array, idx_array
