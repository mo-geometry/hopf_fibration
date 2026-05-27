"""Tests for UV sphere mesh generation."""

from __future__ import annotations

import numpy as np

from hopf_fibration.assets.sphere import create_uv_sphere


class TestCreateUvSphere:
    """Unit sphere mesh generation tests."""

    def test_vertex_count(self) -> None:
        positions, normals, colors, _indices = create_uv_sphere(slices=8, stacks=4)
        expected_verts = (4 + 1) * (8 + 1)  # (stacks+1) * (slices+1)
        assert positions.shape == (expected_verts, 3)
        assert normals.shape == (expected_verts, 3)
        assert colors.shape == (expected_verts, 3)

    def test_index_count(self) -> None:
        _, _, _, indices = create_uv_sphere(slices=8, stacks=4)
        expected_tris = 8 * 4 * 2  # slices * stacks * 2 triangles
        assert len(indices) == expected_tris * 3

    def test_normals_are_unit_length(self) -> None:
        _, normals, _, _ = create_uv_sphere()
        lengths = np.linalg.norm(normals, axis=1)
        np.testing.assert_allclose(lengths, 1.0, atol=1e-5)

    def test_vertices_on_sphere_surface(self) -> None:
        positions, _, _, _ = create_uv_sphere(radius=1.0)
        distances = np.linalg.norm(positions, axis=1)
        np.testing.assert_allclose(distances, 1.0, atol=1e-5)

    def test_custom_radius(self) -> None:
        positions, _, _, _ = create_uv_sphere(radius=2.5)
        distances = np.linalg.norm(positions, axis=1)
        np.testing.assert_allclose(distances, 2.5, atol=1e-5)

    def test_normals_equal_normalised_positions(self) -> None:
        positions, normals, _, _ = create_uv_sphere(radius=1.0)
        np.testing.assert_allclose(normals, positions, atol=1e-5)

    def test_indices_in_range(self) -> None:
        positions, _, _, indices = create_uv_sphere()
        assert indices.max() < len(positions)
        assert indices.min() >= 0

    def test_dtypes(self) -> None:
        positions, normals, colors, indices = create_uv_sphere()
        assert positions.dtype == np.float32
        assert normals.dtype == np.float32
        assert colors.dtype == np.float32
        assert indices.dtype == np.uint32

    def test_custom_color(self) -> None:
        _, _, colors, _ = create_uv_sphere(color=(1.0, 0.0, 0.0))
        np.testing.assert_allclose(colors[:, 0], 1.0)
        np.testing.assert_allclose(colors[:, 1], 0.0)
        np.testing.assert_allclose(colors[:, 2], 0.0)
