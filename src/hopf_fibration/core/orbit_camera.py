"""Mouse-driven orbit camera controller.

Provides arcball-style orbit controls: click-drag to rotate around a
target point, scroll to zoom in/out.
"""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt

from hopf_fibration.core.transforms import look_at


class OrbitCamera:
    """Orbit camera that rotates around a target point.

    Attributes:
        target: The point the camera orbits around.
        distance: Distance from the target.
        azimuth: Horizontal angle in radians (around Y axis).
        elevation: Vertical angle in radians (from XZ plane).
        min_distance: Minimum zoom distance.
        max_distance: Maximum zoom distance.
    """

    def __init__(
        self,
        target: npt.NDArray[np.float64] | None = None,
        distance: float = 5.0,
        azimuth: float = 0.0,
        elevation: float = 0.3,
    ) -> None:
        self.target = target if target is not None else np.zeros(3, dtype=np.float64)
        self.distance = distance
        self.azimuth = azimuth
        self.elevation = elevation
        self.min_distance = 1.0
        self.max_distance = 50.0

        # Mouse state for drag tracking
        self._dragging = False
        self._last_x = 0.0
        self._last_y = 0.0
        self._sensitivity = 0.005

    @property
    def eye(self) -> npt.NDArray[np.float64]:
        """Compute camera position from spherical coordinates."""
        x = self.distance * math.cos(self.elevation) * math.sin(self.azimuth)
        y = self.distance * math.sin(self.elevation)
        z = self.distance * math.cos(self.elevation) * math.cos(self.azimuth)
        return self.target + np.array([x, y, z], dtype=np.float64)

    @property
    def view_matrix(self) -> npt.NDArray[np.float64]:
        """Compute the 4x4 view matrix for the current orbit state."""
        return look_at(self.eye, self.target)

    def on_mouse_button(self, button: int, action: int, _mods: int) -> None:
        """Handle mouse button events (GLFW callback signature).

        Args:
            button: GLFW mouse button code (0 = left).
            action: GLFW action code (1 = press, 0 = release).
            _mods: Modifier keys (unused).
        """
        if button == 0:  # Left mouse button
            self._dragging = action == 1

    def on_cursor_pos(self, x: float, y: float) -> None:
        """Handle cursor movement events (GLFW callback signature).

        Args:
            x: Cursor x position in window coordinates.
            y: Cursor y position in window coordinates.
        """
        if self._dragging:
            dx = x - self._last_x
            dy = y - self._last_y
            self.azimuth -= dx * self._sensitivity
            self.elevation += dy * self._sensitivity
            # Clamp elevation to avoid gimbal lock at poles
            self.elevation = max(-math.pi / 2 + 0.01, min(math.pi / 2 - 0.01, self.elevation))

        self._last_x = x
        self._last_y = y

    def on_scroll(self, _x_offset: float, y_offset: float) -> None:
        """Handle scroll events for zoom (GLFW callback signature).

        Args:
            _x_offset: Horizontal scroll offset (unused).
            y_offset: Vertical scroll offset (positive = zoom in).
        """
        self.distance -= y_offset * 0.5
        self.distance = max(self.min_distance, min(self.max_distance, self.distance))
