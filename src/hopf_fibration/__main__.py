"""Hopf Fibration — GLFW windowed renderer.

Renders a unit sphere with Blinn-Phong shading, a ground-plane grid,
and mouse-driven orbit camera controls.

Run with:
    python -m hopf_fibration
"""

from __future__ import annotations

import sys

import glfw  # type: ignore[import-untyped]
import moderngl
import numpy as np

from hopf_fibration.assets.grid import create_grid
from hopf_fibration.assets.sphere import create_uv_sphere
from hopf_fibration.core.orbit_camera import OrbitCamera
from hopf_fibration.rendering.shaders import (
    FRAGMENT_SHADER,
    LINE_FRAGMENT_SHADER,
    LINE_VERTEX_SHADER,
    VERTEX_SHADER,
)

# Window dimensions
WIDTH = 1280
HEIGHT = 720
TITLE = "Hopf Fibration — ModernGL"


def _perspective(fov_y: float, aspect: float, near: float, far: float) -> np.ndarray:  # type: ignore[type-arg]
    """Build a 4x4 perspective projection matrix (row-major)."""
    f = 1.0 / np.tan(fov_y / 2.0)
    proj = np.zeros((4, 4), dtype=np.float32)
    proj[0, 0] = f / aspect
    proj[1, 1] = f
    proj[2, 2] = (far + near) / (near - far)
    proj[3, 2] = -1.0
    proj[2, 3] = (2.0 * far * near) / (near - far)
    return proj


def main() -> None:
    """Launch the windowed renderer."""
    # --- GLFW initialisation ---
    if not glfw.init():
        sys.exit("Failed to initialise GLFW")

    glfw.window_hint(glfw.CONTEXT_VERSION_MAJOR, 3)
    glfw.window_hint(glfw.CONTEXT_VERSION_MINOR, 3)
    glfw.window_hint(glfw.OPENGL_PROFILE, glfw.OPENGL_CORE_PROFILE)
    glfw.window_hint(glfw.OPENGL_FORWARD_COMPAT, glfw.TRUE)
    glfw.window_hint(glfw.SAMPLES, 4)

    window = glfw.create_window(WIDTH, HEIGHT, TITLE, None, None)
    if not window:
        glfw.terminate()
        sys.exit("Failed to create GLFW window")

    glfw.make_context_current(window)
    glfw.swap_interval(1)  # Vsync

    # --- ModernGL context ---
    ctx = moderngl.create_context()
    ctx.enable(moderngl.DEPTH_TEST)
    ctx.enable(moderngl.BLEND)
    ctx.blend_func = moderngl.SRC_ALPHA, moderngl.ONE_MINUS_SRC_ALPHA

    # Query actual framebuffer size (Retina displays are 2x)
    fb_width, fb_height = glfw.get_framebuffer_size(window)
    ctx.viewport = (0, 0, fb_width, fb_height)

    print(f"OpenGL: {ctx.info['GL_VERSION']}")
    print(f"Renderer: {ctx.info['GL_RENDERER']}")
    print(f"Framebuffer: {fb_width}x{fb_height}")

    # --- Compile shaders ---
    prog = ctx.program(vertex_shader=VERTEX_SHADER, fragment_shader=FRAGMENT_SHADER)
    line_prog = ctx.program(
        vertex_shader=LINE_VERTEX_SHADER, fragment_shader=LINE_FRAGMENT_SHADER
    )

    # --- Create sphere mesh ---
    positions, normals, colors, indices = create_uv_sphere(
        radius=1.0, slices=48, stacks=24, color=(0.4, 0.6, 0.85)
    )

    vbo_pos = ctx.buffer(positions.tobytes())
    vbo_norm = ctx.buffer(normals.tobytes())
    vbo_col = ctx.buffer(colors.tobytes())
    ibo = ctx.buffer(indices.tobytes())

    vao = ctx.vertex_array(
        prog,
        [
            (vbo_pos, "3f", "in_position"),
            (vbo_norm, "3f", "in_normal"),
            (vbo_col, "3f", "in_color"),
        ],
        index_buffer=ibo,
        index_element_size=4,
    )

    # --- Create ground grid ---
    grid_positions, grid_colors = create_grid(size=5.0, divisions=10)
    grid_vbo_pos = ctx.buffer(grid_positions.tobytes())
    grid_vbo_col = ctx.buffer(grid_colors.tobytes())

    grid_vao = ctx.vertex_array(
        line_prog,
        [
            (grid_vbo_pos, "3f", "in_position"),
            (grid_vbo_col, "3f", "in_color"),
        ],
    )

    # --- Camera setup ---
    camera = OrbitCamera(distance=4.0, azimuth=0.5, elevation=0.4)

    # Projection matrix
    aspect = fb_width / fb_height
    fov_y = float(np.radians(60.0))
    near, far = 0.1, 100.0
    proj = _perspective(fov_y, aspect, near, far)

    # Model matrix (identity — sphere at origin)
    model = np.eye(4, dtype=np.float32)

    # --- GLFW callbacks ---
    def _mouse_button_cb(
        _win: object, button: int, action: int, mods: int
    ) -> None:
        camera.on_mouse_button(button, action, mods)

    def _cursor_pos_cb(_win: object, x: float, y: float) -> None:
        camera.on_cursor_pos(x, y)

    def _scroll_cb(_win: object, x_offset: float, y_offset: float) -> None:
        camera.on_scroll(x_offset, y_offset)

    def _framebuffer_size_cb(_win: object, width: int, height: int) -> None:
        nonlocal proj, aspect
        if height == 0:
            return
        ctx.viewport = (0, 0, width, height)
        aspect = width / height
        proj = _perspective(fov_y, aspect, near, far)

    glfw.set_mouse_button_callback(window, _mouse_button_cb)
    glfw.set_cursor_pos_callback(window, _cursor_pos_cb)
    glfw.set_scroll_callback(window, _scroll_cb)
    glfw.set_framebuffer_size_callback(window, _framebuffer_size_cb)

    # --- Light setup ---
    light_pos = np.array([5.0, 8.0, 5.0], dtype=np.float32)
    light_color = np.array([1.0, 1.0, 1.0], dtype=np.float32)

    # --- Render loop ---
    while not glfw.window_should_close(window):
        glfw.poll_events()

        if glfw.get_key(window, glfw.KEY_ESCAPE) == glfw.PRESS:
            glfw.set_window_should_close(window, True)

        ctx.clear(0.10, 0.10, 0.12, 1.0)

        # Build column-major matrices for OpenGL
        view = camera.view_matrix.astype(np.float32).T.copy()
        eye_pos = camera.eye.astype(np.float32)
        model_cm = model.T.copy()
        proj_cm = proj.T.copy()

        # --- Draw grid (behind sphere) ---
        line_prog["view"].write(view.tobytes())  # type: ignore[index]
        line_prog["projection"].write(proj_cm.tobytes())  # type: ignore[index]
        grid_vao.render(moderngl.LINES)

        # --- Draw sphere ---
        prog["model"].write(model_cm.tobytes())  # type: ignore[index]
        prog["view"].write(view.tobytes())  # type: ignore[index]
        prog["projection"].write(proj_cm.tobytes())  # type: ignore[index]
        prog["light_pos"].write(light_pos.tobytes())  # type: ignore[index]
        prog["light_color"].write(light_color.tobytes())  # type: ignore[index]
        prog["view_pos"].write(eye_pos.tobytes())  # type: ignore[index]
        prog["alpha"].value = 0.85  # type: ignore[index]

        # Draw back faces first, then front faces for correct transparency
        ctx.cull_face = "front"
        vao.render(moderngl.TRIANGLES)
        ctx.cull_face = "back"
        vao.render(moderngl.TRIANGLES)

        glfw.swap_buffers(window)

    # --- Cleanup ---
    vao.release()
    vbo_pos.release()
    vbo_norm.release()
    vbo_col.release()
    ibo.release()
    grid_vao.release()
    grid_vbo_pos.release()
    grid_vbo_col.release()
    prog.release()
    line_prog.release()
    ctx.release()
    glfw.terminate()


if __name__ == "__main__":
    main()
