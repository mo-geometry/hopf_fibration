# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Versioning Scheme

- **Major (X)**: Breaking API changes, architectural rewrites
- **Minor (Y)**: New features, backwards-compatible additions
- **Patch (Z)**: Bug fixes, minor improvements, dependency updates

## [Unreleased]

## [0.1.0] - 2026-05-24

### Added
- Project scaffolding: src layout, pyproject.toml with hatchling build backend
- UV sphere mesh generator with per-vertex normals and configurable resolution
- Blinn-Phong vertex/fragment shaders with per-vertex colour and alpha transparency
- Line shaders for path and wireframe rendering
- Ground-plane grid with colour-coded axes (red X, blue Z)
- Mouse-driven orbit camera controller (drag to rotate, scroll to zoom)
- Geometric transforms: Euler-to-rotation, axis-angle (Rodrigues), look-at view matrix
- GLFW windowed application rendering a lit unit sphere on a grid
- GitHub Actions CI/CD: lint (ruff), typecheck (mypy strict), test matrix (3.11/3.12/3.13)
- Pre-commit hooks: ruff lint + format, mypy, trailing whitespace, YAML/TOML checks
- Makefile with dev, lint, format, typecheck, test, clean targets
- Test suite: sphere mesh (9 tests), orbit camera (8 tests), transforms (10 tests)
- Empty `geometry/` subpackage for parallel transport and Hopf fibration work

[Unreleased]: https://github.com/mo-geometry/hopf_fibration/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/mo-geometry/hopf_fibration/releases/tag/v0.1.0
