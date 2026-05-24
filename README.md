# Hopf Fibration

Interactive visualisation of Hopf fibrations, parallel transport, and
fibre bundles on the unit sphere, built with ModernGL and GLFW.

## Features (planned)

- Unit sphere rendering with wireframe and solid modes
- Path generation on S2 (great circles, loxodromes, arbitrary curves)
- Parallel transport frame animation along paths
- Hopf fibration: S3 to S2 fibre bundle visualisation
- Interactive orbit camera controls

## Quick Start

```bash
python -m venv .venv
source .venv/bin/activate
make dev
make test
python -m hopf_fibration
```

## References

- Hopf, H. (1931). "Uber die Abbildungen der dreidimensionalen Sphare auf die Kugelflache"
- Geometrically: the Hopf map sends each point of S3 to a point on S2,
  with each fibre being a great circle in S3.
