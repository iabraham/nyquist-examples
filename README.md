# nyquist-examples

Interactive Pluto notebook for teaching how **Bode magnitude/phase data maps directly to Nyquist curves**.

## What changed in this revision

- Replaced `ControlSystems.jl` with a lightweight in-notebook transfer-function evaluator (`RationalTF`) built from polynomial coefficients.
- Switched plotting from `Plots.jl`/`GR` to `CairoMakie` and `PolarAxis`.
- Added per-section viewport controls (radial and angular zoom) so students can focus on specific regions of the Nyquist polar plane.

## Notebook

- File: `nyquist_examples.jl`
- Topics covered:
  - Section 1: \(G_1(s)=1/((s+1)(s+2))\)
  - Section 2: \(G_2(s)=1/((s-1)(s^2+2s+3))\)
  - Section 3: \(G_3(s)=(s-1)/((s+2)(s^2-s+1))\)

Each section includes:

- Frequency slider for linked Bode/Nyquist marker movement.
- Polar zoom controls:
  - `rmax` for radial viewport.
  - `θzoom` for angular window (`[-θzoom, +θzoom]` degrees).

## Run locally

1. Install Julia 1.10+ (1.12.x recommended).
2. Start Pluto:

```julia
using Pkg
Pkg.add("Pluto")
using Pluto
Pluto.run()
```

3. Open `nyquist_examples.jl` in Pluto.
4. Pluto will resolve notebook dependencies from embedded project metadata.

## Dependencies (embedded in notebook)

- `CairoMakie`
- `PlutoUI`

No `ControlSystems.jl` dependency is required.
