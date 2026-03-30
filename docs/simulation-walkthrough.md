# Stable Fluid Simulation Walkthrough (Chronological + Technical)

This document explains **how your simulation progresses in time**, mapping each phase to key functions and to the abstract math behind Stable Fluids.

---

## 1) Process startup: scripts, compilation, and runtime inputs

- You launch via `start.sh` (Linux/macOS) or `start.ps1` (Windows).
  - Each script accepts the same ordered arguments:
    - `gridWidth`
    - `gridHeight`
    - `emitterCount`
    - `simulationSteps`
    - `exportVideo`
    - `logEveryStep`
  - If omitted, defaults are supplied by the script itself.
- The script compiles `src/*.java` into `out/classes` using `javac`.
- The JVM starts `Main`, passing those six arguments.

### Why this matters technically
- Grid resolution controls numerical cost and smallest representable flow features.
- Number of steps controls total simulated time:
  - `simulatedTime ≈ simulationSteps × TIMESTEP`.
- Video/log flags only affect I/O; they do not change solver physics.

---

## 2) `Main.main(...)`: building the simulation state graph

Chronologically, `Main` does this:

- Parse CLI args into a typed `SimulationConfig`.
  - Validates positive integers and strict booleans.
- Construct `FluidGrid(width, height, cellSize)`.
  - `cellSize = 1 / max(width, height)`.
  - This normalizes the domain so derivatives scale with resolution.
- Construct `SimulationParameters` from `Constants`:
  - `TIMESTEP` (`Δt`)
  - `VISCOSITY` (`ν`)
  - `DIFFUSION_RATE` (`κ`)
  - `SOLVER_ITERATIONS`
  - `VORTICITY_CONFINEMENT` (`ε`)
- Create injectors and constraints for this specific setup:
  - Directional edge emitters (currently empty in this configuration).
  - Multiple `RadialFluidEmitter` objects for localized colored momentum + density injection.
  - A central `Vortex` influence object.
  - Four vertical wall segments with gaps, forming obstacle channels.
- Create `FluidSolver` with all of the above.

### Why `cellSize` is important
Most spatial operators are finite differences of the form:

- `∂f/∂x ≈ (f(x+1) - f(x-1)) / (2Δx)`

So `Δx = cellSize` appears in diffusion, divergence, pressure projection, and advection scaling.

---

## 3) `FluidSolver` initialization: field allocation and invariants

When `FluidSolver` is constructed, it creates these persistent fields:

- **Velocity field** (`VectorField`)
  - `u` = x-velocity
  - `v` = y-velocity
  - Stored with read/write buffers (ping-pong swapping)
- **Density fields** (`ScalarField`)
  - Separate channels: red, green, blue
  - Also double-buffered
- **Pressure and divergence** (`ScalarField` each)
  - Temporary but persistent work buffers reused each step

It also initializes:

- `LinearSolver` (iterative relaxation for implicit solves / Poisson)
- `BoundaryHandler` (domain-edge boundary conditions)
- `wallMask` (boolean occupancy map for interior obstacles)

Then it validates all sources/emitters/walls are in bounds and builds wall occupancy from geometric segments.

### Conceptual invariant
At the end of construction:

- All solver-owned arrays have consistent sizes (`grid.totalCellCount`).
- Inputs are safe (validated positions).
- Obstacles are pre-rasterized to avoid expensive geometry checks inside each numerical kernel.

---

## 4) Main simulation loop (`for step = 1..N`)

For each step in `Main`:

- Call `solver.step()` (core fluid update).
- Optionally render current fields to an image.
- Optionally save frame for MP4 generation.
- Optionally save high-resolution intermittent snapshots.

After the loop:

- Save final still image.
- Optionally encode MP4 using `ffmpeg`.

---

## 5) `FluidSolver.step()` in exact execution order (with math intuition)

Your implementation follows a classic Stable Fluids split-step pipeline:

1. `addSources()`
2. `applyWallCollisions()`
3. `diffuseVelocity()`
4. `projectVelocity()`
5. `advectVelocity()`
6. `applyVorticity()`
7. `projectVelocity()`
8. `applyWallCollisions()`
9. `diffuseDensity(red/green/blue)`
10. `advectDensity(red/green/blue)`
11. `applyWallCollisions()`

Below is what each does technically.

### 5.1 `addSources()` — explicit forcing/injection

- Point sources add scalar density directly to RGB channels.
- Emitters inject:
  - density (color mass)
  - velocity (momentum)
- Vortex objects apply custom influence fields to velocity/density.

Abstractly this corresponds to source terms in PDE form:

- `∂u/∂t = ... + f`
- `∂ρ/∂t = ... + s`

where `f` is momentum forcing and `s` is density production.

### 5.2 `applyWallCollisions()` — interior obstacle constraints

- For any wall cell:
  - zero velocity and scalar quantities.
- For neighboring fluid cells:
  - remove velocity components that point into walls.

This approximates a no-penetration condition at obstacle surfaces.

### 5.3 `diffuseVelocity()` — implicit viscosity solve

Physical model (continuous):

- `∂u/∂t = ν ∇²u`

Discrete implicit update solved iteratively:

- `(I - a∇²) u^{n+1} = u^n`, where `a = Δt·ν/Δx²`

Implementation notes:

- Uses `LinearSolver.solve(...)` for x and y components separately.
- Implicit solve is unconditionally stable for large enough `Δt` relative to explicit alternatives.

### 5.4 `projectVelocity()` — enforce incompressibility

Goal: make velocity divergence-free (`∇·u = 0`).

Pipeline in function:

- Compute divergence via centered differences.
- Solve pressure Poisson equation:
  - `∇²p = ∇·u` (scaled by `Δx²` in discrete form)
- Subtract pressure gradient:
  - `u ← u - ∇p`

Result:

- Removes local compression/expansion artifacts introduced by forcing/advection.
- Produces incompressible-looking flow.

### 5.5 `advectVelocity()` — semi-Lagrangian transport of momentum

For each cell center `x`:

- Backtrace:
  - `x_src = x - (Δt/Δx)·u(x)`
- Sample old velocity at `x_src` with bilinear interpolation.

Why semi-Lagrangian:

- Very stable even for larger time steps.
- Trade-off: introduces numerical dissipation (smoothing), later countered partly by vorticity confinement.

### 5.6 `applyVorticity()` — restore small-scale rotational detail

The function computes 2D scalar curl:

- `ω = ∂v/∂x - ∂u/∂y`

Then:

- Compute gradient of `|ω|` to find edges of vortical regions.
- Normalize gradient to `N`.
- Apply confinement force perpendicular to `N`, scaled by `ω`:
  - conceptually `F_conf ~ ε (N × ω)` in 2D form.

Effect:

- Re-amplifies visually important swirling detail that semi-Lagrangian advection tends to damp.

### 5.7 second `projectVelocity()` — re-divergence cleanup

- Vorticity force can reintroduce divergence.
- Running projection again restores near-incompressibility before density transport.

### 5.8 second `applyWallCollisions()` — post-velocity constraint pass

- Ensures any velocity moved/created into obstacle directions is removed before density steps.

### 5.9 `diffuseDensity(...)` — implicit scalar diffusion per channel

For each channel `ρ ∈ {R,G,B}`:

- Solve `(I - a∇²) ρ^{n+1} = ρ^n`, with `a = Δt·κ/Δx²`.

Interpretation:

- Mimics molecular/turbulent mixing at unresolved scales.
- Larger `κ` causes faster plume broadening.

### 5.10 `advectDensity(...)` — semi-Lagrangian scalar transport

For each channel and cell:

- Backtrace with current velocity.
- Bilinear sample old scalar at traced position.

This applies the advection equation:

- `∂ρ/∂t + u·∇ρ = 0` (plus sources/diffusion handled in other substeps).

### 5.11 final `applyWallCollisions()` — scalar/velocity cleanup near solids

- Removes scalar/velocity contamination inside walls after transport.
- Re-applies edge and obstacle consistency constraints.

---

## 6) Boundary handling and ghost-cell behavior (abstract view)

`BoundaryHandler` applies boundary-type-specific rules:

- Scalar fields: copied/adjusted to avoid invalid edge access and maintain stable stencil behavior.
- Velocity components: sign/axis-aware handling so normal components behave correctly at boundaries.

At a high level, this is equivalent to ghost-cell conditions used in finite-difference PDE solvers.

---

## 7) Linear solver role in the pipeline

`LinearSolver` is reused in two major places:

- Implicit diffusion solves
- Pressure Poisson solve for projection

Why iterative relaxation is used:

- Sparse grid Laplacian systems are large but structured.
- Iterative methods are straightforward, memory-efficient, and adequate for graphics-oriented fluid solves.
- `SOLVER_ITERATIONS` controls quality/cost tradeoff.

Practical interpretation:

- Too few iterations: softer/inexact pressure solve, more residual divergence.
- More iterations: cleaner incompressibility and diffusion fidelity, higher CPU time.

---

## 8) Rendering pipeline (`ImageRenderer`) and what it visualizes

The renderer visualizes **density**, not velocity vectors directly.

- It reads current RGB density fields.
- Finds per-channel maxima over the grid.
- Normalizes each channel independently.
- Bilinearly samples onto output pixel grid.
- Writes ARGB pixels.

Why per-channel normalization:

- Prevents a dominant channel from washing out weaker channels.
- Improves color separation for qualitative analysis.

Output modes:

- Final upscaled still image.
- Optional per-step frames and MP4 assembly.

---

## 9) Parameter sensitivity: what each knob does numerically

- `TIMESTEP (Δt)`
  - Larger: faster apparent evolution; can increase transport error and visual instability.
  - Smaller: more stable/accurate temporal integration; more steps needed.
- `VISCOSITY (ν)`
  - Increases momentum diffusion; suppresses fine velocity features.
- `DIFFUSION_RATE (κ)`
  - Increases scalar smoothing/mixing rate.
- `SOLVER_ITERATIONS`
  - Improves convergence of implicit/Poisson solves.
- `VORTICITY_CONFINEMENT (ε)`
  - Adds rotational energy at small scales for crisper curls.

---

## 10) End-to-end mental model (technical shorthand)

Each frame is a fractional-step integration of incompressible flow with scalar transport:

- **Inject/force** → **viscous diffuse** → **project** → **advect** → **swirl correction** → **project** → **scalar diffuse+advect** → **obstacle/boundary enforcement**.

This ordering is why the solver is both stable and visually expressive:

- Stability from implicit diffusion + semi-Lagrangian advection + repeated projection.
- Detail retention from vorticity confinement.
- Physical plausibility from divergence control and wall constraints.
