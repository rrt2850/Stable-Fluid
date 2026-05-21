# Stable Fluid Simulation Writeup

---

## Running the code

- Run using `start.sh` (Linux/macOS) or `start.ps1` (Windows).
  - Both scripts accept these arguments:
    - `gridWidth`, the width of the simulation
    - `gridHeight`, the height of the simulation
    -  `gridWidth` and `gridHeight` are the low-res resolutions for the current frame, all frames are upscaled to 
      - `FINAL_STILL_WIDTH = 2400` and `FINAL_STILL_HEIGHT = 1200`
    - `emitterCount`, the number of emitters to generate in emitter generation mode. Enter a random number here if you're running in the other mode.
    - `simulationSteps`, the number of steps to progress the simulation before it ends
    - `exportVideo`, 1 or 0, whether you want the video to be exported in addition to the final png
    - `logEveryStep`, 1 or 0, whether you want a console log marking every iteration. I generally just keep this on.
  - If omitted, defaults are supplied by the script itself. 
- The script starts `Main`, passing these six arguments.

---

## `Main.main()`: Setting up the simulation

Chronologically, `Main` does this:

- Parse command line args into a `SimulationConfig` object
  - Validates the inputs in `parseConfig`
- Construct `FluidGrid(width, height, cellSize)`
  - `cellSize = 1 / max(width, height)`
  - This normalizes the domain so derivatives scale with resolution
- Construct `SimulationParameters` from `Constants`:
  - `TIMESTEP`, how many seconds pass every step of the simulation
    - Smaller is more stable but slower. Larger is faster, but much less stable because advection relies on timestep.
  - `VISCOSITY`, how flowy the grid liquid is
    - Not accurate to real viscosity values
  - `DIFFUSION_RATE`, how fast fluid sources diffuse
    - For a more oily fluid, increase this
  - `SOLVER_ITERATIONS`, the number of Jacobi iterations to do per step (explained more later)
  - `VORTICITY_CONFINEMENT`, how much swirl to add back in after diffusion (explained more later)
- Do different things depending on which function is commented out:
  - If `InitializeRegionDivision` is enabled
    - Two walls are made with four `RadialFluidEmitters`, two on each side
    - Two more radial fluid emitters are made in the middle.
    - There is a gap in the wall that the `RadialFluidEmitter` fluid flows through
    - All the fluids meet in the middle as they're sucked into a `Vortex` in the middle
  - If `InitializeEdgeShooterDemo` is enabled
    - `EmitterMaker` is used to generate `emitterCount` edge emitters all pointing within +- `EMITTER_ANGLE_VARIATION_DEGREES` of the center of the grid.
- Open a window so you can view simulation frames real-time
- Step through the simulation `simulationSteps` times:
  - Calls `FluidSolver.step()` to progress the simulation
  - Saves low-res frames for the mp4 in a temporary directory to be compiled at the end of the simulation
  - Every `SNAPSHOT_INTERVAL` frames, it saves a high-res upscaled version of the current frame to the results directory
- After the simulation is complete, save the resulting image to Density.png
- Compiles all the frames into one mp4 if saving to video is enabled

---

## `FluidSolver` initialization

When `FluidSolver` is constructed, it allocates and initializes the following objects in the constructor in this order:

1. **Core references**
   - `FluidGrid grid` and `SimulationParameters parameters`
   - Purpose: keep a shared reference to domain geometry and all constants

2. **Velocity field object**
   - `velocityField = new VectorField(grid.totalCellCount)`
   - Creates 4 arrays (inside `VectorField`) sized to `grid.totalCellCount`:
     - `readVelocityX`, `writeVelocityX`, `readVelocityY`, `writeVelocityY`
   - Initial values: all `0.0f`
   - Purpose:
     - Read arrays hold the current step state
     - Write arrays are scratch/output buffers for diffusion/advection/projection
     - The solver swaps read/write roles after each substep (ping-pong buffering)

3. **RGB density field objects**
   - `redDensityField   = new ScalarField(grid.totalCellCount)`
   - `greenDensityField = new ScalarField(grid.totalCellCount)`
   - `blueDensityField  = new ScalarField(grid.totalCellCount)`
   - Each `ScalarField` contains double buffers (`readValues`, `writeValues`) of length
     `grid.totalCellCount`, initialized to `0.0f`
   - Purpose:
     - Store per-cell concentration for each color channel
     - Keeps independent buffers so each scalar solve/advect stage can write safely

4. **Pressure and divergence work fields**
   - `pressureField   = new ScalarField(grid.totalCellCount)`
   - `divergenceField = new ScalarField(grid.totalCellCount)`
   - Initial values: all zero in both read/write arrays.
   - Purpose:
     - `divergenceField` stores divergence values computed during projection
       - Divergence represents how much the velocity field is expanding or compressing. Divergence is negative the flow is converging to that point. if it's positive, the flow is diverging away from the point
       - For incompressible fluids, the goal is to get the divergence values to 0.
     - `pressureField` stores the Poisson solution used to subtract the gradient at a given point
       - By subtracting the gradient at a given point from the divergence field, we reduce divergence and enforce incompressibility by adjusting the divergence closer to 0
     - Persistently allocated once to avoid per-frame allocations // TODO: Figure out what this means

5. **Numerical helper objects**
   - `linearSolver = new LinearSolver(parameters.getLinearSolverIterations())`
     - Initialized with the Jacobi iteration count from `SimulationParameters`.
     - Purpose: iterative solve backend for velocity/density diffusion and pressure Poisson steps.
     - // TODO: reword this
   - `boundaryHandler = new BoundaryHandler()`
     - Purpose: apply edge boundary conditions after the linear solver runs to clamp vectors that are going towards boundaries
     - // TODO: Check if this is correct

6. **Wall occupancy map**
   - `wallMask = new boolean[grid.totalCellCount]`
   - Initial values: all `false` (no blocked cells yet).
   - Purpose: fast O(1) “is this cell a wall?” checks in collision and source/vortex logic.
   - Applies boundary rules to fluids interacting with wall cells
     - // TODO: Check if this is correct

7. **Configured injector/collider lists (null-safe copies)**
   - `redensitySources`, `emitters`, `radialEmitters`, `vortexes`, and `walls`
   - Initialization behavior:
     - If the incoming list is non-null, constructor copies it into a new `ArrayList<>(...)`.
     - If null, constructor uses `new ArrayList<>()` (empty list).
   - Purpose:
     - Prevent external mutation from affecting solver internals.
     - Guarantee all loops in `step()` can run without null checks.
   - // TODO: Check if we need to make a new copy or if we can just use the existing objects

8. **Input validation + wall rasterization**
   - Validates every source/emitter/vortex/wall endpoint is non-null and in bounds.
   - Calls `rebuildWallMask()` to rasterize all wall segments into `wallMask`.
     - // TODO: Look into specifics of how this works
   - Result:
     - Constructor fails fast for invalid setup.
     - Simulation loop can assume valid initial state.

---

## `FluidSolver.step()` in execution order

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

### 1 `addSources()`
- Each emitter adds it's density and velocity to the simulation cells using its own emitter-specific function
  - density (color mass)
  - velocity (momentum)
- Vortex objects apply custom influence fields to velocity/density
- Horizontal and vertical velocity arrays receive a force term each step
  - // TODO: Reword force term to be normal
- RGB density arrays (`readValues` in each `ScalarField`) receive injected density each step

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

## 10) End-to-end mental model (technical shorthand)

Each frame is a fractional-step integration of incompressible flow with scalar transport:

- **Inject/force** → **viscous diffuse** → **project** → **advect** → **swirl correction** → **project** → **scalar diffuse+advect** → **obstacle/boundary enforcement**.

This ordering is why the solver is both stable and visually expressive:

- Stability from implicit diffusion + semi-Lagrangian advection + repeated projection.
- Detail retention from vorticity confinement.
- Physical plausibility from divergence control and wall constraints.
