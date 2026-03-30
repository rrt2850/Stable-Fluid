# Stable Fluid Simulation Walkthrough (Chronological + Technical)

This document explains **how the simulation progresses in time**, maps each phase to key functions/classes, and describes both:

- the **mathematical model** (for technical readers), and
- the **plain-English behavior** (for readers who just want to know what the code is doing in practical terms).

---

## 1) Process startup: scripts, compilation, and runtime inputs

You launch the simulation with either:

- `start.sh` (Linux/macOS), or
- `start.ps1` (Windows).

Both scripts accept the same ordered arguments:

1. `gridWidth`
2. `gridHeight`
3. `emitterCount`
4. `simulationSteps`
5. `exportVideo`
6. `logEveryStep`

If arguments are missing, each script falls back to default values.

Then the startup flow is:

1. Compile `src/*.java` into `out/classes` using `javac`.
2. Run `Main` with those six arguments.

### Technical meaning

- Grid size (`gridWidth`, `gridHeight`) controls numerical resolution, cost, and smallest flow detail that can be represented.
- Total simulated time is approximately:
  - `simulatedTime ≈ simulationSteps × TIMESTEP`.
- `exportVideo` and `logEveryStep` only change output behavior (I/O), not physics.

### Plain-English meaning

- Bigger grid = sharper details, but slower runtime.
- More steps = longer simulation.
- Video/log flags only affect what gets saved or printed; they do not change how the fluid moves.

---

## 2) `Main.main(...)`: building the simulation setup

`Main` performs these steps in order:

1. Parse CLI arguments into a typed configuration and validate them.
2. Build a `FluidGrid(width, height, cellSize)` where:
   - `cellSize = 1 / max(width, height)`.
3. Build `SimulationParameters` from constants:
   - `TIMESTEP (Δt)`
   - `VISCOSITY (ν)`
   - `DIFFUSION_RATE (κ)`
   - `SOLVER_ITERATIONS`
   - `VORTICITY_CONFINEMENT (ε)`
4. Build scene content:
   - emitters/sources,
   - a central vortex,
   - vertical walls with gaps.
5. Construct `FluidSolver` with grid, parameters, sources, and walls.

### Technical meaning

- `cellSize` appears in finite-difference derivatives and scaling terms throughout diffusion, projection, and advection.
- Normalizing domain scale this way helps behavior stay consistent across resolutions.

### Plain-English meaning

- `Main` is the “wiring” phase: it decides the size of the world, how fast physics runs, where fluid gets added, and where obstacles exist.
- `cellSize` is the physical size of one grid square; many equations use it to decide how strongly neighboring cells interact.

---

## 3) `FluidSolver` initialization: what data is stored

When `FluidSolver` is created, it allocates and keeps these fields:

- **Velocity (`VectorField`)**
  - `u`: horizontal speed,
  - `v`: vertical speed,
  - with read/write buffers (ping-pong swap each substep).
- **Density (`ScalarField`)**
  - separate red/green/blue channels,
  - also double-buffered.
- **Pressure and divergence (`ScalarField`)**
  - reusable work buffers for incompressibility projection.

It also creates:

- `LinearSolver` for iterative solves,
- `BoundaryHandler` for domain edges,
- `wallMask` occupancy map for internal obstacles.

Then it validates source/wall positions and rasterizes walls into the mask.

### Technical meaning

- Pre-allocation avoids per-step memory churn.
- Rasterized wall masks avoid repeated geometry tests during kernels.

### Plain-English meaning

- Before time starts, the solver builds all arrays it will need and marks which cells are “solid wall” vs “fluid.”
- This front-loads setup so each simulation step can run quickly.

---

## 4) Main simulation loop (`for step = 1..N`)

Every loop iteration does:

1. `solver.step()` (the fluid update).
2. Optional rendering for preview/snapshots.
3. Optional frame saving for video generation.

After the loop:

- Save final image.
- Optionally invoke `ffmpeg` to create MP4.

### Plain-English meaning

- The simulation advances one small time slice at a time.
- After each slice, you can visualize it, and optionally store that frame.

---

## 5) `FluidSolver.step()` exact order and what each part does

Pipeline order:

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

This is a Stable Fluids split-step integration: each stage handles one physical effect while preserving numerical stability.

### 5.1 `addSources()` — inject mass and momentum

**Technical:** Adds source terms to velocity/density equations:

- `∂u/∂t = ... + f`
- `∂ρ/∂t = ... + s`

where `f` is forcing and `s` is scalar generation.

**Plain English:** This is where “new smoke/color and push” enter the scene from emitters/vortices.

### 5.2 `applyWallCollisions()` — enforce interior obstacles

**Technical:**

- Zero velocity/scalar in wall cells.
- Remove velocity components that point into walls near boundaries.

Approximates no-penetration constraints.

**Plain English:** Fluid cannot live inside walls, and nearby flow is redirected so it doesn’t tunnel through solid cells.

### 5.3 `diffuseVelocity()` — viscosity step (implicit)

**Technical:** Solves

- `(I - a∇²) u^{n+1} = u^n`, with `a = Δt·ν/Δx²`.

Done per velocity component via iterative linear solve.

**Plain English:** Velocity gets “smoothed” by neighbor mixing. High viscosity makes motion feel thicker and less sharp.

### 5.4 `projectVelocity()` — enforce incompressibility

**Technical:**

1. Compute divergence `∇·u`.
2. Solve pressure Poisson equation `∇²p = ∇·u`.
3. Subtract gradient `u ← u - ∇p`.

Goal: approximately divergence-free velocity.

**Plain English:** This removes fake compression/expansion so flow looks like liquid/smoke volume is conserved instead of magically appearing/disappearing.

### 5.5 `advectVelocity()` — move velocity with flow (semi-Lagrangian)

**Technical:** Backtrace each cell center:

- `x_src = x - (Δt/Δx)·u(x)`

Then bilinearly sample previous velocity at `x_src`.

**Plain English:** Follow where each parcel came from last step, and copy that motion forward. This is very stable but slightly blurs detail.

### 5.6 `applyVorticity()` — add swirl detail back

**Technical:**

1. Compute curl `ω = ∂v/∂x - ∂u/∂y`.
2. Compute direction from `∇|ω|`.
3. Apply confinement force scaled by `ε` and `ω`.

**Plain English:** Re-energizes small eddies/curls lost to numerical smoothing, making motion look more lively and less “mushy.”

### 5.7 second `projectVelocity()` — clean up divergence again

**Technical:** Vorticity confinement can introduce divergence, so projection runs a second time.

**Plain English:** After adding swirl, we re-balance flow so it still obeys incompressible behavior.

### 5.8 second `applyWallCollisions()` — re-apply obstacle constraints

**Technical:** Removes any newly introduced into-wall components.

**Plain English:** A safety pass that keeps walls respected before density transport.

### 5.9 `diffuseDensity(...)` — scalar mixing (implicit)

For each color channel `ρ ∈ {R,G,B}` solve:

- `(I - a∇²) ρ^{n+1} = ρ^n`, with `a = Δt·κ/Δx²`.

**Plain English:** Color spreads out over time. Higher diffusion makes plumes broaden faster.

### 5.10 `advectDensity(...)` — move color with velocity

**Technical:** Same semi-Lagrangian backtrace/sampling, but applied to scalar channels.

**Plain English:** The color follows the fluid motion.

### 5.11 final `applyWallCollisions()` — final consistency pass

**Technical:** Clears any transported values that ended up in wall cells and reasserts boundary consistency.

**Plain English:** Final cleanup: no fluid or color remains inside solids.

---

## 6) Domain boundaries and ghost-cell behavior

`BoundaryHandler` handles outer edges of the simulation domain.

### Technical meaning

- Scalar and vector fields receive edge-specific updates compatible with finite-difference stencils.
- Velocity edge handling is axis-aware (normal/tangential behavior).

### Plain-English meaning

- Edge rules decide what fluid does at the borders of the grid (bounce/slide/stop behavior depending on component and boundary type).
- Without this, edge cells would behave erratically and destabilize the simulation.

---

## 7) `LinearSolver`: why it exists and where it is used

Used in:

- implicit diffusion,
- pressure Poisson projection.

### Technical meaning

- Grid Laplacian systems are sparse and structured.
- Iterative relaxation is memory-light and good for graphics simulation.
- `SOLVER_ITERATIONS` trades speed for convergence quality.

### Plain-English meaning

- The solver repeatedly “nudges” values toward the right answer.
- More iterations = cleaner physics but slower frames.

---

## 8) Rendering (`ImageRenderer`): what you are actually seeing

Renderer visualizes **density**, not arrows of velocity.

Steps:

1. Read RGB density grids.
2. Compute per-channel maxima.
3. Normalize each channel independently.
4. Bilinear sample onto output pixel resolution.
5. Write ARGB image.

### Technical meaning

- Per-channel normalization improves dynamic range per color.
- Bilinear sampling reduces blockiness when scaling.

### Plain-English meaning

- You are seeing “where colored fluid exists,” not direct speed vectors.
- Normalization helps weaker colors remain visible instead of being drowned out.

---

## 9) Parameter guide: math and intuition together

- `TIMESTEP (Δt)`
  - **Technical:** temporal integration size.
  - **Plain English:** how much simulated time each step advances.
- `VISCOSITY (ν)`
  - **Technical:** momentum diffusion coefficient.
  - **Plain English:** higher values make motion thicker and less sharp.
- `DIFFUSION_RATE (κ)`
  - **Technical:** scalar diffusion coefficient.
  - **Plain English:** higher values make color/smoke spread faster.
- `SOLVER_ITERATIONS`
  - **Technical:** iterative convergence budget.
  - **Plain English:** more loops improve quality, cost more CPU.
- `VORTICITY_CONFINEMENT (ε)`
  - **Technical:** strength of confinement force.
  - **Plain English:** boosts small swirls/curls.

---

## 10) End-to-end mental model

One frame is essentially:

1. Add fluid and force.
2. Make velocity physically plausible (diffuse + project).
3. Move velocity through itself (advect).
4. Recover swirl detail (vorticity) and project again.
5. Move and mix color fields.
6. Re-apply wall/boundary constraints.

### In normal terms

- **Inject → smooth → balance pressure → move → add swirl → move color → clean up walls**.

This order gives a practical balance:

- **Stable numerics** from implicit solves + semi-Lagrangian transport + projection.
- **Visual richness** from vorticity confinement.
- **Physical plausibility** from divergence control and obstacle enforcement.
