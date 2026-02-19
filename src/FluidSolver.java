import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.stream.IntStream;

/**
 * High-level controller that advances the fluid simulation forward in time.
 *
 * Orchestrates velocity and density updates using stable fluids methodology.
 */
public class FluidSolver {
    // External injectors
    private List<FluidSource> densitySources = new ArrayList<>();
    private List<FluidEmitter> emitters = new ArrayList<>();
    private List<RadialFluidEmitter> radialEmitters = new ArrayList<>();
    private List<Vortex> vortexes = new ArrayList<>();

    // The fluid grid and parameters for it
    public final FluidGrid grid;
    public final SimulationParameters parameters;

    // The velocity and density fields
    public final VectorField velocityField;
    public final ScalarField redDensityField;
    public final ScalarField greenDensityField;
    public final ScalarField blueDensityField;

    public final ScalarField pressureField;
    public final ScalarField divergenceField;

    // The linear solver and boundary handler
    public final LinearSolver linearSolver;
    public final BoundaryHandler boundaryHandler;

    /**
     * Convenience constructor when no radial emitters or vortexes are needed.
     */
    public FluidSolver(
            FluidGrid grid,
            SimulationParameters parameters,
            List<FluidSource> densitySources,
            List<FluidEmitter> emitters
    ) {
        this(grid, parameters, densitySources, emitters, List.of(), List.of());
    }

    /**
     * Builds a solver and validates that all injectors are inside the grid
     *
     * <p>The solver owns the simulation fields and applies the Stable Fluids update
     * sequence each time {@link #step()} is called</p>
     */
    public FluidSolver(
            FluidGrid grid,
            SimulationParameters parameters,
            List<FluidSource> densitySources,
            List<FluidEmitter> emitters,
            List<RadialFluidEmitter> radialEmitters,
            List<Vortex> vortexes
    ) {
        this.grid = Objects.requireNonNull(grid, "grid must not be null");
        this.parameters = Objects.requireNonNull(parameters, "parameters must not be null");

        this.velocityField = new VectorField(grid.totalCellCount);
        this.redDensityField = new ScalarField(grid.totalCellCount);
        this.greenDensityField = new ScalarField(grid.totalCellCount);
        this.blueDensityField = new ScalarField(grid.totalCellCount);

        this.pressureField = new ScalarField(grid.totalCellCount);
        this.divergenceField = new ScalarField(grid.totalCellCount);

        this.linearSolver = new LinearSolver(parameters.getLinearSolverIterations());
        this.boundaryHandler = new BoundaryHandler();

        // Null-safe initialization
        this.densitySources = (densitySources != null)
                ? new ArrayList<>(densitySources)
                : new ArrayList<>();

        this.emitters = (emitters != null)
                ? new ArrayList<>(emitters)
                : new ArrayList<>();

        this.radialEmitters = (radialEmitters != null)
                ? new ArrayList<>(radialEmitters)
                : new ArrayList<>();

        this.vortexes = (vortexes != null)
                ? new ArrayList<>(vortexes)
                : new ArrayList<>();

        for (FluidSource source : this.densitySources) {
            Objects.requireNonNull(source, "density source must not be null");
            validateInBounds(source.gridX, source.gridY, "density source");
        }
        for (FluidEmitter emitter : this.emitters) {
            Objects.requireNonNull(emitter, "emitter must not be null");
            validateInBounds(emitter.gridX(), emitter.gridY(), "emitter");
        }
        for (RadialFluidEmitter radialEmitter : this.radialEmitters) {
            Objects.requireNonNull(radialEmitter, "radial emitter must not be null");
            validateInBounds(radialEmitter.gridX(), radialEmitter.gridY(), "radial emitter");
        }
        for (Vortex vortex : this.vortexes) {
            Objects.requireNonNull(vortex, "vortex must not be null");
            validateInBounds(vortex.gridX(), vortex.gridY(), "vortex");
        }
    }


    /**
     * Advances the simulation by one time step.
     *
     * <p>The order is: add sources, diffuse/project velocity, advect velocity, apply
     * swirl recovery, re-project, then diffuse/advect each color density channel</p>
     */
    public void step() {
        addSources();

        diffuseVelocity();
        projectVelocity();

        advectVelocity();
        applyVorticity();
        projectVelocity();

        diffuseDensity(redDensityField);
        diffuseDensity(greenDensityField);
        diffuseDensity(blueDensityField);

        advectDensity(redDensityField);
        advectDensity(greenDensityField);
        advectDensity(blueDensityField);
    }

    /** Injects density and momentum from configured sources, emitters, and vortexes. */
    private void addSources() {

        float dt = parameters.getTimeStep();

        // Loop through all the fluid sources and apply density for each of them
        for (FluidSource source : densitySources) {
            if (!grid.inBounds(source.gridX, source.gridY)) {
                continue;
            }
            int index = grid.index(source.gridX, source.gridY);
            redDensityField.readValues[index] += dt * source.strength;
            greenDensityField.readValues[index] += dt * source.strength;
            blueDensityField.readValues[index] += dt * source.strength;
        }

        // Loop through all the fluid emitters and have them apply their own density and
        // velocity functions
        for (FluidEmitter emitter : emitters) {
            if (!grid.inBounds(emitter.gridX(), emitter.gridY())) {
                continue;
            }
            emitter.applyDensity(redDensityField, greenDensityField, blueDensityField, grid, dt);
            emitter.applyVelocity(velocityField, grid);
        }

        // Loop through all the radial fluid emitters and have them apply their own density and
        // velocity functions
        for (RadialFluidEmitter radialEmitter : radialEmitters) {
            if (!grid.inBounds(radialEmitter.gridX(), radialEmitter.gridY())) {
                continue;
            }
            radialEmitter.applyDensity(redDensityField, greenDensityField, blueDensityField, grid, dt);
            radialEmitter.applyVelocity(velocityField, grid);
        }

        // Loop through all the vortexes and have them apply all their absorbtion fucnctions
        for (Vortex vortex : vortexes) {
            if (!grid.inBounds(vortex.gridX(), vortex.gridY())) {
                continue;
            }
            vortex.applyVelocity(velocityField, grid, dt);
            vortex.absorbDensity(redDensityField, grid, dt);
            vortex.absorbDensity(greenDensityField, grid, dt);
            vortex.absorbDensity(blueDensityField, grid, dt);
        }
    }

    /** Smooths velocity to model viscosity (thicker fluid spreads momentum faster) */
    private void diffuseVelocity() {
        float timeStep = parameters.getTimeStep();
        float viscosity = parameters.getViscosity();

        // Stable Fluids diffusion constant:
        // spreadStrength = (secondsPerStep * fluidViscosity) / (cellWidth^2)
        // Larger values make momentum blur across neighbors more each step
        float spreadStrength = timeStep * viscosity / (grid.cellSize * grid.cellSize);

        // Each Jacobi update divides by (centerWeight + 4 neighbor weights).
        float normalizationValue = 1.0f + 4.0f * spreadStrength;

        // Solve horizontal and vertical velocities so each cell is smoothed with it's neighbors
        linearSolver.solve(
                BoundaryHandler.BoundaryType.H_VELOCITY,
                velocityField.writeVelocityX,
                velocityField.readVelocityX,
                spreadStrength,
                normalizationValue,
                grid,
                boundaryHandler
        );

        linearSolver.solve(
                BoundaryHandler.BoundaryType.V_VELOCITY,
                velocityField.writeVelocityY,
                velocityField.readVelocityY,
                spreadStrength,
                normalizationValue,
                grid,
                boundaryHandler
        );

        // Swap the read and write buffers so the latest read buffer is accurate
        velocityField.swapBuffers();
    }

    /**
     * Enforces incompressibility by removing divergence from velocity
     *
     * <p>Makes sure flow neither creates or destroys volume in a cell</p>
     *
     * <p>Math:</p>
     * <ul>
     *   <li><b>velocityDivergence</b> = -0.5/cellWidth * ((rightVelocityX - leftVelocityX)
     *       + (topVelocityY - bottomVelocityY))</li>
     *   <li>Solve Poisson equation for <b>pressure</b> so pressure gradients cancel divergence</li>
     *   <li>Correct velocity with pressure gradient:
     *       correctedVelocityX = velocityX - 0.5/cellWidth * (pressureRight - pressureLeft)</li>
     * </ul>
     */
    private void projectVelocity() {
        // Precompute 1 / (2 * cellSize) for centered finite differences
        // The denominator for (b-a) / (2 * cellSize)
        // This is used to approximate spatial derivatives
        float centeredDifferenceDenom = 0.5f / grid.cellSize;


        // Compute velocity divergence at cell centers
        // Divergence measures local volume expansion or compression
        // A divergence-free field conserves volume (incompressible flow)
        IntStream.rangeClosed(1, grid.height).parallel().forEach(y -> {
            for (int x = 1; x <= grid.width; x++) {
                int i = grid.index(x, y);

                int left = grid.index(x - 1, y);
                int right = grid.index(x + 1, y);
                int down = grid.index(x, y - 1);
                int up = grid.index(x, y + 1);

                // Sample neighboring velocity components around this cell
                // Using centered differences gives second-order accuracy
                float hvRight = velocityField.readVelocityX[right];
                float hvLeft  = velocityField.readVelocityX[left];
                float vvUp    = velocityField.readVelocityY[up];
                float vvDown  = velocityField.readVelocityY[down];

                // Discrete divergence:
                // (du/dx + dv/dy) evaluated at the cell center
                // The negative sign matches the Poisson equation convention used later
                divergenceField.writeValues[i] =
                        -centeredDifferenceDenom * ((hvRight - hvLeft) + (vvUp - vvDown));

                // Initialize pressure guess to zero for the solver.
                pressureField.writeValues[i] = 0.0f;
            }
        });

        // Swap buffers so divergence and pressure become readable
        divergenceField.swapBuffers();
        pressureField.swapBuffers();

        // Scale divergence to match the discrete Poisson equation
        // The linear solver assumes:
        //     laplacian(pressure) = divergence * cellSize^2
        // Parallel so this doesn't take years
        IntStream.rangeClosed(1, grid.height).parallel().forEach(y -> {
            for (int x = 1; x <= grid.width; x++) {
                int i = grid.index(x, y);
                divergenceField.writeValues[i] =
                        divergenceField.readValues[i] * grid.cellSize * grid.cellSize;
            }
        });

        // Apply scalar boundary conditions to divergence and pressure
        boundaryHandler.applyBoundaries(
                BoundaryHandler.BoundaryType.SCALAR,
                divergenceField.writeValues,
                grid
        );

        boundaryHandler.applyBoundaries(
                BoundaryHandler.BoundaryType.SCALAR,
                pressureField.readValues,
                grid
        );

        // Solve Poisson equation for pressure
        // This finds the pressure field whose gradient cancels divergence
        linearSolver.solve(
                BoundaryHandler.BoundaryType.SCALAR,
                pressureField.readValues,     // unknown (pressure)
                divergenceField.writeValues,  // right-hand side (scaled divergence)
                1.0f,                         // center coefficient
                4.0f,                         // neighbor sum coefficient (2D grid)
                grid,
                boundaryHandler
        );

        // Subtract pressure gradient from velocity
        // This removes divergence and enforces incompressibility
        IntStream.rangeClosed(1, grid.height).parallel().forEach(y -> {
            for (int x = 1; x <= grid.width; x++) {
                int i = grid.index(x, y);

                int left = grid.index(x - 1, y);
                int right = grid.index(x + 1, y);
                int down = grid.index(x, y - 1);
                int up = grid.index(x, y + 1);

                // Centered finite differences of pressure
                float pRight = pressureField.readValues[right];
                float pLeft  = pressureField.readValues[left];
                float pUp    = pressureField.readValues[up];
                float pDown  = pressureField.readValues[down];

                // Velocity minus pressure gradient:  (b-a) / (2 * cellSize)
                velocityField.readVelocityX[i] -= centeredDifferenceDenom * (pRight - pLeft);
                velocityField.readVelocityY[i] -= centeredDifferenceDenom * (pUp - pDown);
            }
        });

        // Enforce velocity boundary conditions after projection
        boundaryHandler.applyBoundaries(
                BoundaryHandler.BoundaryType.H_VELOCITY,
                velocityField.readVelocityX,
                grid
        );
        boundaryHandler.applyBoundaries(
                BoundaryHandler.BoundaryType.V_VELOCITY,
                velocityField.readVelocityY,
                grid
        );
    }

    /**
     * Applies vorticity confinement to re-introduce small-scale swirling motion
     * that is otherwise damped out by numerical diffusion
     *
     *  - Compute scalar vorticity (curl of the velocity field)
     *  - Compute the gradient of |vorticity| to find where spin intensity increases
     *  - Apply a force perpendicular to that gradient, scaled by the local spin
     *
     * This creates tight, visually pretty vortices without affecting large-scale flow.
     */
    private void applyVorticity() {
        float confinementStrength = parameters.getVorticityConfinement();

        // Exit early if there's no vorticity confinement
        if (confinementStrength <= 0.0f) {
            return;
        }

        float timeStep = parameters.getTimeStep();

        // Precompute 1 / (2 * cellSize) for centered finite differences
        float centeredDifferenceDenom = 0.5f / grid.cellSize;

        // Stores |curl| at each grid cell (scalar field)
        float[] curlMagnitude = new float[grid.totalCellCount];

        /*
         * Compute vorticity magnitude at each cell
         *
         * Vorticity (scalar, 2D):
         *   vorticity = d(velocityY)/dx − d(velocityX)/dy
         *
         * We store only the absolute value here because we only need its spatial
         * gradient direction (toward stronger rotational motion)
         */
        IntStream.rangeClosed(1, grid.height).parallel().forEach(y -> {
            for (int x = 1; x <= grid.width; x++) {
                int i = grid.index(x, y);

                // Neighbor indices for centered differences
                int left  = grid.index(x - 1, y);
                int right = grid.index(x + 1, y);
                int down  = grid.index(x, y - 1);
                int up    = grid.index(x, y + 1);

                // Partial derivative of horizontal velocity with respect to Y
                float dVxDy =
                        (velocityField.readVelocityX[up] -
                                velocityField.readVelocityX[down]) * centeredDifferenceDenom;

                // Partial derivative of vertical velocity with respect to X
                float dVyDx =
                        (velocityField.readVelocityY[right] -
                                velocityField.readVelocityY[left]) * centeredDifferenceDenom;

                // Store absolute vorticity magnitude
                curlMagnitude[i] = Math.abs(dVyDx - dVxDy);
            }
        });

        /*
         * PASS 2: Apply vorticity confinement force
         *
         * Confinement force:
         *   confinementForce = strength * cellSize *
         *                      (normalizedGradientOfVorticityMagnitude
         *                       x
         *                       vorticity)
         *
         * where:
         *   - strength  = user-defined confinement strength
         *   - cellSize  = grid resolution
         *   - normalizedGradientOfVorticityMagnitude points toward stronger rotation
         *   - vorticity is the signed scalar curl at the cell
         *
         * In 2D, the cross product reduces to a perpendicular vector scaled
         * by the vorticity value
         */
        IntStream.rangeClosed(1, grid.height).parallel().forEach(y -> {
            for (int x = 1; x <= grid.width; x++) {
                int i = grid.index(x, y);

                // Neighbor indices
                int left  = grid.index(x - 1, y);
                int right = grid.index(x + 1, y);
                int down  = grid.index(x, y - 1);
                int up    = grid.index(x, y + 1);

                // Recompute signed vorticity at this cell
                float dVxDy =
                        (velocityField.readVelocityX[up] -
                                velocityField.readVelocityX[down]) * centeredDifferenceDenom;

                float dVyDx =
                        (velocityField.readVelocityY[right] -
                                velocityField.readVelocityY[left]) * centeredDifferenceDenom;

                float vorticity = dVyDx - dVxDy;

                // Gradient of |vorticityMagnitude| points toward regions of stronger rotation
                // Calculated using (b-a) / (2 * cellSize)
                float gradX =
                        (curlMagnitude[right] -
                                curlMagnitude[left]) * centeredDifferenceDenom;

                float gradY =
                        (curlMagnitude[up] -
                                curlMagnitude[down]) * centeredDifferenceDenom;

                // Normalize gradient (avoid divide-by-zero with tiny value)
                float gradLength =
                        (float) Math.sqrt(gradX * gradX + gradY * gradY) + 1e-6f;

                float normalX = gradX / gradLength;
                float normalY = gradY / gradLength;

                // Scale force by user strength and grid resolution
                float forceScale = confinementStrength * grid.cellSize;

                // Rotate gradient 90 degrees to obtain tangential force direction
                // This causes the force to reinforce rotation rather than cancel it
                float forceX =  forceScale * normalY * vorticity;
                float forceY = -forceScale * normalX * vorticity;

                // Explicitly integrate confinement force into velocity
                velocityField.readVelocityX[i] += timeStep * forceX;
                velocityField.readVelocityY[i] += timeStep * forceY;
            }
        });

        // Enforce the boundary conditions
        boundaryHandler.applyBoundaries(
                BoundaryHandler.BoundaryType.H_VELOCITY,
                velocityField.readVelocityX,
                grid
        );

        boundaryHandler.applyBoundaries(
                BoundaryHandler.BoundaryType.V_VELOCITY,
                velocityField.readVelocityY,
                grid
        );
    }


    /** Diffuses one density channel, softening sharp concentration gradients. */
    private void diffuseDensity(ScalarField densityField) {
        float timeStep = parameters.getTimeStep();
        float diffusionRate = parameters.getDiffusionRate();

        // Diffusion coefficient for implicit solve:
        // spreadStrength = (dt * diffusionRate) / (cellSize squared)
        float spreadStrength = timeStep * diffusionRate / (grid.cellSize * grid.cellSize);

        // Diagonal term of the 5-point Laplacian stencil:
        // (1 + 4 * spreadStrength) * density_center
        float normalizationValue = 1.0f + 4.0f * spreadStrength;

        linearSolver.solve(
                BoundaryHandler.BoundaryType.SCALAR,
                densityField.writeValues,
                densityField.readValues,
                spreadStrength,
                normalizationValue,
                grid,
                boundaryHandler
        );

        densityField.swapBuffers();
    }

    /**
     * Moves velocity through the grid by tracing each cell backward along the flow
     *
     * <p>Backtrace equation:</p>
     * <p><b>sourcePosition</b> = <b>cellCenterPosition</b> - (secondsPerStep / cellWidth) * <b>cellVelocity</b></p>
     * <p>We then bilinearly sample velocity at that sourcePosition</p>
     */
    private void advectVelocity() {
        // Duration of the simulation step in seconds
        float timeStepSeconds = parameters.getTimeStep();
        float dtOverDx = timeStepSeconds / grid.cellSize;

        // Read buffers: velocity field from the previous frame
        float[] currentVelocityX = velocityField.readVelocityX;
        float[] currentVelocityY = velocityField.readVelocityY;

        // Source fields to sample from (same as read buffers for semi-Lagrangian advection)
        float[] sourceVelocityX = velocityField.readVelocityX;
        float[] sourceVelocityY = velocityField.readVelocityY;

        // Iterate over interior grid cells (parallelized over rows)
        IntStream.rangeClosed(1, grid.height).parallel().forEach(gridY -> {
            for (int gridX = 1; gridX <= grid.width; gridX++) {

                // Convert 2D grid coordinates into 1D array index
                int cellIndex = grid.index(gridX, gridY);

                // Semi-Lagrangian backtrace
                // Trace the velocity at this cell backward in time to find
                // where the fluid originated in the previous frame
                float sourceX = gridX - dtOverDx * currentVelocityX[cellIndex];
                float sourceY = gridY - dtOverDx * currentVelocityY[cellIndex];

                // Clamp source position to valid sampling range
                // (cell centers are offset by 0.5)
                sourceX = clamp(sourceX, 0.5f, grid.width + 0.5f);
                sourceY = clamp(sourceY, 0.5f, grid.height + 0.5f);

                // Sample the X component of velocity at the backtraced position
                velocityField.writeVelocityX[cellIndex] =
                        bilinearSample(sourceVelocityX, sourceX, sourceY);

                // Sample the Y component of velocity at the backtraced position
                velocityField.writeVelocityY[cellIndex] =
                        bilinearSample(sourceVelocityY, sourceX, sourceY);
            }
        });

        // Enforce horizontal velocity boundary conditions
        boundaryHandler.applyBoundaries(
                BoundaryHandler.BoundaryType.H_VELOCITY,
                velocityField.writeVelocityX,
                grid
        );

        // Enforce vertical velocity boundary conditions
        boundaryHandler.applyBoundaries(
                BoundaryHandler.BoundaryType.V_VELOCITY,
                velocityField.writeVelocityY,
                grid
        );

        // Update the read buffer for the next step
        velocityField.swapBuffers();
    }

    /**
     * Transports one density channel along the current velocity field
     *
     * <p>Uses semi-Lagrangian advection:</p>
     * <p><b>sourceCell</b> = <b>currentCell</b> - (secondsPerStep / cellWidth) * <b>cellVelocity</b></p>
     * <p>The density value is then bilinearly sampled at that source position.</p>
     */
    private void advectDensity(ScalarField densityField) {

        // Duration of the simulation step in seconds
        float timeStepSeconds = parameters.getTimeStep();
        float dtOverDx = timeStepSeconds / grid.cellSize;

        // Read buffers: velocity field from the previous frame
        float[] velocityX = velocityField.readVelocityX;
        float[] velocityY = velocityField.readVelocityY;

        // Source density field to sample from
        float[] sourceDensity = densityField.readValues;

        // Iterate over interior grid cells (parallelized over rows)
        IntStream.rangeClosed(1, grid.height).parallel().forEach(gridY -> {
            for (int gridX = 1; gridX <= grid.width; gridX++) {

                // Convert 2D grid coordinates into 1D array index
                int cellIndex = grid.index(gridX, gridY);

                // Semi-Lagrangian backtrace
                // Trace this density sample backward through the velocity field
                float sourceX = gridX - dtOverDx * velocityX[cellIndex];
                float sourceY = gridY - dtOverDx * velocityY[cellIndex];

                // Clamp source position to valid sampling range
                // (cell centers are offset by 0.5)
                sourceX = clamp(sourceX, 0.5f, grid.width + 0.5f);
                sourceY = clamp(sourceY, 0.5f, grid.height + 0.5f);

                // Sample density at the backtraced position
                densityField.writeValues[cellIndex] =
                        bilinearSample(sourceDensity, sourceX, sourceY);
            }
        });

        // Enforce scalar boundary conditions
        boundaryHandler.applyBoundaries(
                BoundaryHandler.BoundaryType.SCALAR,
                densityField.writeValues,
                grid
        );

        // Update the read buffer for the next step
        densityField.swapBuffers();
    }

    /** Clamps a value to [min, max] */
    private static float clamp(float value, float min, float max) {
        if (value < min) return min;
        if (value > max) return max;
        return value;
    }

    /**
     * Samples a grid field at a non-integer location using bilinear interpolation
     *
     * <p>interpolatedValue = lerpY( lerpX(bottomLeft, bottomRight),
     * lerpX(topLeft, topRight) ) using fractional offsets inside the cell.</p>
     */
    private float bilinearSample(float[] field, float sampleX, float sampleY) {
        int x0 = (int) sampleX;
        int y0 = (int) sampleY;

        int x1 = x0 + 1;
        int y1 = y0 + 1;

        float sx = sampleX - x0;
        float sy = sampleY - y0;

        int i00 = grid.index(x0, y0);
        int i10 = grid.index(x1, y0);
        int i01 = grid.index(x0, y1);
        int i11 = grid.index(x1, y1);

        float v00 = field[i00];
        float v10 = field[i10];
        float v01 = field[i01];
        float v11 = field[i11];

        float lerpX0 = v00 + sx * (v10 - v00);
        float lerpX1 = v01 + sx * (v11 - v01);

        return lerpX0 + sy * (lerpX1 - lerpX0);
    }

    /** Registers an additional point density source during runtime */
    public void addDensitySource(FluidSource source) {
        Objects.requireNonNull(source, "density source must not be null");
        validateInBounds(source.gridX, source.gridY, "density source");
        densitySources.add(source);
    }

    /** Registers an additional directional emitter during runtime */
    public void addEmitter(FluidEmitter emitter) {
        Objects.requireNonNull(emitter, "emitter must not be null");
        validateInBounds(emitter.gridX(), emitter.gridY(), "emitter");
        emitters.add(emitter);
    }

    /** Throws if a user-provided object sits outside the active simulation area */
    private void validateInBounds(int gridX, int gridY, String objectType) {
        if (!grid.inBounds(gridX, gridY)) {
            throw new IllegalArgumentException(objectType + " out of bounds: (" + gridX + ", " + gridY + ")");
        }
    }
}
