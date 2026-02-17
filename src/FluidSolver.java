import java.util.ArrayList;
import java.util.List;

/**
 * High-level controller that advances the fluid simulation forward in time.
 *
 * Orchestrates velocity and density updates using stable fluids methodology.
 */
public class FluidSolver {
    // External injectors
    private List<FluidSource> densitySources = new ArrayList<>();
    private List<FluidEmitter> emitters = new ArrayList<>();

    public final FluidGrid grid;
    public final SimulationParameters parameters;

    public final VectorField velocityField;
    public final ScalarField densityField;

    public final ScalarField pressureField;
    public final ScalarField divergenceField;

    public final LinearSolver linearSolver;
    public final BoundaryHandler boundaryHandler;

    public FluidSolver(
            FluidGrid grid,
            SimulationParameters parameters,
            List<FluidSource> densitySources,
            List<FluidEmitter> emitters
    ) {
        this.grid = grid;
        this.parameters = parameters;

        this.velocityField = new VectorField(grid.totalCellCount);
        this.densityField = new ScalarField(grid.totalCellCount);

        this.pressureField = new ScalarField(grid.totalCellCount);
        this.divergenceField = new ScalarField(grid.totalCellCount);

        this.linearSolver = new LinearSolver(parameters.linearSolverIterations);
        this.boundaryHandler = new BoundaryHandler();

        // Null-safe initialization
        this.densitySources = (densitySources != null)
                ? densitySources
                : new ArrayList<>();

        this.emitters = (emitters != null)
                ? emitters
                : new ArrayList<>();
    }

    /**
     * Advances the simulation by one time step.
     *
     * Stable Fluids canonical order:
     *
     * 1) Add sources
     * 2) Diffuse velocity (viscosity)
     * 3) Project velocity (remove compression)
     * 4) Advect velocity
     * 5) Project velocity again
     * 6) Diffuse density
     * 7) Advect density
     */
    public void step() {
        addSources();

        diffuseVelocity();
        projectVelocity();

        advectVelocity();
        projectVelocity();

        diffuseDensity();
        advectDensity();
    }

    /**
     * Adds external sources (emitters, forces, density injection).
     *
     * English equation:
     *     newValue = oldValue + timeStep * addedAmount
     *
     * Applies to:
     * - density
     * - velocity
     */
    private void addSources() {

        float dt = parameters.timeStep;

        // --- Density-only sources ---
        for (FluidSource source : densitySources) {
            int index = grid.index(source.gridX, source.gridY);
            densityField.readValues[index] += dt * source.strength;
        }

        // --- Emitters: density + momentum ---
        for (FluidEmitter emitter : emitters) {

            int index = grid.index(emitter.gridX, emitter.gridY);

            // Density injection
            densityField.readValues[index] += dt * emitter.densityRate;

            // Momentum injection
            float vx = (float) Math.cos(emitter.angleRadians) * emitter.emissionSpeed;
            float vy = (float) Math.sin(emitter.angleRadians) * emitter.emissionSpeed;

            velocityField.readVelocityX[index] += dt * vx;
            velocityField.readVelocityY[index] += dt * vy;
        }
    }


    /**
     * Diffuses velocity using viscosity (implicit solve).
     *
     * English equation (per cell):
     *
     * newVelocity =
     * (
     *     oldVelocity
     *     + spreadStrength * (left + right + up + down)
     * )
     * ÷ normalizationValue
     *
     * spreadStrength = timeStep * viscosity * gridWidth * gridHeight
     * normalizationValue = 1 + 4 * spreadStrength
     */
    private void diffuseVelocity() {
        float timeStep = parameters.timeStep;
        float viscosity = parameters.viscosity;

        // Matches your LinearSolver docs: "timeStep * diffusion * gridSize * gridSize"
        // Since your grid may be non-square, use width * height as a reasonable scale factor.
        float spreadStrength = timeStep * viscosity * grid.width * grid.height;
        float normalizationValue = 1.0f + 4.0f * spreadStrength;

        // Diffuse X velocity: writeVelocityX becomes the diffused result
        linearSolver.solve(
                BoundaryHandler.BoundaryType.H_VELOCITY,
                velocityField.writeVelocityX,
                velocityField.readVelocityX,
                spreadStrength,
                normalizationValue,
                grid,
                boundaryHandler
        );

        // Diffuse Y velocity
        linearSolver.solve(
                BoundaryHandler.BoundaryType.V_VELOCITY,
                velocityField.writeVelocityY,
                velocityField.readVelocityY,
                spreadStrength,
                normalizationValue,
                grid,
                boundaryHandler
        );

        // Commit diffused velocity
        velocityField.swapBuffers();
    }

    /**
     * Projects velocity to enforce incompressibility.
     * 
     * Projection removes compression by computing a pressure field
     * whose gradient exactly cancels any local inflow or outflow of fluid
     *
     * Step A: compute "compression" (divergence) in each cell:
     *
     * compression =
     *   -0.5 * (
     *       rightVelocityX - leftVelocityX
     *     + topVelocityY   - bottomVelocityY
     *   )
     *
     * Step B: solve for pressure using Poisson relaxation:
     *
     * pressure =
     * (
     *     compression
     *     + leftPressure + rightPressure + topPressure + bottomPressure
     * )
     * ÷ 4
     *
     * Step C: subtract pressure gradient from velocity:
     *
     * velocityX -= 0.5 * (rightPressure - leftPressure)
     * velocityY -= 0.5 * (topPressure   - bottomPressure)
     */
    private void projectVelocity() {
        // 1) Compute divergence (compression)
        for (int y = 1; y <= grid.height; y++) {
            for (int x = 1; x <= grid.width; x++) {
                int i = grid.index(x, y);

                int left = grid.index(x - 1, y);
                int right = grid.index(x + 1, y);
                int down = grid.index(x, y - 1);
                int up = grid.index(x, y + 1);

                float uRight = velocityField.readVelocityX[right];
                float uLeft  = velocityField.readVelocityX[left];
                float vUp    = velocityField.readVelocityY[up];
                float vDown  = velocityField.readVelocityY[down];

                // "compression" per your English equation
                divergenceField.writeValues[i] = -0.5f * ((uRight - uLeft) + (vUp - vDown));

                // Initialize pressure guess to 0 each projection step
                pressureField.writeValues[i] = 0.0f;
            }
        }

        divergenceField.swapBuffers();
        pressureField.swapBuffers();

        boundaryHandler.applyBoundaries(BoundaryHandler.BoundaryType.SCALAR, divergenceField.readValues, grid);
        boundaryHandler.applyBoundaries(BoundaryHandler.BoundaryType.SCALAR, pressureField.readValues, grid);

        // 2) Solve for pressure: pressure = (divergence + neighbors) / 4
        // This matches your LinearSolver explanation exactly: spreadStrength=1, normalizationValue=4
        linearSolver.solve(
                BoundaryHandler.BoundaryType.SCALAR,
                pressureField.readValues,        // resultField updated in-place
                divergenceField.readValues,      // sourceField
                1.0f,                            // spreadStrength
                4.0f,                            // normalizationValue
                grid,
                boundaryHandler
        );

        // 3) Subtract pressure gradient from velocity (in-place on authoritative velocity)
        for (int y = 1; y <= grid.height; y++) {
            for (int x = 1; x <= grid.width; x++) {
                int i = grid.index(x, y);

                int left = grid.index(x - 1, y);
                int right = grid.index(x + 1, y);
                int down = grid.index(x, y - 1);
                int up = grid.index(x, y + 1);

                float pRight = pressureField.readValues[right];
                float pLeft  = pressureField.readValues[left];
                float pUp    = pressureField.readValues[up];
                float pDown  = pressureField.readValues[down];

                velocityField.readVelocityX[i] -= 0.5f * (pRight - pLeft);
                velocityField.readVelocityY[i] -= 0.5f * (pUp - pDown);
            }
        }

        boundaryHandler.applyBoundaries(BoundaryHandler.BoundaryType.H_VELOCITY, velocityField.readVelocityX, grid);
        boundaryHandler.applyBoundaries(BoundaryHandler.BoundaryType.V_VELOCITY, velocityField.readVelocityY, grid);
    }

    /**
     * For each grid cell, we trace backward through the current velocity field
     * to find where the fluid occupying this cell came from during the previous
     * timestep. We then sample the previous velocity field at that location
     * and assign it as the new velocity for this cell.
     *
     * This backward-tracing approach is unconditionally stable and allows
     * large timesteps without numerical explosion
     */
    private void advectVelocity() {

        // Duration of one simulation step in seconds
        float timeStepSeconds = parameters.timeStep;

        /*
         * Velocity in this solver is expressed in "domain units per second"
         * where a value of 1.0 means "cross the entire simulation domain in one second"
         *
         * Grid coordinates are expressed in discrete grid cells
         *
         * These factors convert a velocity value into a displacement measured
         * in grid cells for the current timestep
         */
        float velocityToGridCellsX = timeStepSeconds * grid.width;
        float velocityToGridCellsY = timeStepSeconds * grid.height;

        /*
         * Current velocity field after diffusion and projection
         * These velocities determine how we trace backward through the flow
         */
        float[] currentVelocityX = velocityField.readVelocityX;
        float[] currentVelocityY = velocityField.readVelocityY;

        /*
         * Source velocity field from the previous timestep
         * We sample from this field when determining the new velocity values
         *
         * At this point in the solver, the previous and current velocity buffers
         * refer to the same arrays, but they represent distinct states:
         *  - "currentVelocity" defines where fluid moves
         *  - "sourceVelocity" defines what velocity the fluid had
         */
        float[] sourceVelocityX = velocityField.readVelocityX;
        float[] sourceVelocityY = velocityField.readVelocityY;

        // Compute new velocity values for every interior grid cell
        for (int gridY = 1; gridY <= grid.height; gridY++) {
            for (int gridX = 1; gridX <= grid.width; gridX++) {

                int cellIndex = grid.index(gridX, gridY);

                /*
                 * Backtrace: determine where the fluid currently in this cell
                 * was located one timestep ago by moving backward along the
                 * velocity vector.
                 */
                float sourceX = gridX - velocityToGridCellsX * currentVelocityX[cellIndex];
                float sourceY = gridY - velocityToGridCellsY * currentVelocityY[cellIndex];

                /*
                 * Clamp the sampling position to the valid interior of the grid.
                 * Cell centers lie at integer coordinates; the 0.5 offset ensures
                 * we never sample outside the domain or into invalid ghost regions.
                 */
                sourceX = clamp(sourceX, 0.5f, grid.width  + 0.5f);
                sourceY = clamp(sourceY, 0.5f, grid.height + 0.5f);

                /*
                 * Sample the previous velocity field at the backtraced position
                 * using bilinear interpolation. This smooth interpolation avoids
                 * grid-aligned artifacts and numerical noise
                 */
                velocityField.writeVelocityX[cellIndex] =
                        bilinearSample(sourceVelocityX, sourceX, sourceY);

                velocityField.writeVelocityY[cellIndex] =
                        bilinearSample(sourceVelocityY, sourceX, sourceY);
            }
        }

        /*
         * Apply physical boundary conditions to enforce no-penetration at the
         * simulation domain boundaries:
         *  - Horizontal velocity is negated at vertical walls
         *  - Vertical velocity is negated at horizontal walls
         */
        boundaryHandler.applyBoundaries(
                BoundaryHandler.BoundaryType.H_VELOCITY,
                velocityField.writeVelocityX,
                grid
        );
        boundaryHandler.applyBoundaries(
                BoundaryHandler.BoundaryType.V_VELOCITY,
                velocityField.writeVelocityY,
                grid
        );

        /*
         * Commit the newly computed velocity field
         * After this swap, the advected velocities become the authoritative state
         */
        velocityField.swapBuffers();
    }


    /**
     * Diffuses density.
     *
     * Same equation as velocity diffusion, but uses diffusionRate.
     *
     * English equation (per cell):
     *
     * newDensity =
     * (
     *     oldDensity
     *     + spreadStrength * (left + right + up + down)
     * )
     * ÷ normalizationValue
     *
     * spreadStrength =
     *     timeStep * diffusionRate * gridWidth * gridHeight
     */
    private void diffuseDensity() {
        float timeStep = parameters.timeStep;
        float diffusionRate = parameters.diffusionRate;

        // Same scaling approach as velocity diffusion (handles non-square grids reasonably)
        float spreadStrength = timeStep * diffusionRate * grid.width * grid.height;
        float normalizationValue = 1.0f + 4.0f * spreadStrength;

        // Write buffer becomes the diffused result
        linearSolver.solve(
                BoundaryHandler.BoundaryType.SCALAR,
                densityField.writeValues,
                densityField.readValues,
                spreadStrength,
                normalizationValue,
                grid,
                boundaryHandler
        );

        // Commit diffused density
        densityField.swapBuffers();
    }

    /**
     * Advects density through the velocity field.
     *
     * Answers:
     * "Where did the density at this cell come from?"
     *
     * English algorithm (per cell):
     *
     * previousX = currentX - timeStep * velocityX
     * previousY = currentY - timeStep * velocityY
     *
     * newDensity =
     *     bilinearInterpolation(
     *         previousDensityField,
     *         previousX,
     *         previousY
     *     )
     */
    private void advectDensity() {

        float timeStepSeconds = parameters.timeStep;

        // Convert "domain units per second" velocity into "grid cells per timestep"
        float velocityToGridCellsX = timeStepSeconds * grid.width;
        float velocityToGridCellsY = timeStepSeconds * grid.height;

        // Velocity field drives the backtrace
        float[] velocityX = velocityField.readVelocityX;
        float[] velocityY = velocityField.readVelocityY;

        // Density source is the previous (authoritative) density buffer
        float[] sourceDensity = densityField.readValues;

        // Compute new density values for every interior cell
        for (int gridY = 1; gridY <= grid.height; gridY++) {
            for (int gridX = 1; gridX <= grid.width; gridX++) {

                int cellIndex = grid.index(gridX, gridY);

                // Backtrace to find where this cell's density came from
                float sourceX = gridX - velocityToGridCellsX * velocityX[cellIndex];
                float sourceY = gridY - velocityToGridCellsY * velocityY[cellIndex];

                // Clamp into valid sampling region (includes ghost cells)
                sourceX = clamp(sourceX, 0.5f, grid.width  + 0.5f);
                sourceY = clamp(sourceY, 0.5f, grid.height + 0.5f);

                // Bilinear sample the previous density field at that location
                densityField.writeValues[cellIndex] = bilinearSample(sourceDensity, sourceX, sourceY);
            }
        }

        // Enforce scalar boundary conditions on density
        boundaryHandler.applyBoundaries(
                BoundaryHandler.BoundaryType.SCALAR,
                densityField.writeValues,
                grid
        );

        // Commit advected density
        densityField.swapBuffers();
    }


    /* ------------------------------ Helpers ------------------------------ */

    private static float clamp(float value, float min, float max) {
        if (value < min) return min;
        if (value > max) return max;
        return value;
    }

    /**
     * Samples a scalar field at a non-integer position using bilinear interpolation
     *
     * This method treats the scalar field as a continuous function defined over
     * the discrete grid by linearly blending the values of the four nearest
     * grid cell centers
     *
     * Coordinate system:
     * - The grid is indexed in "grid coordinates"
     * - Cell centers lie at integer positions: (1,1), (2,1), ..., (width,height)
     * - sampleX and sampleY may be fractional and represent a point inside the grid
     *   where a value is requested
     *
     * Valid sampling range:
     * - sampleX from [0.5, width  + 0.5]
     * - sampleY from [0.5, height + 0.5]
     *
     * These bounds ensure that the four neighboring sample points always exist,
     * including near the domain boundaries where ghost cells are present
     *
     * This is used in semi-Lagrangian advection:
     * after tracing backward through the velocity field, we use bilinear
     * interpolation to smoothly sample the field at the backtraced position
     *
     * @param field   scalar field values stored in a 1D array indexed by the grid
     * @param sampleX x-coordinate in grid space where the field should be sampled
     * @param sampleY y-coordinate in grid space where the field should be sampled
     * @return interpolated scalar value at (sampleX, sampleY)
     */
    private float bilinearSample(float[] field, float sampleX, float sampleY) {

        /*
         * Identify the integer grid cell containing the sampling point
         * Casting to int truncates toward zero computing floor()
         * for positive coordinates
         *
         * (x0, y0) is the lower-left corner of the interpolation cell
         */
        int x0 = (int) sampleX;
        int y0 = (int) sampleY;

        // The neighboring cell one unit to the right and one unit above
        int x1 = x0 + 1;
        int y1 = y0 + 1;

        /*
         * Compute interpolation weights in each dimension.
         *
         * sx = how far sampleX lies between x0 and x1
         * sy = how far sampleY lies between y0 and y1
         *
         * Both values lie in the range [0, 1].
         */
        float sx = sampleX - x0;
        float sy = sampleY - y0;

        /*
         * Convert 2D grid coordinates into 1D array indices
         *
         * The four indices correspond to the corners of the cell
         * surrounding the sampling point:
         *
         *   (x0, y1) ---- (x1, y1)
         *      |            |
         *      |            |
         *   (x0, y0) ---- (x1, y0)
         */
        int i00 = grid.index(x0, y0); // lower left
        int i10 = grid.index(x1, y0); // lower right
        int i01 = grid.index(x0, y1); // upper left
        int i11 = grid.index(x1, y1); // upper right

        /*
         * Fetch the scalar values stored at each corner of the cell
         */
        float v00 = field[i00];
        float v10 = field[i10];
        float v01 = field[i01];
        float v11 = field[i11];

        /*
         * Perform linear interpolation along the x-axis at the bottom and top
         * edges of the cell.
         *
         * lerpX0 interpolates between (x0, y0) and (x1, y0)
         * lerpX1 interpolates between (x0, y1) and (x1, y1)
         */
        float lerpX0 = v00 + sx * (v10 - v00);
        float lerpX1 = v01 + sx * (v11 - v01);

        /*
         * Perform linear interpolation along the y-axis between the two
         * x-interpolated values to obtain the final result.
         */
        return lerpX0 + sy * (lerpX1 - lerpX0);
    }

    public void addDensitySource(FluidSource source) {
        densitySources.add(source);
    }

    public void addEmitter(FluidEmitter emitter) {
        emitters.add(emitter);
    }

}
