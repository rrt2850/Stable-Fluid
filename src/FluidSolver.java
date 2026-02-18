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

    public final FluidGrid grid;
    public final SimulationParameters parameters;

    public final VectorField velocityField;
    public final ScalarField redDensityField;
    public final ScalarField greenDensityField;
    public final ScalarField blueDensityField;

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
        this(grid, parameters, densitySources, emitters, List.of(), List.of());
    }

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

    public void step() {
        addSources();

        diffuseVelocity();
        projectVelocity();

        advectVelocity();
        applyVorticityConfinement();
        projectVelocity();

        diffuseDensity(redDensityField);
        diffuseDensity(greenDensityField);
        diffuseDensity(blueDensityField);

        advectDensity(redDensityField);
        advectDensity(greenDensityField);
        advectDensity(blueDensityField);
    }

    private void addSources() {

        float dt = parameters.getTimeStep();

        for (FluidSource source : densitySources) {
            if (!grid.inBounds(source.gridX, source.gridY)) {
                continue;
            }
            int index = grid.index(source.gridX, source.gridY);
            redDensityField.readValues[index] += dt * source.strength;
            greenDensityField.readValues[index] += dt * source.strength;
            blueDensityField.readValues[index] += dt * source.strength;
        }

        for (FluidEmitter emitter : emitters) {
            if (!grid.inBounds(emitter.gridX(), emitter.gridY())) {
                continue;
            }
            emitter.applyDensity(redDensityField, greenDensityField, blueDensityField, grid, dt);
            emitter.applyVelocity(velocityField, grid);
        }

        for (RadialFluidEmitter radialEmitter : radialEmitters) {
            if (!grid.inBounds(radialEmitter.gridX(), radialEmitter.gridY())) {
                continue;
            }
            radialEmitter.applyDensity(redDensityField, greenDensityField, blueDensityField, grid, dt);
            radialEmitter.applyVelocity(velocityField, grid);
        }

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

    private void diffuseVelocity() {
        float timeStep = parameters.getTimeStep();
        float viscosity = parameters.getViscosity();

        float spreadStrength = timeStep * viscosity / (grid.cellSize * grid.cellSize);
        float normalizationValue = 1.0f + 4.0f * spreadStrength;

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

        velocityField.swapBuffers();
    }

    private void projectVelocity() {
        float invTwoCellSize = 0.5f / grid.cellSize;

        IntStream.rangeClosed(1, grid.height).parallel().forEach(y -> {
            for (int x = 1; x <= grid.width; x++) {
                int i = grid.index(x, y);

                int left = grid.index(x - 1, y);
                int right = grid.index(x + 1, y);
                int down = grid.index(x, y - 1);
                int up = grid.index(x, y + 1);

                float uRight = velocityField.readVelocityX[right];
                float uLeft = velocityField.readVelocityX[left];
                float vUp = velocityField.readVelocityY[up];
                float vDown = velocityField.readVelocityY[down];

                divergenceField.writeValues[i] = -invTwoCellSize * ((uRight - uLeft) + (vUp - vDown));
                pressureField.writeValues[i] = 0.0f;
            }
        });

        divergenceField.swapBuffers();
        pressureField.swapBuffers();

        IntStream.rangeClosed(1, grid.height).parallel().forEach(y -> {
            for (int x = 1; x <= grid.width; x++) {
                int i = grid.index(x, y);
                divergenceField.writeValues[i] = divergenceField.readValues[i] * grid.cellSize * grid.cellSize;
            }
        });
        boundaryHandler.applyBoundaries(BoundaryHandler.BoundaryType.SCALAR, divergenceField.writeValues, grid);

        boundaryHandler.applyBoundaries(BoundaryHandler.BoundaryType.SCALAR, pressureField.readValues, grid);

        linearSolver.solve(
                BoundaryHandler.BoundaryType.SCALAR,
                pressureField.readValues,
                divergenceField.writeValues,
                1.0f,
                4.0f,
                grid,
                boundaryHandler
        );

        IntStream.rangeClosed(1, grid.height).parallel().forEach(y -> {
            for (int x = 1; x <= grid.width; x++) {
                int i = grid.index(x, y);

                int left = grid.index(x - 1, y);
                int right = grid.index(x + 1, y);
                int down = grid.index(x, y - 1);
                int up = grid.index(x, y + 1);

                float pRight = pressureField.readValues[right];
                float pLeft = pressureField.readValues[left];
                float pUp = pressureField.readValues[up];
                float pDown = pressureField.readValues[down];

                velocityField.readVelocityX[i] -= invTwoCellSize * (pRight - pLeft);
                velocityField.readVelocityY[i] -= invTwoCellSize * (pUp - pDown);
            }
        });

        boundaryHandler.applyBoundaries(BoundaryHandler.BoundaryType.H_VELOCITY, velocityField.readVelocityX, grid);
        boundaryHandler.applyBoundaries(BoundaryHandler.BoundaryType.V_VELOCITY, velocityField.readVelocityY, grid);
    }

    private void advectVelocity() {

        float timeStepSeconds = parameters.getTimeStep();

        float velocityToGridCells = timeStepSeconds / grid.cellSize;

        float[] currentVelocityX = velocityField.readVelocityX;
        float[] currentVelocityY = velocityField.readVelocityY;

        float[] sourceVelocityX = velocityField.readVelocityX;
        float[] sourceVelocityY = velocityField.readVelocityY;

        IntStream.rangeClosed(1, grid.height).parallel().forEach(gridY -> {
            for (int gridX = 1; gridX <= grid.width; gridX++) {

                int cellIndex = grid.index(gridX, gridY);

                float sourceX = gridX - velocityToGridCells * currentVelocityX[cellIndex];
                float sourceY = gridY - velocityToGridCells * currentVelocityY[cellIndex];

                sourceX = clamp(sourceX, 0.5f, grid.width + 0.5f);
                sourceY = clamp(sourceY, 0.5f, grid.height + 0.5f);

                velocityField.writeVelocityX[cellIndex] =
                        bilinearSample(sourceVelocityX, sourceX, sourceY);

                velocityField.writeVelocityY[cellIndex] =
                        bilinearSample(sourceVelocityY, sourceX, sourceY);
            }
        });

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

        velocityField.swapBuffers();
    }

    private void applyVorticityConfinement() {
        float confinementStrength = parameters.getVorticityConfinement();
        if (confinementStrength <= 0.0f) {
            return;
        }

        float dt = parameters.getTimeStep();
        float invTwoCellSize = 0.5f / grid.cellSize;
        float[] curlMagnitude = new float[grid.totalCellCount];

        IntStream.rangeClosed(1, grid.height).parallel().forEach(y -> {
            for (int x = 1; x <= grid.width; x++) {
                int i = grid.index(x, y);

                int left = grid.index(x - 1, y);
                int right = grid.index(x + 1, y);
                int down = grid.index(x, y - 1);
                int up = grid.index(x, y + 1);

                float dVxDy = (velocityField.readVelocityX[up] - velocityField.readVelocityX[down]) * invTwoCellSize;
                float dVyDx = (velocityField.readVelocityY[right] - velocityField.readVelocityY[left]) * invTwoCellSize;

                curlMagnitude[i] = Math.abs(dVyDx - dVxDy);
            }
        });

        IntStream.rangeClosed(1, grid.height).parallel().forEach(y -> {
            for (int x = 1; x <= grid.width; x++) {
                int i = grid.index(x, y);

                int left = grid.index(x - 1, y);
                int right = grid.index(x + 1, y);
                int down = grid.index(x, y - 1);
                int up = grid.index(x, y + 1);

                float dVxDy = (velocityField.readVelocityX[up] - velocityField.readVelocityX[down]) * invTwoCellSize;
                float dVyDx = (velocityField.readVelocityY[right] - velocityField.readVelocityY[left]) * invTwoCellSize;
                float vorticity = dVyDx - dVxDy;

                float gradX = (curlMagnitude[right] - curlMagnitude[left]) * invTwoCellSize;
                float gradY = (curlMagnitude[up] - curlMagnitude[down]) * invTwoCellSize;

                float gradLength = (float) Math.sqrt(gradX * gradX + gradY * gradY) + 1e-6f;
                float normalX = gradX / gradLength;
                float normalY = gradY / gradLength;

                float forceScale = confinementStrength * grid.cellSize;
                float forceX = forceScale * normalY * vorticity;
                float forceY = -forceScale * normalX * vorticity;

                velocityField.readVelocityX[i] += dt * forceX;
                velocityField.readVelocityY[i] += dt * forceY;
            }
        });

        boundaryHandler.applyBoundaries(BoundaryHandler.BoundaryType.H_VELOCITY, velocityField.readVelocityX, grid);
        boundaryHandler.applyBoundaries(BoundaryHandler.BoundaryType.V_VELOCITY, velocityField.readVelocityY, grid);
    }

    private void diffuseDensity(ScalarField densityField) {
        float timeStep = parameters.getTimeStep();
        float diffusionRate = parameters.getDiffusionRate();

        float spreadStrength = timeStep * diffusionRate / (grid.cellSize * grid.cellSize);
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

    private void advectDensity(ScalarField densityField) {

        float timeStepSeconds = parameters.getTimeStep();

        float velocityToGridCells = timeStepSeconds / grid.cellSize;

        float[] velocityX = velocityField.readVelocityX;
        float[] velocityY = velocityField.readVelocityY;

        float[] sourceDensity = densityField.readValues;

        IntStream.rangeClosed(1, grid.height).parallel().forEach(gridY -> {
            for (int gridX = 1; gridX <= grid.width; gridX++) {

                int cellIndex = grid.index(gridX, gridY);

                float sourceX = gridX - velocityToGridCells * velocityX[cellIndex];
                float sourceY = gridY - velocityToGridCells * velocityY[cellIndex];

                sourceX = clamp(sourceX, 0.5f, grid.width + 0.5f);
                sourceY = clamp(sourceY, 0.5f, grid.height + 0.5f);

                densityField.writeValues[cellIndex] = bilinearSample(sourceDensity, sourceX, sourceY);
            }
        });

        boundaryHandler.applyBoundaries(
                BoundaryHandler.BoundaryType.SCALAR,
                densityField.writeValues,
                grid
        );

        densityField.swapBuffers();
    }


    private static float clamp(float value, float min, float max) {
        if (value < min) return min;
        if (value > max) return max;
        return value;
    }

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

    public void addDensitySource(FluidSource source) {
        Objects.requireNonNull(source, "density source must not be null");
        validateInBounds(source.gridX, source.gridY, "density source");
        densitySources.add(source);
    }

    public void addEmitter(FluidEmitter emitter) {
        Objects.requireNonNull(emitter, "emitter must not be null");
        validateInBounds(emitter.gridX(), emitter.gridY(), "emitter");
        emitters.add(emitter);
    }

    private void validateInBounds(int gridX, int gridY, String objectType) {
        if (!grid.inBounds(gridX, gridY)) {
            throw new IllegalArgumentException(objectType + " out of bounds: (" + gridX + ", " + gridY + ")");
        }
    }

    /**
     * Computes the root-mean-square (RMS) velocity divergence across active cells.
     *
     * <p>For an incompressible flow this value should remain close to zero after
     * projection. It is useful as a quick correctness signal for a Stable Fluids step.</p>
     */
    public float computeVelocityDivergenceRms() {
        float invTwoCellSize = 0.5f / grid.cellSize;
        float sumSquares = 0.0f;
        int activeCells = grid.width * grid.height;

        for (int y = 1; y <= grid.height; y++) {
            for (int x = 1; x <= grid.width; x++) {
                int left = grid.index(x - 1, y);
                int right = grid.index(x + 1, y);
                int down = grid.index(x, y - 1);
                int up = grid.index(x, y + 1);

                float duDx = (velocityField.readVelocityX[right] - velocityField.readVelocityX[left]) * invTwoCellSize;
                float dvDy = (velocityField.readVelocityY[up] - velocityField.readVelocityY[down]) * invTwoCellSize;

                float divergence = duDx + dvDy;
                sumSquares += divergence * divergence;
            }
        }

        return (float) Math.sqrt(sumSquares / activeCells);
    }
}
