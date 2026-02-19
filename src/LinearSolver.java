import java.util.stream.IntStream;

/**
 * Iterative solver for grid-based fluid equations.
 *
 * <p>This class solves the discrete Poisson / diffusion equations that arise in
 * fluid simulation when applying viscosity (diffusion) or enforcing
 * incompressibility (pressure projection)</p>
 */
public class LinearSolver {
    /**
     * Number of solver iterations per call
     */
    private final int iterationCount;

    /**
     * Reused scratch buffer used by the Jacobi update
     */
    private float[] scratchField;

    /**
     * Constructs a linear solver with a fixed iteration count
     *
     * @param iterationCount number of relaxation iterations
     */
    public LinearSolver(int iterationCount) {
        if (iterationCount <= 0) {
            throw new IllegalArgumentException("iterationCount must be > 0");
        }
        this.iterationCount = iterationCount;
    }

    /**
     * Returns the configured iteration count
     */
    public int getIterationCount() {
        return iterationCount;
    }

    /**
     * Iteratively solves a diffusion/projection linear system for one field
     *
     * <p>This method repeatedly averages each cell with neighbors
     * until the field becomes consistent with the requested physical rule</p>
     *
     * <p>Jacobi update written with descriptive names:</p>
     * <p><b>nextCellValue</b> =
     * (originalSourceValue + neighborInfluenceStrength * (leftNeighbor + rightNeighbor + topNeighbor + bottomNeighbor))
     * / diagonalNormalization</p>
     *
     * <p>Where:</p>
     * <ul>
     *   <li><b>neighborInfluenceStrength</b> is {@code spreadStrength}</li>
     *   <li><b>diagonalNormalization</b> is {@code normalizationValue}</li>
     * </ul>
     * @param boundaryType The boundary condition to enforce when updating ghost cells
     * @param resultField The destination array that will contain the converged solution.
     *     - This array is also used as the initial guess for the solver.
     * @param sourceField The right-hand-side (source) term of the linear system
     * @param spreadStrength The coefficient controlling how strongly neighboring cells influence
     *                       each cell during diffusion
     * @param normalizationValue The diagonal normalization term used to divide the Jacobi update
     * @param grid The grid describing spatial dimensions, indexing, and ghost-cell layout.
     * @param boundaryHandler The handler responsible for applying boundary conditions to the field
     */
    public void solve(BoundaryHandler.BoundaryType boundaryType,
                      float[] resultField,
                      float[] sourceField,
                      float spreadStrength,
                      float normalizationValue,
                      FluidGrid grid,
                      BoundaryHandler boundaryHandler) {

        // check the scratch buffer is large enough to hold a full grid
        checkScratchCapacity(resultField.length);

        // Initialize the solution with the source field (Jacobi initial guess)
        System.arraycopy(sourceField, 0, resultField, 0, sourceField.length);

        // Buffers for Jacobi iterations
        float[] currentField = resultField;
        float[] nextField = scratchField;

        // Perform Jacobi relaxation iterations
        for (int iteration = 0; iteration < iterationCount; iteration++) {

            // Local references for lambda capture
            final float[] iterationCurrent = currentField;
            final float[] iterationNext = nextField;

            // Parallelize over rows because each cell update is independent
            IntStream.rangeClosed(1, grid.height).parallel().forEach(yIndex -> {
                for (int xIndex = 1; xIndex <= grid.width; xIndex++) {

                    // Linear index of the current cell
                    int index = grid.index(xIndex, yIndex);

                    // Indices of 4-connected neighbors (ghost cells included)
                    int leftIndex  = grid.index(xIndex - 1, yIndex);
                    int rightIndex = grid.index(xIndex + 1, yIndex);
                    int upIndex    = grid.index(xIndex, yIndex + 1);
                    int downIndex  = grid.index(xIndex, yIndex - 1);

                    // Sum neighbor values for the Laplacian stencil
                    float neighborSum = iterationCurrent[leftIndex]
                            + iterationCurrent[rightIndex]
                            + iterationCurrent[upIndex]
                            + iterationCurrent[downIndex];

                    // (source + spread * neighborSum) / normalization
                    iterationNext[index] =
                            (sourceField[index] + spreadStrength * neighborSum)
                                    / normalizationValue;
                }
            });

            // Enforce boundary conditions so ghost cells are valid
            // before the next iteration reads from them
            boundaryHandler.applyBoundaries(boundaryType, nextField, grid);

            // Swap buffers
            float[] temp = currentField;
            currentField = nextField;
            nextField = temp;
        }

        // make sure the results get copied over if something went wrong
        if (currentField != resultField) {
            System.arraycopy(currentField, 0, resultField, 0, resultField.length);
        }
    }


    /** Checks the reusable temporary array is large enough for the current grid. */
    private void checkScratchCapacity(int minLength) {
        if (scratchField == null || scratchField.length < minLength) {
            scratchField = new float[minLength];
        }
    }
}
