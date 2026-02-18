import java.util.stream.IntStream;

/**
 * Solves linear systems arising from diffusion and pressure projection.
 *
 * Uses iterative relaxation to converge toward the solution.
 */
public class LinearSolver {
    /**
     * Number of solver iterations per call.
     */
    private final int iterationCount;

    /**
     * Reused scratch buffer used by the Jacobi update.
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
     * Returns the configured iteration count.
     */
    public int getIterationCount() {
        return iterationCount;
    }

    /**
     * Iteratively solves a diffusion/projection linear system for one field.
     *
     * <p>In plain language: this method repeatedly averages each cell with neighbors
     * until the field becomes consistent with the requested physical rule.</p>
     *
     * <p>Jacobi update written with descriptive names:</p>
     * <p><b>nextCellValue</b> =
     * (originalSourceValue + neighborInfluenceStrength * (leftNeighbor + rightNeighbor + topNeighbor + bottomNeighbor))
     * / diagonalNormalization.</p>
     *
     * <p>Where:</p>
     * <ul>
     *   <li><b>neighborInfluenceStrength</b> is {@code spreadStrength}</li>
     *   <li><b>diagonalNormalization</b> is {@code normalizationValue}</li>
     * </ul>
     */
    public void solve(BoundaryHandler.BoundaryType boundaryType,
                      float[] resultField,
                      float[] sourceField,
                      float spreadStrength,
                      float normalizationValue,
                      FluidGrid grid,
                      BoundaryHandler boundaryHandler) {
        ensureScratchCapacity(resultField.length);

        System.arraycopy(sourceField, 0, resultField, 0, sourceField.length);
        float[] currentField = resultField;
        float[] nextField = scratchField;

        for (int iteration = 0; iteration < iterationCount; iteration++) {
            final float[] iterationCurrent = currentField;
            final float[] iterationNext = nextField;

            IntStream.rangeClosed(1, grid.height).parallel().forEach(yIndex -> {
                for (int xIndex = 1; xIndex <= grid.width; xIndex++) {
                    int index = grid.index(xIndex, yIndex);

                    int leftIndex = grid.index(xIndex - 1, yIndex);
                    int rightIndex = grid.index(xIndex + 1, yIndex);
                    int upIndex = grid.index(xIndex, yIndex + 1);
                    int downIndex = grid.index(xIndex, yIndex - 1);

                    // Add the 4-connected neighbor values used by the stencil.
                    float neighborSum = iterationCurrent[leftIndex]
                            + iterationCurrent[rightIndex]
                            + iterationCurrent[upIndex]
                            + iterationCurrent[downIndex];

                    // One Jacobi relaxation update for this cell.
                    iterationNext[index] = (sourceField[index] + spreadStrength * neighborSum) / normalizationValue;
                }
            });

            // Keep ghost cells valid after each iteration so the next stencil read is correct at edges.
            boundaryHandler.applyBoundaries(boundaryType, nextField, grid);

            float[] temp = currentField;
            currentField = nextField;
            nextField = temp;
        }

        if (currentField != resultField) {
            System.arraycopy(currentField, 0, resultField, 0, resultField.length);
        }
    }

    /** Ensures the reusable temporary array is large enough for the current grid. */
    private void ensureScratchCapacity(int minLength) {
        if (scratchField == null || scratchField.length < minLength) {
            scratchField = new float[minLength];
        }
    }
}
