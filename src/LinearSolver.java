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
     * Solves a linear equation of the form:
     *     x - a * Laplacian(x) = b
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

                    float neighborSum = iterationCurrent[leftIndex]
                            + iterationCurrent[rightIndex]
                            + iterationCurrent[upIndex]
                            + iterationCurrent[downIndex];

                    iterationNext[index] = (sourceField[index] + spreadStrength * neighborSum) / normalizationValue;
                }
            });

            boundaryHandler.applyBoundaries(boundaryType, nextField, grid);

            float[] temp = currentField;
            currentField = nextField;
            nextField = temp;
        }

        if (currentField != resultField) {
            System.arraycopy(currentField, 0, resultField, 0, resultField.length);
        }
    }

    private void ensureScratchCapacity(int minLength) {
        if (scratchField == null || scratchField.length < minLength) {
            scratchField = new float[minLength];
        }
    }
}
