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
        int stride = grid.width + 2;
        int width = grid.width;
        int height = grid.height;
        float invNormalizationValue = 1.0f / normalizationValue;

        for (int iteration = 0; iteration < iterationCount; iteration++) {
            final float[] iterationCurrent = currentField;
            final float[] iterationNext = nextField;

            for (int yIndex = 1; yIndex <= height; yIndex++) {
                int rowBase = yIndex * stride;
                for (int xIndex = 1; xIndex <= width; xIndex++) {
                    int index = rowBase + xIndex;

                    int leftIndex = index - 1;
                    int rightIndex = index + 1;
                    int upIndex = index + stride;
                    int downIndex = index - stride;

                    float neighborSum = iterationCurrent[leftIndex]
                            + iterationCurrent[rightIndex]
                            + iterationCurrent[upIndex]
                            + iterationCurrent[downIndex];

                    iterationNext[index] = (sourceField[index] + spreadStrength * neighborSum) * invNormalizationValue;
                }
            }

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
