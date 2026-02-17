/**
 * Solves linear systems arising from diffusion and pressure projection
 *
 * Uses an iterative Gauss–Seidel relaxation method
 */
public class LinearSolver {
    /**
     * Number of solver iterations per call.
     */
    private final int iterationCount;

    /**
     * Constructs a linear solver with a fixed iteration count
     *
     * @param iterationCount number of Gauss–Seidel iterations
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
     *
     * Answers: Given what each cell had last frame, and how much it wants to
     * spread out, what should each cell’s new value be so that everything distributes how it should
     *
     * aka My new value should be close to my old value, but also close to the average of my neighbors
     *
     * New value = ( old value + spread strength * (left + right + up + down) ) / normalization value
     *
     * spreadStrength = timeStep * diffusion * gridSize * gridSize (how much the neighbors influence each cell)
     * normalizationValue = 1 + 4 * spreadStrength (aka total influence acting on the current cell)
     *  - 1 = influence of the cell’s own previous value
     *  - 4 = number of directions
     *
     *  resultField(row, col) =  (sourceField(row, col) + spreadStrength * neighbors) / c
     *
     * @param boundaryType  field type for boundary handling
     * @param resultField  output array storing the solution (the values for the next frame)
     * @param sourceField  input array containing the right-hand side (the values from the last frame)
     * @param spreadStrength diffusion or viscosity coefficient (how much the neighbors influence each cell)
     * @param normalizationValue normalization constant (how much total influence there is)
     * @param grid         simulation grid
     * @param boundaryHandler boundary condition handler
     */
    public void solve(BoundaryHandler.BoundaryType boundaryType,
                      float[] resultField,
                      float[] sourceField,
                      float spreadStrength,
                      float normalizationValue,
                      FluidGrid grid,
                      BoundaryHandler boundaryHandler) {

        int index, leftIndex, rightIndex, upIndex, downIndex;

        for(int iteration = 0; iteration < iterationCount; iteration++){
            for(int yIndex = 1; yIndex <= grid.height; yIndex++){
                for(int xIndex = 1; xIndex <= grid.width; xIndex++){
                    index = grid.index(xIndex, yIndex);

                    leftIndex = grid.index(xIndex - 1, yIndex);
                    rightIndex = grid.index(xIndex + 1, yIndex);
                    upIndex = grid.index(xIndex, yIndex + 1);
                    downIndex = grid.index(xIndex, yIndex - 1);

                    // Sum up the most recent surrounding fields (sourceField has the original values)
                    float neighorSum = resultField[leftIndex] +
                            resultField[rightIndex] +
                            resultField[upIndex] +
                            resultField[downIndex];

                    resultField[index] = (sourceField[index] + spreadStrength * neighorSum) / normalizationValue;
                }
            }

            // Update the ghost cells using the new values
            boundaryHandler.applyBoundaries(boundaryType, resultField, grid);
        }
    }
}
