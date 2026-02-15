import java.util.Arrays;

/**
 * Represents a scalar quantity defined on the fluid grid
 *
 * Typical uses include density, pressure, and divergence.
 * Double buffering for better computation complexity
 */
public class ScalarField {

    /**
     * Authoritative scalar values for each grid cell
     *
     * Solver steps must only read from this buffer
     */
    public float[] readValues;

    /**
     * Scalar write buffer
     *
     * Solver steps write newly computed scalar values into this buffer
     * Its contents are undefined until a solver step completes
     */
    public float[] writeValues;

    /**
     * Allocates scalar buffers sized to the grid
     *
     * @param totalCellCount total number of grid cells including ghost cells
     */
    public ScalarField(int totalCellCount) {
        readValues = new float[totalCellCount];
        writeValues = new float[totalCellCount];
    }

    /**
     * Sets all scalar values to zero.
     */
    public void clear() {
        Arrays.fill(readValues, 0f);
    }

    /**
     * Swaps the read and write scalar buffers.
     *
     * After calling this method:
     * - The newly computed scalar values become authoritative
     * - The old read buffer becomes available for reuse as a write buffer
     *
     * This operation runs in constant time and performs no allocation
     */
    public void swapBuffers() {
        float[] temp = readValues;
        readValues = writeValues;
        writeValues = temp;
    }
}
