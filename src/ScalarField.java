import java.util.Arrays;

/**
 * Represents a scalar quantity defined on the fluid grid
 *
 * Typical uses include density, pressure, and divergence.
 * Double buffering for better computation complexity
 */
public class ScalarField {

    public float[] readValues;
    public float[] writeValues;

    /**
     * Creates a scalar field with double buffers sized to the grid.
     *
     * @param totalCellCount number of active + ghost cells in the underlying grid
     */
    public ScalarField(int totalCellCount) {
        readValues = new float[totalCellCount];
        writeValues = new float[totalCellCount];
    }

    /**
     * Resets both buffers to zero
     *
     * <p>This is useful when restarting a simulation or clearing temporary fields</p>
     */
    public void clear() {
        Arrays.fill(readValues, 0f);
        Arrays.fill(writeValues, 0f);
    }

    /**
     * Swaps read/write buffers after an update pass.
     *
     * <p>The solver writes into {@code writeValues}, then this method makes those new
     * values become the next frame's {@code readValues}.</p>
     */
    public void swapBuffers() {
        float[] temp = readValues;
        readValues = writeValues;
        writeValues = temp;
    }
}
