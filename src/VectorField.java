/**
 * Represents a 2D velocity field defined on the fluid grid.
 *
 * Each cell stores a horizontal (x) and vertical (y) velocity component.
 * There's a double buffer for time step stuff
 */
public class VectorField {

    /**
     * Horizontal (x-axis) velocity for each grid cell
     *
     * This buffer always contains the authoritative velocity values
     * that solver steps are allowed to read from
     */
    public float[] readVelocityX;

    /**
     * Vertical (y-axis) velocity for each grid cell
     *
     * This buffer always contains the authoritative velocity values
     * that solver steps are allowed to read from
     */
    public float[] readVelocityY;

    /**
     * Horizontal (x-axis) velocity write buffer
     *
     * Solver steps write newly computed velocity values into this buffer
     * Its contents are not considered valid until buffers are swapped
     */
    public float[] writeVelocityX;

    /**
     * Vertical (y-axis) velocity write buffer.
     *
     * Solver steps write newly computed velocity values into this buffer.
     * Its contents are not considered valid until buffers are swapped.
     */
    public float[] writeVelocityY;

    /**
     * Allocates velocity arrays sized to the grid
     *
     * @param totalCellCount total number of grid cells including ghost cells
     */
    public VectorField(int totalCellCount) {
        readVelocityX = new float[totalCellCount];
        readVelocityY = new float[totalCellCount];
        writeVelocityX = new float[totalCellCount];
        writeVelocityY = new float[totalCellCount];
    }

    /**
     * Sets all velocities to 0
     * Resets the simulation, clearing all values
     */
    public void clear() {
        for (int i = 0; i < readVelocityX.length; i++) {
            readVelocityX[i] = 0f;
            readVelocityY[i] = 0f;
        }
    }

    /**
     * Swaps the read and write velocity buffers.
     *
     * After calling this method:
     * - The newly computed velocities become the authoritative state
     * - The old read buffers become available for reuse as write buffers
     *
     * You may think "let's just copy it over" but that actually adds
     * more computation and time complexity
     */
    public void swapBuffers() {
        float[] tempX = readVelocityX;
        readVelocityX = writeVelocityX;
        writeVelocityX = tempX;

        float[] tempY = readVelocityY;
        readVelocityY = writeVelocityY;
        writeVelocityY = tempY;
    }
}
