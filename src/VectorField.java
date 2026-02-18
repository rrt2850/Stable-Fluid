import java.util.Arrays;

/**
 * Represents a 2D velocity field defined on the fluid grid.
 *
 * Each cell stores a horizontal (x) and vertical (y) velocity component.
 * There's a double buffer for time step stuff
 */
public class VectorField {

    public float[] readVelocityX;
    public float[] readVelocityY;
    public float[] writeVelocityX;
    public float[] writeVelocityY;

    /**
     * Creates a velocity field with separate read/write buffers for X and Y velocity.
     *
     * @param totalCellCount number of active + ghost cells in the grid
     */
    public VectorField(int totalCellCount) {
        readVelocityX = new float[totalCellCount];
        readVelocityY = new float[totalCellCount];
        writeVelocityX = new float[totalCellCount];
        writeVelocityY = new float[totalCellCount];
    }

    /**
     * Clears both velocity components in both buffers.
     */
    public void clear() {
        Arrays.fill(readVelocityX, 0f);
        Arrays.fill(readVelocityY, 0f);
        Arrays.fill(writeVelocityX, 0f);
        Arrays.fill(writeVelocityY, 0f);
    }

    /**
     * Swaps read/write buffers for both velocity components.
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
