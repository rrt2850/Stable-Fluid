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

    public ScalarField(int totalCellCount) {
        readValues = new float[totalCellCount];
        writeValues = new float[totalCellCount];
    }

    public void clear() {
        Arrays.fill(readValues, 0f);
        Arrays.fill(writeValues, 0f);
    }

    public void swapBuffers() {
        float[] temp = readValues;
        readValues = writeValues;
        writeValues = temp;
    }
}
