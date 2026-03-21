import java.io.File;
import java.util.Random;

public class Utils {
    /**
     * Generates a uniformly distributed random float in a closed-open interval.
     *
     * @param random random number generator
     * @param min inclusive lower bound
     * @param max exclusive upper bound
     * @return random value in the range {@code [min, max)}
     */
    public static float randomRange(Random random, float min, float max) {
        return min + random.nextFloat() * (max - min);
    }

    /**
     * Normalizes an angle to the range {@code [0, 360)} degrees.
     *
     * @param a angle in degrees
     * @return normalized angle in degrees
     */
    public static float normalizeAngle(float a) {
        a %= 360f;
        return a < 0 ? a + 360f : a;
    }

    /**
     * Tests whether two angles are within a given angular tolerance.
     *
     * <p>Comparison is performed on a circular domain.</p>
     *
     * @param a first angle in degrees
     * @param target target angle in degrees
     * @param eps tolerance in degrees
     * @return {@code true} if the angular difference is less than {@code eps}
     */
    public static boolean near(float a, float target, float eps) {
        float d = Math.abs((a - target + 180f) % 360f - 180f);
        return d < eps;
    }

    /**
     * Recursively deletes a directory and all files and subdirectories it contains.
     *
     * @param directory root directory to delete
     */
    public static void deleteDirectoryRecursively(File directory) {
        File[] files = directory.listFiles();
        if (files != null) {
            for (File file : files) {
                if (file.isDirectory()) {
                    deleteDirectoryRecursively(file);
                } else {
                    file.delete();
                }
            }
        }
        directory.delete();
    }
}
