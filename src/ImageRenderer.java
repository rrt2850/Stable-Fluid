import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import static java.lang.Math.clamp;

public class ImageRenderer {
    /**
     * Converts simulation density fields into an RGBA image.
     *
     * <p>Each color channel (red, green, blue) is sampled from the solver's
     * density fields and scaled by a shared global maximum.
     * This preserves hue relationships and avoids artificial white/yellow shifts.</p>
     *
     * <p>The image is generated at an arbitrary output resolution using
     * bilinear interpolation over the simulation grid.</p>
     *
     * @param grid simulation grid defining valid cell indices
     * @param solver solver providing density fields
     * @param outputWidth desired output image width in pixels
     * @param outputHeight desired output image height in pixels
     * @return a newly allocated {@link BufferedImage} containing the rendered density
     */

    public static BufferedImage createDensityImage(
            FluidGrid grid,
            FluidSolver solver,
            int outputWidth,
            int outputHeight
    ) {
        BufferedImage image = new BufferedImage(outputWidth, outputHeight, BufferedImage.TYPE_INT_ARGB);

        float[] red = solver.redDensityField.readValues;
        float[] green = solver.greenDensityField.readValues;
        float[] blue = solver.blueDensityField.readValues;

        // Normalize all channels by one shared peak density.
        // Per-channel normalization shifts hue and can create white/yellow artifacts.
        float maxDensity = 0.0f;
        for (int y = 1; y <= grid.height; y++) {
            for (int x = 1; x <= grid.width; x++) {
                int index = grid.index(x, y);
                maxDensity = Math.max(maxDensity, red[index]);
                maxDensity = Math.max(maxDensity, green[index]);
                maxDensity = Math.max(maxDensity, blue[index]);
            }
        }

        float densityNormalization = maxDensity > 0.0f ? maxDensity : 1.0f;

        for (int outY = 0; outY < outputHeight; outY++) {
            float simY = 1.0f + (outY / (float) outputHeight) * (grid.height - 1);

            for (int outX = 0; outX < outputWidth; outX++) {
                float simX = 1.0f + (outX / (float) outputWidth) * (grid.width - 1);

                float normalizedRed = clamp(bilinearSample(grid, red, simX, simY) / densityNormalization, 0.0f, 1.0f);
                float normalizedGreen = clamp(bilinearSample(grid, green, simX, simY) / densityNormalization, 0.0f, 1.0f);
                float normalizedBlue = clamp(bilinearSample(grid, blue, simX, simY) / densityNormalization, 0.0f, 1.0f);

                int r = Math.round(normalizedRed * 255.0f);
                int g = Math.round(normalizedGreen * 255.0f);
                int b = Math.round(normalizedBlue * 255.0f);

                int argb = (255 << 24) | (r << 16) | (g << 8) | b;
                image.setRGB(outX, outY, argb);
            }
        }

        return image;
    }

    /**
     * Samples a scalar field at a fractional grid position using bilinear interpolation.
     *
     * <p>Coordinates are clamped to the valid grid domain before sampling.</p>
     *
     * @param grid simulation grid used for index mapping
     * @param values scalar field values stored in grid-indexed layout
     * @param x continuous x-coordinate in grid space
     * @param y continuous y-coordinate in grid space
     * @return interpolated scalar value at the given position
     */
    public static float bilinearSample(FluidGrid grid, float[] values, float x, float y) {
        float clampedX = clamp(x, 1.0f, grid.width);
        float clampedY = clamp(y, 1.0f, grid.height);

        int x0 = (int) Math.floor(clampedX);
        int y0 = (int) Math.floor(clampedY);
        int x1 = Math.min(grid.width, x0 + 1);
        int y1 = Math.min(grid.height, y0 + 1);

        float tx = clampedX - x0;
        float ty = clampedY - y0;

        float v00 = values[grid.index(x0, y0)];
        float v10 = values[grid.index(x1, y0)];
        float v01 = values[grid.index(x0, y1)];
        float v11 = values[grid.index(x1, y1)];

        float top = v00 + tx * (v10 - v00);
        float bottom = v01 + tx * (v11 - v01);
        return top + ty * (bottom - top);
    }

    /**
     * Writes a {@link BufferedImage} to disk as a PNG file.
     *
     * @param image the image to write
     * @param outputPath destination file path
     * @throws RuntimeException if the image cannot be written
     */
    public static void saveImage(BufferedImage image, String outputPath) {
        try {
            ImageIO.write(image, "png", new File(outputPath));
        } catch (IOException exception) {
            throw new RuntimeException("Failed to write PNG to " + outputPath, exception);
        }
    }

    /**
     * Renders the current fluid density fields into a high-resolution PNG image.
     *
     * <p>The simulation grid is upscaled using bilinear interpolation to a
     * fixed output resolution defined by {@code FINAL_STILL_WIDTH} and
     * {@code FINAL_STILL_HEIGHT}.</p>
     *
     * @param grid the simulation grid defining logical cell layout
     * @param solver solver containing the current density fields
     * @param outputPath file path where the PNG image will be written
     */
    public static void saveDensityToPng(FluidGrid grid, FluidSolver solver, String outputPath) {
        // Upscaled final still
        BufferedImage image = createDensityImage(grid, solver, Constants.FINAL_STILL_WIDTH, Constants.FINAL_STILL_HEIGHT);
        saveImage(image, outputPath);
    }
}
