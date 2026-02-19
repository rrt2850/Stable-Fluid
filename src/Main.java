import javax.imageio.ImageIO;
import java.awt.Toolkit;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Entry point for the fluid simulation
 *
 * <p>This class wires together simulation setup, per-step execution, and export helpers
 * for image/video output. Most methods are intentionally small utility methods so readers
 * can follow the pipeline from command-line arguments to rendered output.</p>
 */
public class Main {
    // Fluid properties
    private static final float DEFAULT_DENSITY_RATE = 1.2f;
    private static final float MIN_EMISSION_SPEED = 0.7f;
    private static final float MAX_EMISSION_SPEED = 1.1f;
    private static final float TIMESTEP = 0.020f;
    private static final float VISCOSITY = 0.0000001f;
    private static final float DIFFUSION_RATE = 0.00001f;
    private static final int SOLVER_ITERATIONS = 25;
    private static final float VORTICITY_CONFINEMENT = 10.0f;

    private static final int DEFAULT_SIMULATION_STEPS = 100;
    private static final int DEFAULT_EMITTER_COUNT = 12;

    // Defaults for PowerShell runner
    private static final int DEFAULT_GRID_WIDTH = 800;
    private static final int DEFAULT_GRID_HEIGHT = 400;

    // Final still export resolution (upscaled from sim)
    private static final int FINAL_STILL_WIDTH = 2400;
    private static final int FINAL_STILL_HEIGHT = 1200;

    private static final int MP4_FRAMES_PER_SECOND = 30;
    private static final int INTERMITTENT_SNAPSHOT_INTERVAL = 25;

    private static final long RANDOM_SEED = System.currentTimeMillis();

    // Emitter settings
    private static final float EMITTER_RADIUS_RATIO = 0.012f;
    private static final int MIN_EMITTER_RADIUS = 8;
    private static final int MAX_EMITTER_RADIUS = 60;
    private static final float EMITTER_ANGLE_VARIATION_DEGREES = 60.0f;
    private static final float WALL_TANGENT_EXCLUSION = 10f;
    private static final float CORNER_EXCLUSION = 5f;

    // Good seeds :)
    // 1771389195668L
    // 1771400966868L
    // 1771402741333L
    // 1771453995487L

    // 12 namespace colors represented as RGB triplets in the [0, 1] range.
    private static final float[][] NAMESPACE_COLORS = {
            {0.0f, 0.69f, 0.94f},
            {0.0f, 0.78f, 0.58f},
            {0.55f, 0.80f, 0.26f},
            {0.98f, 0.82f, 0.20f},
            {1.0f, 0.63f, 0.0f},
            {0.95f, 0.36f, 0.20f},
            {0.87f, 0.22f, 0.52f},
            {0.67f, 0.27f, 0.87f},
            {0.31f, 0.42f, 0.94f},
            {0.18f, 0.68f, 0.98f},
            {0.0f, 0.73f, 0.75f},
            {0.42f, 0.75f, 0.45f}
    };

    /**
     * Program entry point
     *
     * <p>Initializes the simulation configuration, constructs the fluid grid,
     * solver, and emitters, then advances the simulation for a fixed number
     * of steps. Depending on command-line flags, this method may export
     * intermediate PNG frames, a final high-resolution still image, and an
     * MP4 video encoded via ffmpeg</p>
     *
     * <p>This method also reports basic diagnostics such as emitter placement,
     * center-cell density, divergence metrics (optional), and total runtime</p>
     *
     * @param args command-line arguments in the following order:
     *             <ol>
     *               <li>grid width (positive integer)</li>
     *               <li>grid height (positive integer)</li>
     *               <li>number of edge emitters (positive integer)</li>
     *               <li>number of simulation steps (positive integer)</li>
     *               <li>whether to export an MP4 video ("true" or "false")</li>
     *               <li>whether to log diagnostics every step ("true" or "false")</li>
     *             </ol>
     */
    public static void main(String[] args) {
        long programStartTime = System.nanoTime();
        SimulationConfig config = parseConfig(args);

        // Makes a new fluid grid
        // TODO: Figure out a way to get cellSize that's not biased towards square resolutions
        FluidGrid grid = new FluidGrid(
                config.gridWidth,
                config.gridHeight,
                1.0f / Math.max(config.gridWidth, config.gridHeight)
        );

        SimulationParameters parameters = new SimulationParameters(
                TIMESTEP,
                VISCOSITY,
                DIFFUSION_RATE,
                SOLVER_ITERATIONS,
                VORTICITY_CONFINEMENT
        );

        List<FluidSource> sources = List.of();

        Random random = new Random(RANDOM_SEED);
        List<FluidEmitter> emitters = generateEdgeEmitters(grid, config.emitterCount, random);

        List<RadialFluidEmitter> radialEmitters = List.of(
                /*
                new RadialFluidEmitter(
                        (grid.width + 1) / 2,
                        (grid.height + 1) / 2,
                        10,
                        4.0f,
                        0.08f,
                        0.9f,
                        0.3f,
                        0.9f
                )*/
        );

        List<Vortex> vortexes = List.of(
                /*
                new Vortex(
                        (grid.width + 1) / 4,
                        (grid.height + 1) / 2,
                        200,
                        4000.0f,
                        20.0f,
                        1000f
                ),
                new Vortex(
                        ((grid.width + 1) / 4) * 3,
                        (grid.height + 1) / 2,
                        200,
                        4000.0f,
                        20.0f,
                        1000f
                )*/
        );

        FluidSolver solver = new FluidSolver(grid, parameters, sources, emitters, radialEmitters, vortexes);

        System.out.println("Seed: " + RANDOM_SEED);
        System.out.println("Generated " + emitters.size() + " edge emitters:");

        // Go into each emitter and convert the rgb values to 255 instead of 0-1
        for (int i = 0; i < emitters.size(); i++) {
            FluidEmitter emitter = emitters.get(i);
            int red = Math.round(clamp(emitter.red(), 0.0f, 1.0f) * 255.0f);
            int green = Math.round(clamp(emitter.green(), 0.0f, 1.0f) * 255.0f);
            int blue = Math.round(clamp(emitter.blue(), 0.0f, 1.0f) * 255.0f);
            System.out.printf(
                    "  #%d at (%d,%d) angle=%.2f deg speed=%.2f color=(%d,%d,%d)%n",
                    i + 1,
                    emitter.gridX(),
                    emitter.gridY(),
                    emitter.angleDegrees(),
                    emitter.emissionSpeed(),
                    red,
                    green,
                    blue
            );
        }

        System.out.println("Configured " + radialEmitters.size() + " radial emitters and "
                + vortexes.size() + " vortex emitters.");

        boolean takeIntermittentSnapshots = config.simulationSteps >= INTERMITTENT_SNAPSHOT_INTERVAL;
        Path tempFramesDirectory = null;
        try {
            if (config.exportVideo) {
                tempFramesDirectory = Files.createTempDirectory("stable-fluid-frames-");
            }

            for (int step = 1; step <= config.simulationSteps; step++) {
                solver.step();

                // ---------------------------------------------------------------------
                // MP4 video frames (simulation resolution, fast)
                // ---------------------------------------------------------------------
                if (config.exportVideo) {
                    BufferedImage videoFrame = createDensityImage(
                            grid,
                            solver,
                            grid.width,
                            grid.height
                    );

                    Path framePath = tempFramesDirectory.resolve(
                            String.format("frame-%05d.png", step - 1)
                    );
                    saveImage(videoFrame, framePath.toString());
                }

                // ---------------------------------------------------------------------
                // Intermittent snapshots (UPSCALED, high quality)
                // ---------------------------------------------------------------------
                if (takeIntermittentSnapshots && step % INTERMITTENT_SNAPSHOT_INTERVAL == 0) {
                    BufferedImage snapshotImage = createDensityImage(
                            grid,
                            solver,
                            FINAL_STILL_WIDTH,
                            FINAL_STILL_HEIGHT
                    );

                    String intermittentPath =
                            String.format("./results/density-step-%05d.png", step);

                    saveImage(snapshotImage, intermittentPath);
                    System.out.println("Saved intermittent density image to " + intermittentPath);
                }

                if (config.logEveryStep) System.out.println("Step " + step);
            }

            int center = grid.index((grid.width + 1) / 2, (grid.height + 1) / 2);
            float centerDensity = solver.redDensityField.readValues[center]
                    + solver.greenDensityField.readValues[center]
                    + solver.blueDensityField.readValues[center];
            System.out.println("Simulation ran. Center density=" + centerDensity);

            String outputPath = "density.png";
            saveDensityToPng(grid, solver, outputPath);
            System.out.println("Saved density image to " + outputPath);

            if (config.exportVideo) {
                String videoOutputPath = "./results/density-diffusion.mp4";
                boolean mp4Exported = saveDensityToMp4(tempFramesDirectory, videoOutputPath);
                if (mp4Exported) {
                    System.out.println("Saved density animation to " + videoOutputPath);
                } else {
                    System.out.println("Skipped MP4 export because ffmpeg is not available on PATH.");
                }
            } else {
                System.out.println("Skipped MP4 export. Pass a 5th argument of 'true' to enable video output.");
            }
        } catch (IOException exception) {
            throw new RuntimeException("Failed to create temporary frame directory.", exception);
        } finally {
            if (tempFramesDirectory != null) {
                deleteDirectoryRecursively(tempFramesDirectory.toFile());
            }
            Toolkit.getDefaultToolkit().beep();
        }

        long programEndTime = System.nanoTime();
        double elapsedSeconds = (programEndTime - programStartTime) / 1_000_000_000.0;
        System.out.printf("Total runtime: %.3f seconds%n", elapsedSeconds);
    }

    /**
     * Parses command-line arguments into an immutable simulation configuration.
     *
     * <p>Each argument is optional. If an argument is missing, a predefined
     * default value is used. Invalid values (non-integer or non-positive)
     * result in an {@link IllegalArgumentException}.</p>
     *
     * @param args raw command-line arguments passed to {@link #main(String[])}
     * @return a {@link SimulationConfig} containing all resolved configuration values
     */
    private static SimulationConfig parseConfig(String[] args) {
        // Defaults match your PowerShell runner, but you can still override via args.
        int gridWidth = parsePositiveInt(args, 0, DEFAULT_GRID_WIDTH, "grid width");
        int gridHeight = parsePositiveInt(args, 1, DEFAULT_GRID_HEIGHT, "grid height");
        int emitterCount = parsePositiveInt(args, 2, DEFAULT_EMITTER_COUNT, "emitter count");
        int simulationSteps = parsePositiveInt(args, 3, DEFAULT_SIMULATION_STEPS, "simulation steps");
        boolean exportVideo = parseBoolean(args, 4, false, "export video");
        boolean logEveryStep = parseBoolean(args, 5, false, "log every step");

        return new SimulationConfig(gridWidth, gridHeight, emitterCount, simulationSteps, exportVideo, logEveryStep);
    }

    /**
     * Parses a strictly positive integer from a given argument index.
     *
     * <p>If the argument is missing, the provided default value is returned.
     * If the argument is present but cannot be parsed as a positive integer,
     * an {@link IllegalArgumentException} is thrown.</p>
     *
     * @param args raw command-line argument array
     * @param index index into {@code args} to read
     * @param defaultValue value used when {@code args.length <= index}
     * @param argumentName human-readable name used in error messages
     * @return the parsed positive integer or {@code defaultValue}
     */
    private static int parsePositiveInt(String[] args, int index, int defaultValue, String argumentName) {
        if (args.length <= index) {
            return defaultValue;
        }

        try {
            int parsed = Integer.parseInt(args[index]);
            if (parsed <= 0) {
                throw new IllegalArgumentException(argumentName + " must be greater than zero.");
            }
            return parsed;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(argumentName + " must be an integer.", exception);
        }
    }

    /**
     * Parses a strict boolean value from a command-line argument.
     *
     * <p>Only the exact strings {@code "true"} or {@code "false"} (case-insensitive)
     * are accepted. Any other value results in an exception.</p>
     *
     * @param args raw command-line argument array
     * @param index index into {@code args} to read
     * @param defaultValue value used when {@code args.length <= index}
     * @param argumentName human-readable name used in error messages
     * @return the parsed boolean or {@code defaultValue}
     */
    private static boolean parseBoolean(String[] args, int index, boolean defaultValue, String argumentName) {
        if (args.length <= index) {
            return defaultValue;
        }

        String parsed = args[index].trim().toLowerCase();
        if ("true".equals(parsed)) {
            return true;
        }

        if ("false".equals(parsed)) {
            return false;
        }

        throw new IllegalArgumentException(argumentName + " must be either 'true' or 'false'.");
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

    private static void saveDensityToPng(FluidGrid grid, FluidSolver solver, String outputPath) {
        // Upscaled final still
        BufferedImage image = createDensityImage(grid, solver, FINAL_STILL_WIDTH, FINAL_STILL_HEIGHT);
        saveImage(image, outputPath);
    }

    /**
     * Writes a {@link BufferedImage} to disk as a PNG file.
     *
     * @param image the image to write
     * @param outputPath destination file path
     * @throws RuntimeException if the image cannot be written
     */
    private static void saveImage(BufferedImage image, String outputPath) {
        try {
            ImageIO.write(image, "png", new File(outputPath));
        } catch (IOException exception) {
            throw new RuntimeException("Failed to write PNG to " + outputPath, exception);
        }
    }

    /**
     * Converts simulation density fields into an RGBA image.
     *
     * <p>Each color channel (red, green, blue) is sampled independently from
     * the solver's density fields and normalized by its own maximum value.
     * This prevents strong channels from washing out weaker ones.</p>
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

    private static BufferedImage createDensityImage(
            FluidGrid grid,
            FluidSolver solver,
            int outputWidth,
            int outputHeight
    ) {
        BufferedImage image = new BufferedImage(outputWidth, outputHeight, BufferedImage.TYPE_INT_ARGB);

        float[] red = solver.redDensityField.readValues;
        float[] green = solver.greenDensityField.readValues;
        float[] blue = solver.blueDensityField.readValues;

        // Normalize each channel by its own max value to avoid washing out weaker colors.
        float maxRedDensity = 0.0f;
        float maxGreenDensity = 0.0f;
        float maxBlueDensity = 0.0f;
        for (int y = 1; y <= grid.height; y++) {
            for (int x = 1; x <= grid.width; x++) {
                int index = grid.index(x, y);
                maxRedDensity = Math.max(maxRedDensity, red[index]);
                maxGreenDensity = Math.max(maxGreenDensity, green[index]);
                maxBlueDensity = Math.max(maxBlueDensity, blue[index]);
            }
        }

        float redNormalization = maxRedDensity > 0.0f ? maxRedDensity : 1.0f;
        float greenNormalization = maxGreenDensity > 0.0f ? maxGreenDensity : 1.0f;
        float blueNormalization = maxBlueDensity > 0.0f ? maxBlueDensity : 1.0f;

        for (int outY = 0; outY < outputHeight; outY++) {
            float simY = 1.0f + (outY / (float) outputHeight) * (grid.height - 1);

            for (int outX = 0; outX < outputWidth; outX++) {
                float simX = 1.0f + (outX / (float) outputWidth) * (grid.width - 1);

                float normalizedRed = clamp(bilinearSample(grid, red, simX, simY) / redNormalization, 0.0f, 1.0f);
                float normalizedGreen = clamp(bilinearSample(grid, green, simX, simY) / greenNormalization, 0.0f, 1.0f);
                float normalizedBlue = clamp(bilinearSample(grid, blue, simX, simY) / blueNormalization, 0.0f, 1.0f);

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
    private static float bilinearSample(FluidGrid grid, float[] values, float x, float y) {
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
     * Encodes a sequence of numbered PNG frames into an MP4 video using ffmpeg.
     *
     * <p>Frames must be named {@code frame-00000.png}, {@code frame-00001.png}, etc.
     * The resulting MP4 is encoded using H.264 with a YUV420 pixel format.</p>
     *
     * @param framesDirectory directory containing PNG frame images
     * @param outputPath destination path for the MP4 file
     * @return {@code true} if encoding succeeded, {@code false} if ffmpeg is unavailable
     * @throws RuntimeException if ffmpeg fails or is interrupted
     */
    private static boolean saveDensityToMp4(Path framesDirectory, String outputPath) {
        if (framesDirectory == null) {
            throw new IllegalArgumentException("framesDirectory must not be null.");
        }

        try {
            Path tempOutput = framesDirectory.resolve("density-diffusion.mp4");
            ProcessBuilder ffmpegBuilder = new ProcessBuilder(
                    "ffmpeg",
                    "-y",
                    "-framerate", String.valueOf(MP4_FRAMES_PER_SECOND),
                    "-i", framesDirectory.resolve("frame-%05d.png").toString(),
                    "-c:v", "libx264",
                    "-pix_fmt", "yuv420p",
                    tempOutput.toString()
            );
            ffmpegBuilder.inheritIO();
            Process ffmpegProcess = ffmpegBuilder.start();
            int exitCode = ffmpegProcess.waitFor();
            if (exitCode != 0) {
                throw new RuntimeException("ffmpeg exited with code " + exitCode + ". Is ffmpeg installed and available on PATH?");
            }

            Files.move(tempOutput, Path.of(outputPath), StandardCopyOption.REPLACE_EXISTING);
            return true;
        } catch (IOException exception) {
            if (exception.getMessage() != null && exception.getMessage().contains("Cannot run program \"ffmpeg\"")) {
                return false;
            }
            throw new RuntimeException("Failed to write MP4 to " + outputPath + ". Is ffmpeg installed and available on PATH?", exception);
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("MP4 generation interrupted.", exception);
        }
    }

    /**
     * Recursively deletes a directory and all files and subdirectories it contains.
     *
     * @param directory root directory to delete
     */
    private static void deleteDirectoryRecursively(File directory) {
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

    /**
     * Generates fluid emitters positioned near the grid edges.
     *
     * <p>Emitters are placed along the four boundaries and aimed roughly toward
     * the grid center, with randomized angular variation, speed, and color.
     * Additional constraints prevent near-tangential wall flow and corner-sniping.</p>
     *
     * @param grid simulation grid used for placement bounds
     * @param emitterCount number of emitters to generate
     * @param random random source used for placement and parameter jitter
     * @return list of configured {@link FluidEmitter} instances
     */
    private static List<FluidEmitter> generateEdgeEmitters(FluidGrid grid, int emitterCount, Random random) {
        List<FluidEmitter> emitters = new ArrayList<>();

        float centerX = (grid.width + 1) / 2.0f;
        float centerY = (grid.height + 1) / 2.0f;

        int radius = computeEmitterRadius(grid);
        int offset = radius + 2;

        while (emitters.size() < emitterCount) {
            int side = random.nextInt(4);
            int x;
            int y;

            switch (side) {
                case 0 -> { // top
                    x = 1 + random.nextInt(grid.width);
                    y = 1 + offset;
                }
                case 1 -> { // bottom
                    x = 1 + random.nextInt(grid.width);
                    y = grid.height - offset;
                }
                case 2 -> { // left
                    x = 1 + offset;
                    y = 1 + random.nextInt(grid.height);
                }
                case 3 -> { // right
                    x = grid.width - offset;
                    y = 1 + random.nextInt(grid.height);
                }
                default -> throw new IllegalStateException();
            }

            float dx = centerX - x;
            float dy = centerY - y;
            float centerAngle = (float) Math.toDegrees(Math.atan2(dy, dx));

            float candidateAngle;
            int attempts = 0;

            do {
                float jitter = randomRange(
                        random,
                        -EMITTER_ANGLE_VARIATION_DEGREES,
                        EMITTER_ANGLE_VARIATION_DEGREES
                );
                candidateAngle = centerAngle + jitter;
                attempts++;
            } while (
                    (!isValidForWall(candidateAngle, side)) &&
                            attempts < 25
            );

            float[] color = NAMESPACE_COLORS[
                    emitters.size() % NAMESPACE_COLORS.length
                    ];

            float speed = randomRange(
                    random,
                    MIN_EMISSION_SPEED,
                    MAX_EMISSION_SPEED
            );

            emitters.add(new FluidEmitter(
                    x,
                    y,
                    radius,
                    DEFAULT_DENSITY_RATE,
                    candidateAngle,
                    speed,
                    color[0],
                    color[1],
                    color[2]
            ));
        }

        return emitters;
    }

    /**
     * Generates a uniformly distributed random float in a closed-open interval.
                *
                * @param random random number generator
     * @param min inclusive lower bound
     * @param max exclusive upper bound
     * @return random value in the range {@code [min, max)}
     */
    private static float randomRange(Random random, float min, float max) {
        return min + random.nextFloat() * (max - min);
    }

    /**
     * Computes an emitter radius based on grid dimensions.
     *
     * <p>The radius scales with grid size but is clamped to a practical
     * minimum and maximum.</p>
     *
     * @param grid simulation grid
     * @return emitter radius in grid cells
     */
    private static int computeEmitterRadius(FluidGrid grid) {
        int radiusFromRatio = Math.round(Math.min(grid.width, grid.height) * EMITTER_RADIUS_RATIO);
        return Math.min(MAX_EMITTER_RADIUS, Math.max(MIN_EMITTER_RADIUS, radiusFromRatio));
    }

    /**
     * Clamps a floating-point value to a closed interval.
     *
     * @param value input value
     * @param min lower bound
     * @param max upper bound
     * @return {@code value} constrained to {@code [min, max]}
     */
    private static float clamp(float value, float min, float max) {
        return Math.max(min, Math.min(max, value));
    }


    private record SimulationConfig(
            int gridWidth,
            int gridHeight,
            int emitterCount,
            int simulationSteps,
            boolean exportVideo,
            boolean logEveryStep
    ) {}

    /**
     * Determines whether an emission angle is valid for a given wall.
     *
     * <p>This method rejects angles that are nearly tangential to the wall
     * or that point directly toward grid corners, which can cause numerical
     * artifacts or degenerate flow patterns.</p>
     *
     * @param angleDeg emission angle in degrees
     * @param side wall identifier (0=top, 1=bottom, 2=left, 3=right)
     * @return {@code true} if the angle is acceptable for this wall
     */
    private static boolean isValidForWall(float angleDeg, int side) {
        float a = normalizeAngle(angleDeg);

        // Reject near-tangential flow along the wall
        float tangent = switch (side) {
            case 0, 1 -> 0f;    // top/bottom → horizontal tangent
            case 2, 3 -> 90f;   // left/right → vertical tangent
            default -> throw new IllegalArgumentException();
        };

        if (near(a, tangent, WALL_TANGENT_EXCLUSION) ||
                near(a, tangent + 180f, WALL_TANGENT_EXCLUSION)) {
            return false;
        }

        // Reject corner-sniping diagonals
        return switch (side) {
            case 0 -> !near(a, 45f, CORNER_EXCLUSION) &&
                    !near(a, 135f, CORNER_EXCLUSION);
            case 1 -> !near(a, 225f, CORNER_EXCLUSION) &&
                    !near(a, 315f, CORNER_EXCLUSION);
            case 2 -> !near(a, 45f, CORNER_EXCLUSION) &&
                    !near(a, 315f, CORNER_EXCLUSION);
            case 3 -> !near(a, 135f, CORNER_EXCLUSION) &&
                    !near(a, 225f, CORNER_EXCLUSION);
            default -> true;
        };
    }

    /**
     * Normalizes an angle to the range {@code [0, 360)} degrees.
     *
     * @param a angle in degrees
     * @return normalized angle in degrees
     */
    private static float normalizeAngle(float a) {
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
    private static boolean near(float a, float target, float eps) {
        float d = Math.abs((a - target + 180f) % 360f - 180f);
        return d < eps;
    }

}
