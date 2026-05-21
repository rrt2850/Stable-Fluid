
import java.awt.GraphicsEnvironment;
import java.awt.Toolkit;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.List;
import java.util.Random;
import static java.lang.Math.clamp;

/**
 * Entry point for the fluid simulation
 *
 * <p>This class wires together simulation setup, per-step execution, and export helpers
 * for image/video output. Most methods are intentionally small utility methods so readers
 * can follow the pipeline from command-line arguments to rendered output.</p>
 */
public class Main {
    /**
     * Program entry point
     *
     * <p>Initializes the simulation configuration, constructs the fluid grid,
     * solver, and emitters, then advances the simulation for a fixed number
     * of steps. Depending on commandline flags, this method may export
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
        FluidGrid grid = new FluidGrid(
                config.gridWidth,
                config.gridHeight,
                1.0f / Math.max(config.gridWidth, config.gridHeight)
        );

        SimulationParameters parameters = new SimulationParameters(
                Constants.TIMESTEP,
                Constants.VISCOSITY,
                Constants.DIFFUSION_RATE,
                Constants.SOLVER_ITERATIONS,
                Constants.VORTICITY_CONFINEMENT
        );

        // Change the commented out function to change modes :)
        //EmitterHolder holder = InitializeRegionDivision(grid);
        EmitterHolder holder = InitializeEdgeShooterDemo(grid, config);

        List<FluidSource> sources = List.of();
        List<FluidEmitter> emitters = holder.fluidEmitters();
        List<RadialFluidEmitter> radialEmitters = holder.radialFluidEmitters();
        List<Vortex> vortexes = holder.vortexes();
        List<Wall> walls = holder.walls();


        FluidSolver solver = new FluidSolver(grid, parameters, sources, emitters, radialEmitters, vortexes, walls);

        System.out.println("Seed: " + Constants.RANDOM_SEED);
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

        boolean takeIntermittentSnapshots = config.simulationSteps >= Constants.SNAPSHOT_INTERVAL;
        Path tempFramesDirectory = null;
        SimulationPreviewWindow previewWindow = null;
        try {
            if (!GraphicsEnvironment.isHeadless()) {
                previewWindow = new SimulationPreviewWindow(grid.width, grid.height);
            } else {
                System.out.println("Skipped live preview window because the environment is headless.");
            }

            if (config.exportVideo) {
                tempFramesDirectory = Files.createTempDirectory("stable-fluid-frames-");
            }

            for (int step = 1; step <= config.simulationSteps; step++) {
                solver.step();

                // Generate a small frame for the mp4 if exportVideo is enabled
                BufferedImage simulationResolutionFrame = null;
                if (config.exportVideo || previewWindow != null) {
                    simulationResolutionFrame = ImageRenderer.createDensityImage(
                            grid,
                            solver,
                            grid.width,
                            grid.height
                    );
                }

                // if exportVideo is enabled, save the current frame in the temporary directory
                if (config.exportVideo) {
                    Path framePath = tempFramesDirectory.resolve(
                            String.format("frame-%05d.png", step - 1)
                    );
                    ImageRenderer.saveImage(simulationResolutionFrame, framePath.toString());
                }

                // if the preview window exists, update it with the current frame
                if (previewWindow != null) {
                    previewWindow.update(simulationResolutionFrame, step);
                }

                // Upscale the current frame every SNAPSHOT_INTERVAL frames
                if (takeIntermittentSnapshots && step % Constants.SNAPSHOT_INTERVAL == 0) {
                    BufferedImage snapshotImage = ImageRenderer.createDensityImage(
                            grid,
                            solver,
                            Constants.FINAL_STILL_WIDTH,
                            Constants.FINAL_STILL_HEIGHT
                    );

                    String intermittentPath =
                            String.format("./results/density-step-%05d.png", step);

                    ImageRenderer.saveImage(snapshotImage, intermittentPath);
                    System.out.println("Saved intermittent density image to " + intermittentPath);
                }

                if (config.logEveryStep) System.out.println("Step " + step);
            }

            String outputPath = "density.png";
            ImageRenderer.saveDensityToPng(grid, solver, outputPath);
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
            if (previewWindow != null) {
                previewWindow.close();
            }
            if (tempFramesDirectory != null) {
                Utils.deleteDirectoryRecursively(tempFramesDirectory.toFile());
            }
            Toolkit.getDefaultToolkit().beep();
        }

        long programEndTime = System.nanoTime();
        double elapsedSeconds = (programEndTime - programStartTime) / 1_000_000_000.0;
        System.out.printf("Total runtime: %.3f seconds%n", elapsedSeconds);
    }

    /**
     * Initializes a series of emitters around the edge of the frame
     *
     * @param grid The FluidGrid the emitters are in
     * @return an EmitterHolder containing the lists of emitters
     */
    private static EmitterHolder InitializeEdgeShooterDemo(FluidGrid grid, SimulationConfig config){
        Random random = new Random(Constants.RANDOM_SEED);
        List<FluidEmitter> emitters = EmitterMaker.generateEdgeEmitters(grid, config.emitterCount, random);

        return new EmitterHolder(emitters, List.of(), List.of(), List.of());
    }

    /**
     * Initializes the emitters in their correct positions for the Region Division assignment
     *
     * @param grid The FluidGrid the emitters are in
     * @return an EmitterHolder containing the lists of emitters
     */
    private static EmitterHolder InitializeRegionDivision(FluidGrid grid){
        List<FluidEmitter> fluidEmitters = List.of();
        List<RadialFluidEmitter> radialEmitters = List.of(
                new RadialFluidEmitter(
                        ((grid.width + 1) / 8) * 7,
                        (grid.height + 1) / 8,
                        20,
                        1.5f,
                        0.5f,
                        248f/255f, 248f/255f, 56f/255f
                ),
                new RadialFluidEmitter(
                        ((grid.width + 1) / 8),
                        ((grid.height + 1) / 8) * 7,
                        20,
                        1.5f,
                        0.5f,
                        221f / 255f, 74f / 255f, 225f/ 255f
                ),
                new RadialFluidEmitter(
                        ((grid.width + 1) / 8) * 7,
                        ((grid.height + 1) / 8) * 7,
                        20,
                        1.5f,
                        0.5f,
                        248f/255f, 248f/255f, 56f/255f
                ),
                new RadialFluidEmitter(
                        ((grid.width + 1) / 8),
                        (grid.height + 1) / 8,
                        20,
                        1.5f,
                        0.5f,
                        221f / 255f, 74f / 255f, 225f/ 255f
                ),
                new RadialFluidEmitter(
                        ((grid.width + 1) / 2),
                        ((grid.height + 1) / 8) * 7,
                        9,
                        0.6f,
                        0.3f,
                        3f / 255f, 25f / 255f, 19f/ 255f
                ),
                new RadialFluidEmitter(
                        ((grid.width + 1) / 2),
                        ((grid.height + 1) / 8),
                        9,
                        0.6f,
                        0.3f,
                        7f / 255f, 231f / 255f, 242f/ 255f
                )

        );

        List<Vortex> vortexes = List.of(

                new Vortex(
                        (grid.width + 1) / 2,
                        (grid.height + 1) / 2,
                        grid.width / 3,
                        0.15f,
                        0.00001f,
                        0.1f,
                        -1
                )/*,
                new Vortex(
                        ((grid.width + 1) / 4) * 3,
                        (grid.height + 1) / 2,
                        200,
                        4000.0f,
                        20.0f,
                        1000f
                )*/
        );

        int leftWallX = (grid.width + 1) / 4;
        int rightWallX = ((grid.width + 1) / 4) * 3;

        int rightGapTop = (grid.height / 10) * 4;
        int rightGapBottom = (grid.height/ 10) * 6;

        int leftGapTop = (grid.height / 10) * 4;
        int leftGapBottom = (grid.height / 10) * 6;

        List<Wall> walls = List.of(
                new Wall(
                        40,
                        new WallPoint(rightWallX, 1),
                        new WallPoint(rightWallX, leftGapTop)
                ),
                new Wall(
                        40,
                        new WallPoint(rightWallX, leftGapBottom),
                        new WallPoint(rightWallX, grid.height)
                ),
                new Wall(
                        40,
                        new WallPoint(leftWallX, 1),
                        new WallPoint(leftWallX, rightGapTop)
                ),
                new Wall(
                        40,
                        new WallPoint(leftWallX, rightGapBottom),
                        new WallPoint(leftWallX, grid.height)
                )
        );

        return new EmitterHolder(fluidEmitters, radialEmitters, vortexes, walls);
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
        int gridWidth = parsePositiveInt(args, 0, Constants.DEFAULT_GRID_WIDTH, "grid width");
        int gridHeight = parsePositiveInt(args, 1, Constants.DEFAULT_GRID_HEIGHT, "grid height");
        int emitterCount = parsePositiveInt(args, 2, Constants.DEFAULT_EMITTER_COUNT, "emitter count");
        int simulationSteps = parsePositiveInt(args, 3, Constants.DEFAULT_SIMULATION_STEPS, "simulation steps");
        boolean exportVideo = parseBoolean(args, 4, false, "export video");
        boolean logEveryStep = parseBoolean(args, 5, false, "log every step");

        return new SimulationConfig(
                gridWidth,
                gridHeight,
                emitterCount,
                simulationSteps,
                exportVideo,
                logEveryStep
        );
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
                    "-framerate", String.valueOf(Constants.MP4_FRAMES_PER_SECOND),
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


    private record SimulationConfig(
            int gridWidth,
            int gridHeight,
            int emitterCount,
            int simulationSteps,
            boolean exportVideo,
            boolean logEveryStep
    ) {}
}
