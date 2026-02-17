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

public class Main {
    private static final float DEFAULT_DENSITY_RATE = 15.0f;
    private static final float MIN_EMISSION_SPEED = 0.5f;
    private static final float MAX_EMISSION_SPEED = 0.75f;
    private static final float TIMESTEP = 0.020f;
    private static final float VISCOSITY = 0.00001f;
    private static final float DIFFUSION_RATE = 0.0001f;
    private static final int SOLVER_ITERATIONS = 25;
    private static final int EMITTER_RADIUS = 20; // Explodes around 35
    private static final float VORTICITY_CONFINEMENT = 2.1f;

    private static final int DEFAULT_SIMULATION_STEPS = 100;
    private static final int DEFAULT_EMITTER_COUNT = 12;
    private static final int MP4_FRAMES_PER_SECOND = 30;
    private static final int INTERMITTENT_SNAPSHOT_INTERVAL = 25;
    private static final int RANDOM_SEED = 67; // TODO: switch to time-based seed after testing


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

    public static void main(String[] args) {
        long programStartTime = System.nanoTime();
        SimulationConfig config = parseConfig(args);

        FluidGrid grid = new FluidGrid(config.gridWidth, config.gridHeight, 1.0f / Math.max(config.gridWidth, config.gridHeight));
        SimulationParameters parameters = new SimulationParameters(TIMESTEP, VISCOSITY, DIFFUSION_RATE, SOLVER_ITERATIONS, VORTICITY_CONFINEMENT);

        List<FluidSource> sources = List.of();

        Random random = new Random(RANDOM_SEED);
        List<FluidEmitter> emitters = generateEdgeEmitters(grid, config.emitterCount, random);


        FluidSolver solver = new FluidSolver(grid, parameters, sources, emitters);

        System.out.println("Generated " + emitters.size() + " edge emitters:");
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

        boolean takeIntermittentSnapshots = config.simulationSteps >= INTERMITTENT_SNAPSHOT_INTERVAL;
        Path tempFramesDirectory = null;
        try {
            if (config.exportVideo) {
                tempFramesDirectory = Files.createTempDirectory("stable-fluid-frames-");
            }

            for (int step = 1; step <= config.simulationSteps; step++) {
                solver.step();

                BufferedImage stepImage = null;
                if (config.exportVideo || (takeIntermittentSnapshots && step % INTERMITTENT_SNAPSHOT_INTERVAL == 0)) {
                    stepImage = createDensityImage(grid, solver);
                }

                if (config.exportVideo) {
                    Path framePath = tempFramesDirectory.resolve(String.format("frame-%05d.png", step - 1));
                    saveImage(stepImage, framePath.toString());
                }

                if (takeIntermittentSnapshots && step % INTERMITTENT_SNAPSHOT_INTERVAL == 0) {
                    String intermittentPath = String.format("./results/density-step-%05d.png", step);
                    saveImage(stepImage, intermittentPath);
                    System.out.println("Saved intermittent density image to " + intermittentPath);
                }

                if (config.logEveryStep) {
                    System.out.println("Step " + step + " divergence RMS=" + solver.computeVelocityDivergenceRms());
                }
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

    private static SimulationConfig parseConfig(String[] args) {
        int gridWidth = parsePositiveInt(args, 0, 128, "grid width");
        int gridHeight = parsePositiveInt(args, 1, 128, "grid height");
        int emitterCount = parsePositiveInt(args, 2, DEFAULT_EMITTER_COUNT, "emitter count");
        int simulationSteps = parsePositiveInt(args, 3, DEFAULT_SIMULATION_STEPS, "simulation steps");
        boolean exportVideo = parseBoolean(args, 4, false, "export video");
        boolean logEveryStep = parseBoolean(args, 5, false, "log every step");

        return new SimulationConfig(gridWidth, gridHeight, emitterCount, simulationSteps, exportVideo, logEveryStep);
    }

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

    private static void saveDensityToPng(FluidGrid grid, FluidSolver solver, String outputPath) {
        BufferedImage image = createDensityImage(grid, solver);
        saveImage(image, outputPath);
    }

    private static void saveImage(BufferedImage image, String outputPath) {
        try {
            ImageIO.write(image, "png", new File(outputPath));
        } catch (IOException exception) {
            throw new RuntimeException("Failed to write PNG to " + outputPath, exception);
        }
    }

    private static BufferedImage createDensityImage(FluidGrid grid, FluidSolver solver) {
        BufferedImage image = new BufferedImage(grid.width, grid.height, BufferedImage.TYPE_INT_ARGB);

        float maxRedDensity = 0.0f;
        float maxGreenDensity = 0.0f;
        float maxBlueDensity = 0.0f;
        for (int y = 1; y <= grid.height; y++) {
            for (int x = 1; x <= grid.width; x++) {
                int index = grid.index(x, y);
                maxRedDensity = Math.max(maxRedDensity, solver.redDensityField.readValues[index]);
                maxGreenDensity = Math.max(maxGreenDensity, solver.greenDensityField.readValues[index]);
                maxBlueDensity = Math.max(maxBlueDensity, solver.blueDensityField.readValues[index]);
            }
        }

        float redNormalization = maxRedDensity > 0.0f ? maxRedDensity : 1.0f;
        float greenNormalization = maxGreenDensity > 0.0f ? maxGreenDensity : 1.0f;
        float blueNormalization = maxBlueDensity > 0.0f ? maxBlueDensity : 1.0f;

        for (int y = 1; y <= grid.height; y++) {
            for (int x = 1; x <= grid.width; x++) {
                int index = grid.index(x, y);
                float normalizedRed = clamp(solver.redDensityField.readValues[index] / redNormalization, 0.0f, 1.0f);
                float normalizedGreen = clamp(solver.greenDensityField.readValues[index] / greenNormalization, 0.0f, 1.0f);
                float normalizedBlue = clamp(solver.blueDensityField.readValues[index] / blueNormalization, 0.0f, 1.0f);

                int red = Math.round(normalizedRed * 255.0f);
                int green = Math.round(normalizedGreen * 255.0f);
                int blue = Math.round(normalizedBlue * 255.0f);

                int argb = (255 << 24) | (red << 16) | (green << 8) | blue;
                image.setRGB(x - 1, y - 1, argb);
            }
        }

        return image;
    }

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

    private static List<FluidEmitter> generateEdgeEmitters(FluidGrid grid, int emitterCount, Random random) {
        List<FluidEmitter> emitters = new ArrayList<>();
        float centerX = (grid.width + 1) / 2.0f;
        float centerY = (grid.height + 1) / 2.0f;

        while (emitters.size() < emitterCount) {
            int side = random.nextInt(4);
            int x;
            int y;

            int offset = EMITTER_RADIUS;

            if (side == 0) { // top
                x = 1 + random.nextInt(grid.width);
                y = 1 + offset;
            } else if (side == 1) { // bottom
                x = 1 + random.nextInt(grid.width);
                y = grid.height - offset;
            } else if (side == 2) { // left
                x = 1 + offset;
                y = 1 + random.nextInt(grid.height);
            } else { // right
                x = grid.width - offset;
                y = 1 + random.nextInt(grid.height);
            }

            float towardCenterX = centerX - x;
            float towardCenterY = centerY - y;

            float candidateAngle = (float) Math.toDegrees(Math.atan2(towardCenterY, towardCenterX));


            float[] emitterColor = NAMESPACE_COLORS[emitters.size() % NAMESPACE_COLORS.length];
            float emissionSpeed = randomRange(random, MIN_EMISSION_SPEED, MAX_EMISSION_SPEED);
            emitters.add(new FluidEmitter(
                    x,
                    y,
                    EMITTER_RADIUS,
                    DEFAULT_DENSITY_RATE,
                    candidateAngle,
                    emissionSpeed,
                    emitterColor[0],
                    emitterColor[1],
                    emitterColor[2]
            ));
        }

        return emitters;
    }

    private static float randomRange(Random random, float min, float max) {
        return min + random.nextFloat() * (max - min);
    }

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
}
