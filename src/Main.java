import javax.imageio.ImageIO;
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
    private static final float DEFAULT_DENSITY_RATE = 40.0f;
    private static final float MIN_EMISSION_SPEED = 10.0f;
    private static final float MAX_EMISSION_SPEED = 15.0f;
    private static final float TIMESTEP = 0.040f;
    private static final float VISCOSITY = 0.00005f;
    private static final float DIFFUSION_RATE = 0.0001f;
    private static final int SOLVER_ITERATIONS = 40;

    private static final int DEFAULT_SIMULATION_STEPS = 100;
    private static final int DEFAULT_EMITTER_COUNT = 12;
    private static final int MP4_FRAMES_PER_SECOND = 30;


    // 12 namespace colors represented as RGB triplets in the [0, 1] range.
    private static final float[][] NAMESPACE_COLORS = {
            {0.0f, 0.69f, 0.94f},  // azure
            {0.0f, 0.78f, 0.58f},  // mint
            {0.55f, 0.80f, 0.26f}, // lime
            {0.98f, 0.82f, 0.20f}, // yellow
            {1.0f, 0.63f, 0.0f},   // orange
            {0.95f, 0.36f, 0.20f}, // vermillion
            {0.87f, 0.22f, 0.52f}, // magenta
            {0.67f, 0.27f, 0.87f}, // purple
            {0.31f, 0.42f, 0.94f}, // indigo
            {0.18f, 0.68f, 0.98f}, // sky
            {0.0f, 0.73f, 0.75f},  // teal
            {0.42f, 0.75f, 0.45f}  // green
    };

    public static void main(String[] args) {
        long programStartTime = System.nanoTime();
        SimulationConfig config = parseConfig(args);

        FluidGrid grid = new FluidGrid(config.gridWidth, config.gridHeight, 1.0f / Math.max(config.gridWidth, config.gridHeight));
        SimulationParameters parameters = new SimulationParameters(TIMESTEP, VISCOSITY, DIFFUSION_RATE, SOLVER_ITERATIONS);

        List<FluidSource> sources = List.of(new FluidSource(grid.width / 2, grid.height / 2, 40.0f));
        List<FluidEmitter> emitters = generateEdgeEmitters(grid, config.emitterCount, new Random());

        FluidSolver solver = new FluidSolver(grid, parameters, sources, emitters);

        System.out.println("Generated " + emitters.size() + " edge emitters:");
        for (int i = 0; i < emitters.size(); i++) {
            FluidEmitter emitter = emitters.get(i);
            int red = Math.round(clamp(emitter.red, 0.0f, 1.0f) * 255.0f);
            int green = Math.round(clamp(emitter.green, 0.0f, 1.0f) * 255.0f);
            int blue = Math.round(clamp(emitter.blue, 0.0f, 1.0f) * 255.0f);
            System.out.printf(
                    "  #%d at (%d,%d) angle=%.2f rad speed=%.2f color=(%d,%d,%d)%n",
                    i + 1,
                    emitter.gridX,
                    emitter.gridY,
                    emitter.angleRadians,
                    emitter.emissionSpeed,
                    red,
                    green,
                    blue
            );
        }

        List<BufferedImage> frames = config.exportVideo ? new ArrayList<>() : List.of();
        for (int step = 1; step <= config.simulationSteps; step++) {
            solver.step();
            if (config.exportVideo) {
                frames.add(createDensityImage(grid, solver));
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
            String videoOutputPath = "density-diffusion.mp4";
            boolean mp4Exported = saveDensityToMp4(frames, videoOutputPath);
            if (mp4Exported) {
                System.out.println("Saved density animation to " + videoOutputPath);
            } else {
                System.out.println("Skipped MP4 export because ffmpeg is not available on PATH.");
            }
        } else {
            System.out.println("Skipped MP4 export. Pass a 5th argument of 'true' to enable video output.");
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

        try {
            ImageIO.write(image, "png", new File(outputPath));
        } catch (IOException exception) {
            throw new RuntimeException("Failed to write density PNG to " + outputPath, exception);
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

    private static boolean saveDensityToMp4(List<BufferedImage> frames, String outputPath) {
        if (frames.isEmpty()) {
            throw new IllegalArgumentException("At least one frame is required to create an MP4.");
        }

        Path tempFramesDirectory = null;

        try {
            tempFramesDirectory = Files.createTempDirectory("stable-fluid-frames-");

            for (int i = 0; i < frames.size(); i++) {
                Path framePath = tempFramesDirectory.resolve(String.format("frame-%05d.png", i));
                ImageIO.write(frames.get(i), "png", framePath.toFile());
            }

            Path tempOutput = tempFramesDirectory.resolve("density-diffusion.mp4");
            ProcessBuilder ffmpegBuilder = new ProcessBuilder(
                    "ffmpeg",
                    "-y",
                    "-framerate", String.valueOf(MP4_FRAMES_PER_SECOND),
                    "-i", tempFramesDirectory.resolve("frame-%05d.png").toString(),
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
        } finally {
            if (tempFramesDirectory != null) {
                deleteDirectoryRecursively(tempFramesDirectory.toFile());
            }
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

            int offset = 3; // same as radius

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
            float towardCenterAngle = (float) Math.atan2(towardCenterY, towardCenterX);

            float maxDeviation = (float) Math.toRadians(30.0);
            float candidateAngle = towardCenterAngle
                    + randomRange(random, -maxDeviation, maxDeviation);
            float directionX = (float) Math.cos(candidateAngle);
            float directionY = (float) Math.sin(candidateAngle);
            float inwardDot = directionX * towardCenterX + directionY * towardCenterY;

            if (inwardDot <= 0.0f) {
                continue;
            }

            float[] emitterColor = NAMESPACE_COLORS[emitters.size() % NAMESPACE_COLORS.length];
            float emissionSpeed = randomRange(random, MIN_EMISSION_SPEED, MAX_EMISSION_SPEED);
            emitters.add(new FluidEmitter(
                    x,
                    y,
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
