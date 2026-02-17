import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Main {
    private static final float DEFAULT_DENSITY_RATE = 25.0f;
    private static final float MIN_EMISSION_SPEED = 0.2f;
    private static final float MAX_EMISSION_SPEED = 1.2f;

    public static void main(String[] args) {
        SimulationConfig config = parseConfig(args);

        FluidGrid grid = new FluidGrid(config.gridWidth, config.gridHeight, 1.0f / Math.max(config.gridWidth, config.gridHeight));
        SimulationParameters parameters = new SimulationParameters(0.016f, 0.0001f, 0.0001f, 20);

        List<FluidSource> sources = List.of(new FluidSource(grid.width / 2, grid.height / 2, 40.0f));
        List<FluidEmitter> emitters = generateEdgeEmitters(grid, config.emitterCount, new Random());

        FluidSolver solver = new FluidSolver(grid, parameters, sources, emitters);

        System.out.println("Generated " + emitters.size() + " edge emitters:");
        for (int i = 0; i < emitters.size(); i++) {
            FluidEmitter emitter = emitters.get(i);
            System.out.printf(
                    "  #%d at (%d,%d) angle=%.2f rad speed=%.2f%n",
                    i + 1,
                    emitter.gridX,
                    emitter.gridY,
                    emitter.angleRadians,
                    emitter.emissionSpeed
            );
        }

        for (int step = 1; step <= 5; step++) {
            solver.step();
            System.out.println("Step " + step + " divergence RMS=" + solver.computeVelocityDivergenceRms());
        }

        int center = grid.index((grid.width + 1) / 2, (grid.height + 1) / 2);
        System.out.println("Simulation ran. Center density=" + solver.densityField.readValues[center]);

        String outputPath = "density.png";
        saveDensityToPng(grid, solver.densityField, outputPath);
        System.out.println("Saved density image to " + outputPath);
    }

    private static SimulationConfig parseConfig(String[] args) {
        int gridWidth = parsePositiveInt(args, 0, 128, "grid width");
        int gridHeight = parsePositiveInt(args, 1, 128, "grid height");
        int emitterCount = parsePositiveInt(args, 2, 8, "emitter count");

        return new SimulationConfig(gridWidth, gridHeight, emitterCount);
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

    private static void saveDensityToPng(FluidGrid grid, ScalarField densityField, String outputPath) {
        BufferedImage image = new BufferedImage(grid.width, grid.height, BufferedImage.TYPE_INT_ARGB);

        float maxDensity = 0.0f;
        for (int y = 1; y <= grid.height; y++) {
            for (int x = 1; x <= grid.width; x++) {
                int index = grid.index(x, y);
                maxDensity = Math.max(maxDensity, densityField.readValues[index]);
            }
        }

        float normalization = maxDensity > 0.0f ? maxDensity : 1.0f;

        for (int y = 1; y <= grid.height; y++) {
            for (int x = 1; x <= grid.width; x++) {
                int index = grid.index(x, y);
                float normalized = clamp(densityField.readValues[index] / normalization, 0.0f, 1.0f);
                int grayscale = Math.round(normalized * 255.0f);
                int argb = (255 << 24) | (grayscale << 16) | (grayscale << 8) | grayscale;
                image.setRGB(x - 1, y - 1, argb);
            }
        }

        try {
            ImageIO.write(image, "png", new File(outputPath));
        } catch (IOException exception) {
            throw new RuntimeException("Failed to write density PNG to " + outputPath, exception);
        }
    }

    private static List<FluidEmitter> generateEdgeEmitters(FluidGrid grid, int emitterCount, Random random) {
        List<FluidEmitter> emitters = new ArrayList<>();
        float centerX = (grid.width + 1) / 2.0f;
        float centerY = (grid.height + 1) / 2.0f;

        while (emitters.size() < emitterCount) {
            int side = random.nextInt(4);
            int x;
            int y;

            if (side == 0) {
                x = 1 + random.nextInt(grid.width);
                y = 1;
            } else if (side == 1) {
                x = 1 + random.nextInt(grid.width);
                y = grid.height;
            } else if (side == 2) {
                x = 1;
                y = 1 + random.nextInt(grid.height);
            } else {
                x = grid.width;
                y = 1 + random.nextInt(grid.height);
            }

            float towardCenterX = centerX - x;
            float towardCenterY = centerY - y;
            float towardCenterAngle = (float) Math.atan2(towardCenterY, towardCenterX);

            float candidateAngle = towardCenterAngle + randomRange(random, -(float) (Math.PI / 2.0), (float) (Math.PI / 2.0));
            float directionX = (float) Math.cos(candidateAngle);
            float directionY = (float) Math.sin(candidateAngle);
            float inwardDot = directionX * towardCenterX + directionY * towardCenterY;

            if (inwardDot <= 0.0f) {
                continue;
            }

            float emissionSpeed = randomRange(random, MIN_EMISSION_SPEED, MAX_EMISSION_SPEED);
            emitters.add(new FluidEmitter(x, y, DEFAULT_DENSITY_RATE, candidateAngle, emissionSpeed));
        }

        return emitters;
    }

    private static float randomRange(Random random, float min, float max) {
        return min + random.nextFloat() * (max - min);
    }

    private static float clamp(float value, float min, float max) {
        return Math.max(min, Math.min(max, value));
    }

    private record SimulationConfig(int gridWidth, int gridHeight, int emitterCount) {}
}
