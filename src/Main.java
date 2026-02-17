import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Main {
    private static final float DEFAULT_DENSITY_RATE = 25.0f;
    private static final float MIN_EMISSION_SPEED = 0.2f;
    private static final float MAX_EMISSION_SPEED = 1.2f;

    public static void main(String[] args) {
        int emitterCount = parseEmitterCount(args);

        FluidGrid grid = new FluidGrid(128, 128, 1.0f / 128.0f);
        SimulationParameters parameters = new SimulationParameters(0.016f, 0.0001f, 0.0001f, 20);

        List<FluidSource> sources = List.of(new FluidSource(grid.width / 2, grid.height / 2, 40.0f));
        List<FluidEmitter> emitters = generateEdgeEmitters(grid, emitterCount, new Random());

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
    }

    private static int parseEmitterCount(String[] args) {
        if (args.length == 0) {
            return 8;
        }

        try {
            int parsed = Integer.parseInt(args[0]);
            if (parsed <= 0) {
                throw new IllegalArgumentException("Emitter count must be greater than zero.");
            }
            return parsed;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException("First argument must be an integer emitter count.", exception);
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
}
