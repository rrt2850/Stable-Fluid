import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class EmitterMaker {
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
    public static List<FluidEmitter> generateEdgeEmitters(FluidGrid grid, int emitterCount, Random random) {
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
                float jitter = Utils.randomRange(
                        random,
                        -Constants.EMITTER_ANGLE_VARIATION_DEGREES,
                        Constants.EMITTER_ANGLE_VARIATION_DEGREES
                );
                candidateAngle = centerAngle + jitter;
                attempts++;
            } while (
                    (!isValidForWall(candidateAngle, side)) &&
                            attempts < 25
            );

            float[] color = Constants.NAMESPACE_COLORS[
                    emitters.size() % Constants.NAMESPACE_COLORS.length
                    ];

            float speed = Utils.randomRange(
                    random,
                    Constants.MIN_EMISSION_SPEED,
                    Constants.MAX_EMISSION_SPEED
            );

            emitters.add(new FluidEmitter(
                    x,
                    y,
                    radius,
                    Constants.DEFAULT_DENSITY_RATE,
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
     * Computes an emitter radius based on grid dimensions.
     *
     * <p>The radius scales with grid size but is clamped to a practical
     * minimum and maximum.</p>
     *
     * @param grid simulation grid
     * @return emitter radius in grid cells
     */
    private static int computeEmitterRadius(FluidGrid grid) {
        int radiusFromRatio = Math.round(Math.min(grid.width, grid.height) * Constants.EMITTER_RADIUS_RATIO);
        return Math.min(Constants.MAX_EMITTER_RADIUS, Math.max(Constants.MIN_EMITTER_RADIUS, radiusFromRatio));
    }

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
        float a = Utils.normalizeAngle(angleDeg);

        // Reject near-tangential flow along the wall
        float tangent = switch (side) {
            case 0, 1 -> 0f;    // top/bottom → horizontal tangent
            case 2, 3 -> 90f;   // left/right → vertical tangent
            default -> throw new IllegalArgumentException();
        };

        if (Utils.near(a, tangent, Constants.WALL_TANGENT_EXCLUSION) ||
                Utils.near(a, tangent + 180f, Constants.WALL_TANGENT_EXCLUSION)) {
            return false;
        }

        // Reject corner-sniping diagonals
        return switch (side) {
            case 0 -> !Utils.near(a, 45f, Constants.CORNER_EXCLUSION) &&
                    !Utils.near(a, 135f, Constants.CORNER_EXCLUSION);
            case 1 -> !Utils.near(a, 225f, Constants.CORNER_EXCLUSION) &&
                    !Utils.near(a, 315f, Constants.CORNER_EXCLUSION);
            case 2 -> !Utils.near(a, 45f, Constants.CORNER_EXCLUSION) &&
                    !Utils.near(a, 315f, Constants.CORNER_EXCLUSION);
            case 3 -> !Utils.near(a, 135f, Constants.CORNER_EXCLUSION) &&
                    !Utils.near(a, 225f, Constants.CORNER_EXCLUSION);
            default -> true;
        };
    }
}
