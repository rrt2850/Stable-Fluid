/**
 * Represents a fluid emitter that injects mass and momentum
 * into the simulation at a fixed position.
 * <p>
 * Conceptually: "fluid shoots out at this angle and speed".
 *
 * @param gridX         X-coordinate of the emitter in grid space.
 * @param gridY         Y-coordinate of the emitter in grid space.
 * @param radius        The radius of the emitter
 * @param densityRate   Amount of density injected per time step.
 * @param angleDegrees  Direction of emitted fluid in degrees.
 * @param emissionSpeed Speed at which fluid is emitted.
 * @param red           Red component of injected density in the [0, 1] range.
 * @param green         Green component of injected density in the [0, 1] range.
 * @param blue          Blue component of injected density in the [0, 1] range.
 */
public record FluidEmitter(int gridX, int gridY, int radius, float densityRate, float angleDegrees, float emissionSpeed,
                           float red, float green, float blue) {

    /** Safety floor so math never divides by zero radius. */
    private static final int MIN_RADIUS = 1;
    /** Controls how much wider emitters reduce per-cell push strength. */
    private static final float VELOCITY_RADIUS_FALLOFF_EXPONENT = 0.5f;

    /**
     * Injects momentum into the fluid over a circular region
     * with radial falloff so the emission behaves like a jet instead
     * of a single-cell impulse. Velocity is damped for larger radii
     * so widening the emitter does not linearly boost total momentum.
     */
    public void applyVelocity(VectorField velocity, FluidGrid grid) {
        if (!grid.inBounds(gridX, gridY)) {
            throw new IllegalArgumentException(
                    "emitter out of bounds: (" + gridX + ", " + gridY + ")");
        }

        int effectiveRadius = Math.max(radius, MIN_RADIUS);
        float radiusSquared = effectiveRadius * effectiveRadius;

        // Base emission direction vector
        float velocityRadiusScale = 1.0f / (float) Math.pow(effectiveRadius, VELOCITY_RADIUS_FALLOFF_EXPONENT);
        double angleRadians = Math.toRadians(angleDegrees);
        float velocityPerWeightX = (float) Math.cos(angleRadians) * emissionSpeed * velocityRadiusScale;
        float velocityPerWeightY = (float) Math.sin(angleRadians) * emissionSpeed * velocityRadiusScale;

        for (int dy = -effectiveRadius; dy <= effectiveRadius; dy++) {
            for (int dx = -effectiveRadius; dx <= effectiveRadius; dx++) {

                int x = gridX + dx;
                int y = gridY + dy;

                if (!grid.inBounds(x, y)) {
                    continue;
                }

                float distSquared = dx * dx + dy * dy;

                if (distSquared > radiusSquared) {
                    continue; // outside circular region
                }

                // Linear falloff (1 at center → 0 at edge)
                float distance = (float) Math.sqrt(distSquared);
                float weight = 1.0f - (distance / effectiveRadius);

                int index = grid.index(x, y);

                velocity.readVelocityX[index] += velocityPerWeightX * weight;
                velocity.readVelocityY[index] += velocityPerWeightY * weight;
            }
        }
    }

    /**
     * Injects density over the same circular region as velocity.
     * Density is not normalized by area, so increasing emitter radius
     * makes the stream visibly wider and denser.
     */
    public void applyDensity(ScalarField redDensityField, ScalarField greenDensityField, ScalarField blueDensityField,
                             FluidGrid grid, float timeStepSeconds) {
        if (!grid.inBounds(gridX, gridY)) {
            throw new IllegalArgumentException(
                    "emitter out of bounds: (" + gridX + ", " + gridY + ")");
        }

        int effectiveRadius = Math.max(radius, MIN_RADIUS);
        float radiusSquared = effectiveRadius * effectiveRadius;

        float densityPerWeight = timeStepSeconds * densityRate;

        for (int dy = -effectiveRadius; dy <= effectiveRadius; dy++) {
            for (int dx = -effectiveRadius; dx <= effectiveRadius; dx++) {
                int x = gridX + dx;
                int y = gridY + dy;

                if (!grid.inBounds(x, y)) {
                    continue;
                }

                float distSquared = dx * dx + dy * dy;
                if (distSquared > radiusSquared) {
                    continue;
                }

                float distance = (float) Math.sqrt(distSquared);
                float weight = 1.0f - (distance / effectiveRadius);
                float weightedDensity = densityPerWeight * weight;

                int index = grid.index(x, y);
                redDensityField.readValues[index] += weightedDensity * red;
                greenDensityField.readValues[index] += weightedDensity * green;
                blueDensityField.readValues[index] += weightedDensity * blue;
            }
        }
    }
}
