/**
 * Represents a fluid emitter that shoots out a fluid at a given location
 * and direction.
 * <p>
 *
 *
 * @param gridX         X-coordinate of the emitter in grid space
 * @param gridY         Y-coordinate of the emitter in grid space
 * @param radius        The radius of the emitter
 * @param densityRate   Amount of density injected per time step
 * @param angleDegrees  Direction of emitted fluid in degrees
 * @param emissionSpeed Speed fluid is emitted
 * @param red           Red component of injected density in the [0, 1] range
 * @param green         Green component of injected density in the [0, 1] range
 * @param blue          Blue component of injected density in the [0, 1] range
 */
public record FluidEmitter(int gridX, int gridY, int radius, float densityRate, float angleDegrees, float emissionSpeed,
                           float red, float green, float blue) {

    /** Safety floor so math never divides by zero radius. */
    private static final int MIN_RADIUS = 1;

    /** Controls how much wider emitters reduce per-cell push strength. */
    private static final float VELOCITY_RADIUS_FALLOFF_EXPONENT = 0.5f;

    /**
     * Injects momentum into the fluid over a circular region
     * with radial falloff. Velocity is damped for larger radiuses
     * so widening the emitter does not linearly boost total momentum
     */
    public void applyVelocity(VectorField velocity, FluidGrid grid) {

        // Make sure the grid is in bounds and the radius isn't too small

        if (!grid.inBounds(gridX, gridY)) {
            throw new IllegalArgumentException(
                    "emitter out of bounds: (" + gridX + ", " + gridY + ")");
        }

        int radius = Math.max(this.radius, MIN_RADIUS);
        float radiusSquared = radius * radius;


        // Scale velocity contribution by radius so that larger emitters
        // do not inject disproportionately more momentum
        // The exponent allows tuning how strongly radius affects strength.
        float velocityRadiusScale =
                1.0f / (float) Math.pow(radius, VELOCITY_RADIUS_FALLOFF_EXPONENT);

        double angleRadians = Math.toRadians(angleDegrees);

        // Unit direction vector (cos, sin) scaled by emission speed
        // and radius normalization. This represents the maximum
        // velocity contribution at the center of the emitter.
        float velocityPerWeightX =
                (float) Math.cos(angleRadians) * emissionSpeed * velocityRadiusScale;
        float velocityPerWeightY =
                (float) Math.sin(angleRadians) * emissionSpeed * velocityRadiusScale;


        // Iterate over a square region [-radius, +radius] in grid space,
        // then later discard cells outside the circular radius.
        for (int dy = -radius; dy <= radius; dy++) {
            for (int dx = -radius; dx <= radius; dx++) {

                // Convert local offsets to absolute grid coordinates.
                int x = gridX + dx;
                int y = gridY + dy;

                // Skip cells that fall outside the simulation domain.
                if (!grid.inBounds(x, y)) {
                    continue;
                }

                // Squared distance from emitter center to this cell
                float distSquared = dx * dx + dy * dy;

                // Reject cells outside the circular influence region
                // trimming the square iteration into a disk
                if (distSquared > radiusSquared) {
                    continue;
                }

                // Actual Euclidean distance from center
                float distance = (float) Math.sqrt(distSquared);

                // Linear radial falloff:
                //   - weight = 1.0 at the center
                //   - weight = 0.0 exactly at the radius
                //
                // This produces a smooth velocity gradient
                float weight = 1.0f - (distance / radius);

                // Convert 2D grid coordinates to the underlying
                // 1D array index used by the velocity field
                int index = grid.index(x, y);

                // Accumulate velocity into the grid cell
                // Add so multiple emitters or sim steps can contribute
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
    public void applyDensity(ScalarField redDensityField,
                             ScalarField greenDensityField,
                             ScalarField blueDensityField,
                             FluidGrid grid,
                             float timeStepSeconds) {

        // Make sure the grid is in bounds and the radius isn't too small
        if (!grid.inBounds(gridX, gridY)) {
            throw new IllegalArgumentException(
                    "emitter out of bounds: (" + gridX + ", " + gridY + ")");
        }

        int radius = Math.max(this.radius, MIN_RADIUS);
        float radiusSquared = radius * radius;

        // Base density contribution per unit weight, scaled by timestep
        // This represents the maximum density injected at the center
        float densityPerWeight = timeStepSeconds * densityRate;

        // Iterate over a square region [-radius, +radius] in grid space,
        // then later discard cells outside the circular radius.
        for (int dy = -radius; dy <= radius; dy++) {
            for (int dx = -radius; dx <= radius; dx++) {

                // Convert local offsets to absolute grid coordinates.
                int x = gridX + dx;
                int y = gridY + dy;

                // Skip cells that fall outside the simulation domain.
                if (!grid.inBounds(x, y)) {
                    continue;
                }

                // Squared distance from emitter center to this cell
                float distSquared = dx * dx + dy * dy;

                // Reject cells outside the circular influence region
                // trimming the square iteration into a disk
                if (distSquared > radiusSquared) {
                    continue;
                }

                // Actual Euclidean distance from center
                float distance = (float) Math.sqrt(distSquared);

                // Linear radial falloff:
                //   - weight = 1.0 at the center
                //   - weight = 0.0 exactly at the radius
                //
                // This produces a smooth density gradient
                float weight = 1.0f - (distance / radius);

                float weightedDensity = densityPerWeight * weight;

                // Convert 2D grid coordinates to the underlying
                // 1D array index used by the scalar fields
                int index = grid.index(x, y);

                // Accumulate density into the grid cell
                // Add so multiple emitters or sim steps can contribute
                redDensityField.readValues[index]   += weightedDensity * red;
                greenDensityField.readValues[index] += weightedDensity * green;
                blueDensityField.readValues[index]  += weightedDensity * blue;
            }
        }
    }
}
