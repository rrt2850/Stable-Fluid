/**
 * Emits fluid density and momentum uniformly in all directions from a center point.
 *
 * @param gridX         X-coordinate of the emitter in grid space.
 * @param gridY         Y-coordinate of the emitter in grid space.
 * @param radius        Radius of the emitter influence.
 * @param densityRate   Amount of density injected per second.
 * @param emissionSpeed Base outward speed added to nearby cells.
 * @param red           Red component of injected density in [0, 1].
 * @param green         Green component of injected density in [0, 1].
 * @param blue          Blue component of injected density in [0, 1].
 */
public record RadialFluidEmitter(int gridX, int gridY, int radius, float densityRate, float emissionSpeed,
                                 float red, float green, float blue) {

    private static final int MIN_RADIUS = 1;

    /**
     * Pushes nearby velocity outward from the emitter center in all directions.
     */
    public void applyVelocity(VectorField velocity, FluidGrid grid) {
        if (!grid.inBounds(gridX, gridY)) {
            throw new IllegalArgumentException("radial emitter out of bounds: (" + gridX + ", " + gridY + ")");
        }

        int effectiveRadius = Math.max(radius, MIN_RADIUS);
        float radiusSquared = effectiveRadius * effectiveRadius;

        for (int dy = -effectiveRadius; dy <= effectiveRadius; dy++) {
            for (int dx = -effectiveRadius; dx <= effectiveRadius; dx++) {
                int x = gridX + dx;
                int y = gridY + dy;

                if (!grid.inBounds(x, y)) {
                    continue;
                }

                float distSquared = dx * dx + dy * dy;
                if (distSquared == 0.0f || distSquared > radiusSquared) {
                    continue;
                }

                float distance = (float) Math.sqrt(distSquared);
                float weight = 1.0f - (distance / effectiveRadius);
                float directionX = dx / distance;
                float directionY = dy / distance;

                int index = grid.index(x, y);
                velocity.readVelocityX[index] += directionX * emissionSpeed * weight;
                velocity.readVelocityY[index] += directionY * emissionSpeed * weight;
            }
        }
    }

    /**
     * Adds colored density around the center with a smooth falloff toward the edge.
     */
    public void applyDensity(ScalarField redDensityField, ScalarField greenDensityField, ScalarField blueDensityField,
                             FluidGrid grid, float timeStepSeconds) {
        if (!grid.inBounds(gridX, gridY)) {
            throw new IllegalArgumentException("radial emitter out of bounds: (" + gridX + ", " + gridY + ")");
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
