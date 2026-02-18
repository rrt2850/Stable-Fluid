/**
 * Pulls nearby fluid toward its center and removes density, acting like a sink.
 *
 * @param gridX             X-coordinate of the vortex center in grid space.
 * @param gridY             Y-coordinate of the vortex center in grid space.
 * @param radius            Radius of influence.
 * @param suctionStrength   Inward pull magnitude applied to velocity.
 * @param absorptionRate    Density removed per second near the center.
 */
public record Vortex(int gridX, int gridY, int radius, float suctionStrength, float absorptionRate) {

    private static final int MIN_RADIUS = 1;

    public void applyVelocity(VectorField velocity, FluidGrid grid, float timeStepSeconds) {
        if (!grid.inBounds(gridX, gridY)) {
            throw new IllegalArgumentException("vortex out of bounds: (" + gridX + ", " + gridY + ")");
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

                float directionToCenterX = -dx / distance;
                float directionToCenterY = -dy / distance;
                float pullStrength = suctionStrength * weight * timeStepSeconds;

                int index = grid.index(x, y);
                velocity.readVelocityX[index] += directionToCenterX * pullStrength;
                velocity.readVelocityY[index] += directionToCenterY * pullStrength;
            }
        }
    }

    public void absorbDensity(ScalarField densityField, FluidGrid grid, float timeStepSeconds) {
        if (!grid.inBounds(gridX, gridY)) {
            throw new IllegalArgumentException("vortex out of bounds: (" + gridX + ", " + gridY + ")");
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
                if (distSquared > radiusSquared) {
                    continue;
                }

                float distance = (float) Math.sqrt(distSquared);
                float weight = 1.0f - (distance / effectiveRadius);
                float removedDensity = absorptionRate * timeStepSeconds * weight;

                int index = grid.index(x, y);
                densityField.readValues[index] = Math.max(0.0f, densityField.readValues[index] - removedDensity);
            }
        }
    }
}
