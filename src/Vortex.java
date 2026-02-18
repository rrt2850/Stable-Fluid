/**
 * Pulls nearby fluid toward its center and removes density, acting like a sink.
 *
 * @param gridX             X-coordinate of the vortex center in grid space.
 * @param gridY             Y-coordinate of the vortex center in grid space.
 * @param radius            Radius of influence.
 * @param suctionStrength   Inward pull magnitude applied to velocity.
 * @param absorptionRate    Density removed per second near the center.
 * @param swirlStrength     Tangential spin strength for whirlpool-like circular flow.
 */
public record Vortex(int gridX, int gridY, int radius, float suctionStrength, float absorptionRate,
                     float swirlStrength) {

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

                // Tangential vector for clockwise rotation: perpendicular to radial direction.
                float tangentX = dy / distance;
                float tangentY = -dx / distance;

                float pullStrength = suctionStrength * weight * timeStepSeconds;
                float spinStrength = swirlStrength * weight * timeStepSeconds;

                int index = grid.index(x, y);
                velocity.readVelocityX[index] += directionToCenterX * pullStrength + tangentX * spinStrength;
                velocity.readVelocityY[index] += directionToCenterY * pullStrength + tangentY * spinStrength;
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
