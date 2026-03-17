/**
 * Pulls nearby fluid toward its center and removes density, acting like a sink
 *
 * @param gridX             X-coordinate of the vortex center in grid space
 * @param gridY             Y-coordinate of the vortex center in grid space
 * @param radius            Radius of influence
 * @param suctionStrength   Inward pull magnitude applied to velocity per second
 * @param absorptionRate    Density removed per second near the center
 * @param swirlStrength     Tangential spin strength for whirlpool-like circular flow per second
 */
public record Vortex(int gridX, int gridY, int radius, float suctionStrength, float absorptionRate,
                     float swirlStrength) {

    private static final int MIN_RADIUS = 1;

    /**
     * Pulls velocity inward while adding sideways spin to create a whirlpool effect
     */
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
                float suctionWeight = suctionWeight(distance, effectiveRadius);
                float swirlWeight = swirlWeight(distance, effectiveRadius);

                float directionToCenterX = -dx / distance;
                float directionToCenterY = -dy / distance;

                // Tangential vector for clockwise rotation: perpendicular to radial direction.
                float tangentX = dy / distance;
                float tangentY = -dx / distance;

                // Scale by timestep so strengths are expressed in "per-second" terms.
                float pullStrength = suctionStrength * suctionWeight * timeStepSeconds;
                float spinStrength = swirlStrength * swirlWeight * timeStepSeconds;

                int index = grid.index(x, y);
                velocity.readVelocityX[index] += directionToCenterX * pullStrength + tangentX * spinStrength;
                velocity.readVelocityY[index] += directionToCenterY * pullStrength + tangentY * spinStrength;
            }
        }
    }

    /**
     * Removes density near the vortex center so material appears to be swallowed
     */
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
                float weight = suctionWeight(distance, effectiveRadius);
                float removedDensity = absorptionRate * timeStepSeconds * weight;

                int index = grid.index(x, y);
                densityField.readValues[index] = Math.max(0.0f, densityField.readValues[index] - removedDensity);
            }
        }
    }


    /**
     * Keeps visible swirl across the vortex while still tapering at the boundary.
     */
    private static float swirlWeight(float distance, int effectiveRadius) {
        float normalizedDistance = distance / effectiveRadius;
        return Math.max(0.0f, 1.0f - normalizedDistance);
    }

    /**
     * Creates a radial falloff where suction starts weak near the edge and grows toward the center.
     */
    private static float suctionWeight(float distance, int effectiveRadius) {
        float normalizedDistance = distance / effectiveRadius;
        float centerProximity = Math.max(0.0f, 1.0f - normalizedDistance);

        // Squared falloff creates a smoother gradient near the outer radius while remaining strongest at center.
        return centerProximity * centerProximity;
    }
}
