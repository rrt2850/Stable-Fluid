/**
 * Pulls nearby fluid toward its center and removes density, acting like a sink
 *
 * @param gridX             X-coordinate of the vortex center in grid space
 * @param gridY             Y-coordinate of the vortex center in grid space
 * @param radius            Radius of influence
 * @param suctionStrength   Inward pull magnitude applied to velocity per second
 * @param absorptionRate    Density removed per second near the center
 * @param swirlStrength     Tangential spin strength for whirlpool-like circular flow per second
 * @param swirlDirection    Spin orientation, use 1 for clockwise and -1 for counterclockwise
 */
public record Vortex(int gridX, int gridY, int radius, float suctionStrength, float absorptionRate,
                     float swirlStrength, int swirlDirection) {

    private static final int MIN_RADIUS = 1;

    /**
     * Applies vortex influence in one pass (velocity + all 3 density channels)
     * to avoid repeated wall-occlusion ray checks.
     */
    public void applyInfluence(
            VectorField velocity,
            ScalarField redDensity,
            ScalarField greenDensity,
            ScalarField blueDensity,
            FluidGrid grid,
            float timeStepSeconds,
            boolean[] wallMask
    ) {
        if (!grid.inBounds(gridX, gridY)) {
            throw new IllegalArgumentException("vortex out of bounds: (" + gridX + ", " + gridY + ")");
        }
        if (swirlDirection != 1 && swirlDirection != -1) {
            throw new IllegalArgumentException("swirlDirection must be 1 (clockwise) or -1 (counterclockwise)");
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

                if (isOccludedByWall(x, y, grid, wallMask)) {
                    continue;
                }

                float distance = (float) Math.sqrt(distSquared);
                float suctionWeight = suctionWeight(distance, effectiveRadius);
                int index = grid.index(x, y);

                float removedDensity = absorptionRate * timeStepSeconds * suctionWeight;
                redDensity.readValues[index] = Math.max(0.0f, redDensity.readValues[index] - removedDensity);
                greenDensity.readValues[index] = Math.max(0.0f, greenDensity.readValues[index] - removedDensity);
                blueDensity.readValues[index] = Math.max(0.0f, blueDensity.readValues[index] - removedDensity);

                if (distSquared == 0.0f) {
                    continue;
                }

                float swirlWeight = swirlWeight(distance, effectiveRadius);
                float directionToCenterX = -dx / distance;
                float directionToCenterY = -dy / distance;

                // Tangential vector for clockwise rotation: perpendicular to radial direction.
                float tangentX = dy / distance;
                float tangentY = -dx / distance;

                // Scale by timestep so strengths are expressed in "per-second" terms.
                float pullStrength = suctionStrength * suctionWeight * timeStepSeconds;
                float spinStrength = swirlStrength * swirlWeight * timeStepSeconds;

                velocity.readVelocityX[index] += directionToCenterX * pullStrength + tangentX * spinStrength * swirlDirection;
                velocity.readVelocityY[index] += directionToCenterY * pullStrength + tangentY * spinStrength * swirlDirection;
            }
        }
    }

    /**
     * True when a wall lies between the vortex center and the target cell.
     *
     * <p>Uses a Bresenham walk across grid cells so vortex forces cannot pass
     * through wall tiles but can still bend around wall gaps via advection.</p>
     */
    private boolean isOccludedByWall(int targetX, int targetY, FluidGrid grid, boolean[] wallMask) {
        if (wallMask == null || wallMask.length == 0) {
            return false;
        }

        int x = gridX;
        int y = gridY;
        int dx = Math.abs(targetX - x);
        int dy = Math.abs(targetY - y);
        int stepX = Integer.compare(targetX, x);
        int stepY = Integer.compare(targetY, y);

        int error = dx - dy;

        while (x != targetX || y != targetY) {
            int doubledError = error * 2;

            if (doubledError > -dy) {
                error -= dy;
                x += stepX;
            }
            if (doubledError < dx) {
                error += dx;
                y += stepY;
            }

            if (!grid.inBounds(x, y)) {
                return true;
            }

            if (wallMask[grid.index(x, y)]) {
                return true;
            }
        }

        return false;
    }

    /**
     * Keeps visible swirl across the vortex while still tapering at the boundary.
     */
    private static float swirlWeight(float distance, int effectiveRadius) {
        float normalizedDistance = distance / effectiveRadius;

        // Keep tangential motion strongest around the middle of the vortex and
        // weaker at the boundary so particles are less likely to orbit outside.
        float edgeToCenter = Math.max(0.0f, 1.0f - normalizedDistance);
        return normalizedDistance * edgeToCenter;
    }

    /**
     * Creates a radial falloff where suction starts weak near the edge and grows toward the center.
     */
    private static float suctionWeight(float distance, int effectiveRadius) {
        float normalizedDistance = distance / effectiveRadius;
        float centerProximity = Math.max(0.0f, 1.0f - normalizedDistance);

        // Keep a non-zero pull near the edge so trajectories spiral inward
        // instead of settling into mostly circular orbits.
        return 0.2f + 0.8f * centerProximity;
    }

}
