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
 * @param angleRadians  Direction of emitted fluid in radians.
 * @param emissionSpeed Speed at which fluid is emitted.
 * @param red           Red component of injected density in the [0, 1] range.
 * @param green         Green component of injected density in the [0, 1] range.
 * @param blue          Blue component of injected density in the [0, 1] range.
 */
public record FluidEmitter(int gridX, int gridY, int radius, float densityRate, float angleRadians, float emissionSpeed,
                           float red, float green, float blue) {

    /**
     * Injects momentum into the fluid by adding velocity
     * in the emission direction.
     */
    /**
     * Injects momentum into the fluid over a small circular region
     * with radial falloff so the emission behaves like a jet instead
     * of a single-cell impulse.
     */
    public void applyVelocity(VectorField velocity, FluidGrid grid) {
        if (!grid.inBounds(gridX, gridY)) {
            throw new IllegalArgumentException(
                    "emitter out of bounds: (" + gridX + ", " + gridY + ")");
        }

        // Radius of injection region (tweak 2–5 depending on grid resolution)
        float radiusSquared = radius * radius;

        // Base emission direction vector
        float baseVX = (float) Math.cos(angleRadians) * emissionSpeed;
        float baseVY = (float) Math.sin(angleRadians) * emissionSpeed;

        for (int dy = -radius; dy <= radius; dy++) {
            for (int dx = -radius; dx <= radius; dx++) {

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
                float weight = 1.0f - (distance / radius);

                int index = grid.index(x, y);

                velocity.readVelocityX[index] += baseVX * weight;
                velocity.readVelocityY[index] += baseVY * weight;
            }
        }
    }
}
