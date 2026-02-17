/**
 * Represents a fluid emitter that injects mass and momentum
 * into the simulation at a fixed position.
 *
 * Conceptually: "fluid shoots out at this angle and speed".
 */
public class FluidEmitter {

    /**
     * X-coordinate of the emitter in grid space.
     */
    public final int gridX;

    /**
     * Y-coordinate of the emitter in grid space.
     */
    public final int gridY;

    /**
     * Amount of density injected per time step.
     */
    public final float densityRate;

    /**
     * Direction of emitted fluid in radians.
     */
    public final float angleRadians;

    /**
     * Speed at which fluid is emitted.
     */
    public final float emissionSpeed;

    /**
     * Red component of injected density in the [0, 1] range.
     */
    public final float red;

    /**
     * Green component of injected density in the [0, 1] range.
     */
    public final float green;

    /**
     * Blue component of injected density in the [0, 1] range.
     */
    public final float blue;

    public FluidEmitter(int gridX,
                        int gridY,
                        float densityRate,
                        float angleRadians,
                        float emissionSpeed,
                        float red,
                        float green,
                        float blue) {
        this.gridX = gridX;
        this.gridY = gridY;
        this.densityRate = densityRate;
        this.angleRadians = angleRadians;
        this.emissionSpeed = emissionSpeed;
        this.red = red;
        this.green = green;
        this.blue = blue;
    }

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
        int radius = 3;
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
