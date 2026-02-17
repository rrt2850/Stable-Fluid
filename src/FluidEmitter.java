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
    public void applyVelocity(VectorField velocity, FluidGrid grid) {
        if (!grid.inBounds(gridX, gridY)) {
            throw new IllegalArgumentException("emitter out of bounds: (" + gridX + ", " + gridY + ")");
        }
        int index = grid.index(gridX, gridY);

        float velocityX = (float) Math.cos(angleRadians) * emissionSpeed;
        float velocityY = (float) Math.sin(angleRadians) * emissionSpeed;

        velocity.readVelocityX[index] += velocityX;
        velocity.readVelocityY[index] += velocityY;
    }
}
