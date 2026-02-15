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

    public FluidEmitter(int gridX,
                        int gridY,
                        float densityRate,
                        float angleRadians,
                        float emissionSpeed) {
        this.gridX = gridX;
        this.gridY = gridY;
        this.densityRate = densityRate;
        this.angleRadians = angleRadians;
        this.emissionSpeed = emissionSpeed;
    }

    /**
     * Injects density into the fluid.
     *
     * This modifies the authoritative density field.
     */
    public void applyDensity(ScalarField density, FluidGrid grid) {
        int index = grid.index(gridX, gridY);
        density.readValues[index] += densityRate;
    }

    /**
     * Injects momentum into the fluid by adding velocity
     * in the emission direction.
     */
    public void applyVelocity(VectorField velocity, FluidGrid grid) {
        int index = grid.index(gridX, gridY);

        float velocityX = (float) Math.cos(angleRadians) * emissionSpeed;
        float velocityY = (float) Math.sin(angleRadians) * emissionSpeed;

        velocity.readVelocityX[index] += velocityX;
        velocity.readVelocityY[index] += velocityY;
    }
}
