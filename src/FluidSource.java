/**
 * Represents a localized source of density or momentum.
 *
 * Used to inject material or energy into the fluid.
 */
public class FluidSource {

    /**
     * X-coordinate of the source in grid space
     */
    public final int gridX;

    /**
     * Y-coordinate of the source in grid space
     */
    public final int gridY;

    /**
     * Strength of the source (amount added per step)
     */
    public final float strength;

    /**
     * Creates a new fluid source
     *
     * @param gridX     x-position on the grid
     * @param gridY     y-position on the grid
     * @param strength  amount of density or velocity added
     */
    public FluidSource(int gridX, int gridY, float strength) {
        this.gridX = gridX;
        this.gridY = gridY;
        this.strength = strength;
    }

    /**
     * Applies this source to a density field.
     *
     * @param density density field to modify
     * @param grid    simulation grid
     */
    public void applyToDensity(ScalarField density, FluidGrid grid) {
        int index = grid.index(gridX, gridY);
        density.readValues[index] += strength;
    }
}
