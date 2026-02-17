/**
 * Represents the discrete simulation domain for the fluid.
 *
 * The grid is a 2D Cartesian lattice with additional "ghost cells"
 * around the edges to support boundary conditions.
 *
 * All scalar and vector fields in the simulation index into this grid.
 */
public class FluidGrid {
    public final int width; // Number of active cells horizontally, not including ghost cells
    public final int height; // Number of simulation cells vertically, not including ghost cells
    public final int totalCellCount; // The total number of cells including ghost cells
    public final float cellSize; // The physical size of each grid cell, used when converting between phyiscal units and grid quantities

    /**
     * Constructs a grid
     *
     * @param width     number of simulation cells along the x-axis
     * @param height    number of simulation cells along the y-axis
     * @param cellSize  physical size of each grid cell
     */
    public FluidGrid(int width, int height, float cellSize) {
        if (width <= 0 || height <= 0) {
            throw new IllegalArgumentException("width and height must be > 0");
        }
        if (cellSize <= 0f) {
            throw new IllegalArgumentException("cellSize must be > 0");
        }
        this.width = width;
        this.height = height;
        this.cellSize = cellSize;
        this.totalCellCount = (width + 2) * (height + 2); // ghost cells
    }

    /**
     * Constructs a grid using default values
     *
     * @param cellSize physical size of each grid cell
     */
    public FluidGrid(float cellSize) {
        this(1920, 1080, cellSize);
    }

    /**
     * Converts 2D grid coordinates into a 1D array index
     *
     * The conversion includes ghost cells, so valid indices start at (1,1)
     * and end at (width, height)
     *
     * @param x grid x-coordinate
     * @param y grid y-coordinate
     * @return index into a 1D field array
     */
    public int index(int x, int y) {
        if (x < 0 || x > width + 1 || y < 0 || y > height + 1) {
            throw new IllegalArgumentException(
                    "grid coordinates out of range: (" + x + ", " + y + ")"
            );
        }
        return x + (width + 2) * y;
    }

    /**
     * Checks whether a grid coordinate lies inside the active simulation domain.
     *
     * Ghost cells are excluded.
     *
     * @param x grid x-coordinate
     * @param y grid y-coordinate
     * @return true if the cell is part of the simulation domain
     */
    public boolean inBounds(int x, int y) {
        return x >= 1 && y >= 1 && x <= width && y <= height;
    }
}
