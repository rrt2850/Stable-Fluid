/**
 * Applies boundary conditions to simulation fields
 *
 * Boundary conditions ensure correct physical behavior
 * at the edges of the simulation domain
 */
public class BoundaryHandler {

    /**
     * Identifies what a field physically represents so edges can be treated correctly.
     */
    public enum BoundaryType{
        SCALAR,
        H_VELOCITY,
        V_VELOCITY
    }

    /**
     * Applies physical boundary conditions to a simulation field by populating
     * the grid's ghost cells (the extra cells around the edge of the grid used for math)
     *
     * <p>This method fills the ghost cells with values so that difference operations
     *     don't have null values when they try to do n+1 at the edge of the grid</p>
     *
     * <p>Boundary behavior depends on the type of field:</p>
     * <ul>
     *   <li><b>SCALAR</b> – Scalar quantities (e.g. density, pressure) are mirrored
     *       at the boundary, enforcing a zero normal gradient (no flux through walls)</li>
     *   <li><b>H_VELOCITY</b> – Horizontal (x-axis) velocity is negated at left and
     *       right walls, enforcing no penetration through vertical boundaries.
     *       Values are mirrored at top and bottom walls</li>
     *   <li><b>V_VELOCITY</b> – Vertical (y-axis) velocity is negated at top and
     *       bottom walls, enforcing no penetration through horizontal boundaries.
     *       Values are mirrored at left and right walls</li>
     * </ul>
     *
     * <p>Corner ghost cells are set to the average of their neighboring edge cells
     * to maintain numerical stability and avoid corner artifacts</p>
     *
     * <p>This method must be called every time a solver step writes to a field
     * that will later be accessed using neighbor stencils</p>
     *
     * @param boundaryType identifies the physical meaning of the field and
     *                     determines how boundary constraints are applied
     * @param field        the field array to constrain, including ghost cells
     * @param grid         simulation grid describing layout and bounds
     */
    public void applyBoundaries(BoundaryType boundaryType,
                                float[] field,
                                FluidGrid grid) {

        final int w = grid.width;
        final int h = grid.height;

        // Top and bottom boundaries
        for (int x = 1; x <= w; x++) {
            int bottom      = grid.index(x, 0);
            int bottomInner = grid.index(x, 1);

            int top      = grid.index(x, h + 1);
            int topInner = grid.index(x, h);

            if (boundaryType == BoundaryType.V_VELOCITY) {
                field[bottom] = -field[bottomInner];
                field[top]    = -field[topInner];
            } else {
                field[bottom] = field[bottomInner];
                field[top]    = field[topInner];
            }
        }

        // Left and right boundaries
        for (int y = 1; y <= h; y++) {
            int left      = grid.index(0, y);
            int leftInner = grid.index(1, y);

            int right      = grid.index(w + 1, y);
            int rightInner = grid.index(w, y);

            if (boundaryType == BoundaryType.H_VELOCITY) {
                field[left]  = -field[leftInner];
                field[right] = -field[rightInner];
            } else {
                field[left]  = field[leftInner];
                field[right] = field[rightInner];
            }
        }

        // Corner cells: average adjacent edge cells
        field[grid.index(0, 0)] =
                0.5f * (field[grid.index(1, 0)] +
                        field[grid.index(0, 1)]);

        field[grid.index(0, h + 1)] =
                0.5f * (field[grid.index(1, h + 1)] +
                        field[grid.index(0, h)]);

        field[grid.index(w + 1, 0)] =
                0.5f * (field[grid.index(w, 0)] +
                        field[grid.index(w + 1, 1)]);

        field[grid.index(w + 1, h + 1)] =
                0.5f * (field[grid.index(w, h + 1)] +
                        field[grid.index(w + 1, h)]);
    }

}
