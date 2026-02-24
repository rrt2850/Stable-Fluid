import java.util.Objects;

/**
 * Represents a wall segment defined by two endpoints and a thickness.
 *
 * <p>The wall spans from {@code startPoint} to {@code endPoint}. Length is
 * derived from those endpoints, while {@code width} controls how thick the
 * wall should be when rasterized into grid cells.</p>
 *
 * @param width wall thickness in grid cells
 * @param startPoint first endpoint of the wall segment
 * @param endPoint second endpoint of the wall segment
 */
public record Wall(int width, WallPoint startPoint, WallPoint endPoint) {

    /**
     * Validates constructor arguments.
     */
    public Wall {
        if (width <= 0) {
            throw new IllegalArgumentException("wall width must be > 0");
        }

        Objects.requireNonNull(startPoint, "startPoint must not be null");
        Objects.requireNonNull(endPoint, "endPoint must not be null");

        if (startPoint.equals(endPoint)) {
            throw new IllegalArgumentException("wall endpoints must be different points");
        }
    }

    /**
     * Computes Euclidean wall length in grid units.
     */
    public float length() {
        int dx = endPoint.x() - startPoint.x();
        int dy = endPoint.y() - startPoint.y();
        return (float) Math.sqrt(dx * dx + dy * dy);
    }

    /**
     * Creates a wall from raw endpoint coordinates.
     */
    public static Wall fromCoordinates(int width, int startX, int startY, int endX, int endY) {
        return new Wall(width, new WallPoint(startX, startY), new WallPoint(endX, endY));
    }
}
