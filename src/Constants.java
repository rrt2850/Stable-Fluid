public class Constants {
    // Fluid properties
    public static final float DEFAULT_DENSITY_RATE = 1.2f;
    public static final float MIN_EMISSION_SPEED = 0.7f;
    public static final float MAX_EMISSION_SPEED = 1.1f;
    public static final float TIMESTEP = 0.030f;
    public static final float VISCOSITY = 0.0000001f;
    public static final float DIFFUSION_RATE = 0.00006f;
    public static final int SOLVER_ITERATIONS = 25;
    public static final float VORTICITY_CONFINEMENT = 10.0f;

    public static final int DEFAULT_SIMULATION_STEPS = 100;
    public static final int DEFAULT_EMITTER_COUNT = 8;

    // Defaults for PowerShell runner
    public static final int DEFAULT_GRID_WIDTH = 800;
    public static final int DEFAULT_GRID_HEIGHT = 400;

    // Final still export resolution (upscaled from sim)
    public static final int FINAL_STILL_WIDTH = 2400;
    public static final int FINAL_STILL_HEIGHT = 1200;

    public static final int MP4_FRAMES_PER_SECOND = 30;
    public static final int INTERMITTENT_SNAPSHOT_INTERVAL = 50;

    public static final long RANDOM_SEED = System.currentTimeMillis();

    // Emitter settings
    public static final float EMITTER_RADIUS_RATIO = 0.012f;
    public static final int MIN_EMITTER_RADIUS = 8;
    public static final int MAX_EMITTER_RADIUS = 60;
    public static final float EMITTER_ANGLE_VARIATION_DEGREES = 60.0f;
    public static final float WALL_TANGENT_EXCLUSION = 10f;
    public static final float CORNER_EXCLUSION = 5f;

    // Good seeds :)
    // 1771389195668L
    // 1771400966868L
    // 1771402741333L
    // 1771453995487L

    // 12 namespace colors represented as RGB triplets in the [0, 1] range.
    public static final float[][] NAMESPACE_COLORS = {
            {0.0f, 0.69f, 0.94f},
            {0.0f, 0.78f, 0.58f},
            {0.55f, 0.80f, 0.26f},
            {0.98f, 0.82f, 0.20f},
            {1.0f, 0.63f, 0.0f},
            {0.95f, 0.36f, 0.20f},
            {0.87f, 0.22f, 0.52f},
            {0.67f, 0.27f, 0.87f},
            {0.31f, 0.42f, 0.94f},
            {0.18f, 0.68f, 0.98f},
            {0.0f, 0.73f, 0.75f},
            {0.42f, 0.75f, 0.45f}
    };
}
