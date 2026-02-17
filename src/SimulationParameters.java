/**
 * Holds all tunable parameters controlling the simulation
 */
public class SimulationParameters {

    /**
     * Time step used for integration
     * Smaller values improve stability at the cost of performance
     */
    private float timeStep;

    /**
     * Viscosity of the fluid, how goopy the fluid is
     * Controls velocity diffusion.
     */
    private float viscosity;

    /**
     * Diffusion rate for scalar quantities like density
     */
    private float diffusionRate;

    /**
     * Number of Gauss–Seidel iterations used in linear solvers
     * Higher values improve accuracy but increase cost
     */
    private int linearSolverIterations;

    /**
     * Strength of vorticity confinement force used to re-inject small-scale swirl.
     * 0 disables the force entirely.
     */
    private float vorticityConfinement;

    /**
     * Constructs a parameter set with explicit values, more detail in
     * SimulationParameters.java
     */
    public SimulationParameters(float timeStep,
                                float viscosity,
                                float diffusionRate,
                                int linearSolverIterations,
                                float vorticityConfinement) {
        setTimeStep(timeStep);
        setViscosity(viscosity);
        setDiffusionRate(diffusionRate);
        setLinearSolverIterations(linearSolverIterations);
        setVorticityConfinement(vorticityConfinement);
    }

    /**
     * Constructs a parameter set with some defaults
     */
    public SimulationParameters() {
        this(0.016f, 0.0001f, 0.0001f, 20, 0.0f);
    }

    public float getTimeStep() {
        return timeStep;
    }

    public void setTimeStep(float timeStep) {
        if (timeStep <= 0f) {
            throw new IllegalArgumentException("timeStep must be > 0");
        }
        this.timeStep = timeStep;
    }

    public float getViscosity() {
        return viscosity;
    }

    public void setViscosity(float viscosity) {
        if (viscosity < 0f) {
            throw new IllegalArgumentException("viscosity must be >= 0");
        }
        this.viscosity = viscosity;
    }

    public float getDiffusionRate() {
        return diffusionRate;
    }

    public void setDiffusionRate(float diffusionRate) {
        if (diffusionRate < 0f) {
            throw new IllegalArgumentException("diffusionRate must be >= 0");
        }
        this.diffusionRate = diffusionRate;
    }

    public int getLinearSolverIterations() {
        return linearSolverIterations;
    }

    public void setLinearSolverIterations(int linearSolverIterations) {
        if (linearSolverIterations <= 0) {
            throw new IllegalArgumentException("linearSolverIterations must be > 0");
        }
        this.linearSolverIterations = linearSolverIterations;
    }

    public float getVorticityConfinement() {
        return vorticityConfinement;
    }

    public void setVorticityConfinement(float vorticityConfinement) {
        if (vorticityConfinement < 0f) {
            throw new IllegalArgumentException("vorticityConfinement must be >= 0");
        }
        this.vorticityConfinement = vorticityConfinement;
    }
}
