/**
 * Holds all tunable parameters controlling the simulation
 */
public class SimulationParameters {

    /**
     * Time step used for integration
     * Smaller values improve stability at the cost of performance
     */
    public float timeStep;

    /**
     * Viscosity of the fluid, how goopy the fluid is
     * Controls velocity diffusion.
     */
    public float viscosity;

    /**
     * Diffusion rate for scalar quantities like density
     */
    public float diffusionRate;

    /**
     * Number of Gauss–Seidel iterations used in linear solvers
     * Higher values improve accuracy but increase cost
     */
    public int linearSolverIterations;

    /**
     * Constructs a parameter set with explicit values, more detail in
     * SimulationParameters.java
     */
    public SimulationParameters(float timeStep,
                                float viscosity,
                                float diffusionRate,
                                int linearSolverIterations) {
        this.timeStep = timeStep;
        this.viscosity = viscosity;
        this.diffusionRate = diffusionRate;
        this.linearSolverIterations = linearSolverIterations;
    }

    /**
     * Constructs a parameter set with some defaults
     */
    public SimulationParameters() {
        this(0.016f, 0.0001f, 0.0001f, 20);
    }
}
