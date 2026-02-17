import java.util.List;

public class Main {
    public static void main(String[] args) {
        FluidGrid grid = new FluidGrid(128, 128, 1.0f / 128.0f);
        SimulationParameters parameters = new SimulationParameters(0.016f, 0.0001f, 0.0001f, 20);

        FluidSource source = new FluidSource(64, 64, 40.0f);
        FluidEmitter emitter = new FluidEmitter(64, 64, 10.0f, 0f, 0.5f);

        FluidSolver solver = new FluidSolver(
                grid,
                parameters,
                List.of(source),
                List.of(emitter)
        );

        for (int step = 1; step <= 5; step++) {
            solver.step();
        }

        int center = grid.index(64, 64);
        System.out.println("Simulation ran. Center density=" + solver.densityField.readValues[center]);
    }
}
