import java.util.List;

public record EmitterHolder(List<FluidEmitter> fluidEmitters,
                            List<RadialFluidEmitter> radialFluidEmitters,
                            List<Vortex> vortexes,
                            List<Wall> walls) {
}


