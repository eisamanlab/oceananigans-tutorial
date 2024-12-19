using Pkg
Pkg.activate(".")

using Dates
using Random
using Printf

using Oceananigans
using Oceananigans.Units: minute, minutes, hour, hours, day, days, meter, meters, kilometer, kilometers

# --------------------------------------------------
# computing architecture
# --------------------------------------------------
architecture = CPU()

# --------------------------------------------------
# grid
# --------------------------------------------------
const Nx = 100  # number of x points 
const Ny = 100  # number of y points 
const Nz = 1    # number of z points 

const Lx = 100kilometers # (m) domain x extent
const Ly = 100kilometers # (m) domain y extent
const Lz = 10meters      # (m) domain depth

grid = RectilinearGrid(
    architecture, 
    topology = (Periodic, Periodic, Bounded),
    size = (Nx, Nx, Nz),
    halo = (3,3,2),
    x = (0, Lx),
    y = (0, Ly),
    z = (-Lz, 0)
)

# --------------------------------------------------
# coriolis
# --------------------------------------------------
coriolis = FPlane(f=1e-4)

# --------------------------------------------------
# closure
# --------------------------------------------------
closure = AnisotropicMinimumDissipation()

# --------------------------------------------------
# seawater buoyancy
# --------------------------------------------------
buoyancy = SeawaterBuoyancy(
    equation_of_state = LinearEquationOfState(
        thermal_expansion = 2e-4,
        haline_contraction = 8e-4)
)

# --------------------------------------------------
# surface stress
# --------------------------------------------------
u₁₀ = 10    # m s⁻¹, average wind velocity 10 meters above the ocean
ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
cᴰ = 2.5e-3 # dimensionless drag coefficient
ρₐ = 1.225  # kg m⁻³, average density of air at sea-level
  
@inline wind_stress_x(x, y, t, p) = (- ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀)) * rand()
@inline wind_stress_y(x, y, t, p) = (- ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀)) * rand()

# --------------------------------------------------
# boundary conditions
# --------------------------------------------------
u_top_bcs = FluxBoundaryCondition(wind_stress_x, parameters=(k=4π, ω=3.0, τ=1e-4))
v_top_bcs = FluxBoundaryCondition(wind_stress_y, parameters=(k=4π, ω=3.0, τ=1e-4))

u_bcs = FieldBoundaryConditions(top = u_top_bcs) 
v_bcs = FieldBoundaryConditions(top = v_top_bcs)

T_bcs = FieldBoundaryConditions()
S_bcs = FieldBoundaryConditions()
c_bcs = FieldBoundaryConditions()

# --------------------------------------------------
# instantiate model
# --------------------------------------------------
model = NonhydrostaticModel(; grid, buoyancy,
                            timestepper = :QuasiAdamsBashforth2,
                            advection = WENO(),
                            tracers = (:T, :S, :c),
                            coriolis = coriolis,
                            closure = closure,
                            boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs, c=c_bcs))


# --------------------------------------------------
# initial conditions
# --------------------------------------------------
const u_initial = 0.0
const v_initial = 0.0
const T_initial = 20.0
const S_initial = 35.0
const c_initial = 10.0

#=
# this code is to set concentration at specific grid cell
# Get the center indices (x, y) 
i_center = Int(grid.Nx / 2)  # Center index in the x-direction
j_center = Int(grid.Ny / 2)  # Center index in the y-direction

# Create the passive tracer and initialize it with some value at central grid cell
c_initial = model.tracers.c .+ rand()

#center_x, center_y = Int(grid.Nx / 2), Int(grid.Ny / 2)
#c_initial[center_x, center_y, grid.Nz] = 1000.0
=#

set!(model, u=u_initial, v=v_initial, T=T_initial, S=S_initial, c=c_initial)

# --------------------------------------------------
# simulation and outputs
# --------------------------------------------------
simulation = Simulation(model, Δt=10minutes, stop_time=30days)

wizard = TimeStepWizard(cfl=0.2)

simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))


progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|u|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim), prettytime(sim), prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.u), prettytime(sim.run_wall_time))

add_callback!(simulation, progress_message, IterationInterval(20))


# --------------------------------------------------
# output writer
# --------------------------------------------------
filename = "output_$(string(today()))"

simulation.output_writers[:surface_slice_writer] = NetCDFOutputWriter(
    model, 
    (; model.velocities.u, model.velocities.v, model.tracers.S, model.tracers.T, model.tracers.c),
    filename = filename * ".nc",
    schedule = AveragedTimeInterval(1hour, window=1hour),
    indices=(:, :, grid.Nz),
    overwrite_existing = true
)

# --------------------------------------------------
# run simulation
# --------------------------------------------------
run!(simulation)
