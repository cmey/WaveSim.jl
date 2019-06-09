using WaveSim
using StaticArrays

# Define simulation parameters (use many default values, see WaveSimParameters).
const fov = @SVector [0.8 * 15e-2, 15e-2]  # [m]
sim_params = WaveSimParameters(
    tx_frequency = 2.75e6,
    pulse_cycles = 5,
    steer_angle = 0.0,  # [deg]
    fov = fov,
    end_simulation_time = 2 * 56.0e-6,  # [s] starts at 0 s
    aperture_size = 0.019968,
    # focus_depth = 0.11,  # [m]
    focus_depth = Inf,  # [m]
    temporal_res = 0.05e-6 ,  # [s]
    spatial_res = [256, 256],  # [pixels]
    dbrange = 30,
    orientation = :vertical
);

# Compute focusing delays for the elements of the phased array.
trans_delays = WaveSim.delays_from_focus_and_steer(sim_params);

# Run the simulation.
images = WaveSim.wavesim(trans_delays, sim_params);
beam_energy_map, transmit_time_map = WaveSim.beam_energy_map_and_transmit_time_map(images, sim_params);

# Display results.
include("src/view.jl")
imshowall(images, beam_energy_map, transmit_time_map, sim_params);

# Save results.
include("src/save.jl")
saveall(images, beam_energy_map, transmit_time_map, sim_params);
