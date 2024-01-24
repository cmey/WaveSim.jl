using WaveSim

# Run the simulator and display results.
function main()
  # Define simulation parameters (use many default values, see WaveSimParameters).
  sim_params = WaveSimParameters(
    focus_depth = 0.03,  # [m]
    steer_angle = 10.0,  # [deg]
    aperture_size = 0.02,  # [m]
    # tx_frequency = 3_000_000.0,  # [Hz]
    # focus_depth = Inf,  # [m]
    # steer_angle = 0.0,  # [deg]
    # transducer_pitch = 205e-6,  # [m]
    # transducer_array_size = 0.01312,  # [m]
    # aperture_size = 0.01312,  # [m]
  )

  # Compute focusing delays for the elements of the phased array.
  trans_delays = WaveSim.delays_from_focus_and_steer(sim_params)

  # Compute optimized "best" spatial and temporal parameters.
  sim_params = WaveSim.autores(sim_params, trans_delays)

  # Run the simulation.
  images = WaveSim.wavesim(trans_delays, sim_params)
  beam_energy_map, transmit_time_map = WaveSim.beam_energy_map_and_transmit_time_map(images, sim_params)

  return beam_energy_map, transmit_time_map, images, sim_params
end

beam_energy_map, transmit_time_map, images, sim_params = main()

# Display results.
include("view.jl")
imshowall(images, beam_energy_map, transmit_time_map, sim_params)

# Save results.
include("save.jl")
saveall(images, beam_energy_map, transmit_time_map, sim_params, "images")
