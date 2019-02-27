include("view.jl")
using WaveSim

# Run the simulator and display results.
function main()
  # Define simulation parameters (use many default values, see WaveSimParameters).
  sim_params = WaveSimParameters(
    focus_depth = 0.03,  # [m]
    steer_angle = 10.0,  # [deg]
    aperture_size = 0.02,  # [m]
    temporal_res = 0.1e-6,  # [s]
    spatial_res = [256, 512]  # [pixels]
  )

  # Compute focusing delays for the elements of the phased array.
  trans_delays = WaveSim.delays_from_focus_and_steer(sim_params)

  # Run the simulation.
  images = WaveSim.wavesim(trans_delays, sim_params)
  beam_energy_map, transmit_time_map = WaveSim.beam_energy_map_and_transmit_time_map(images, sim_params)

  # Display results.
  imshowall(images, beam_energy_map, transmit_time_map, sim_params)

  return beam_energy_map, transmit_time_map, images
end

beam_energy_map, transmit_time_map, images = main()
