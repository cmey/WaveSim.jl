include("wavesim.jl")
include("view.jl")
using WaveSim

# Run the simulator and display results.
function main()
  focus = 0.03  # [m]
  steer = 10.0  # [deg]
  aperture_size = 0.02  # [m]
  dbrange = 40

  # Compute focusing delays for the elements of the phased array.
  trans_delays = WaveSim.delays_from_focus_and_steer(focus, steer, aperture_size)

  # Run the simulation.
  sim_params, images = WaveSim.wavesim(trans_delays)
  beam_energy_map, transmit_time_map = WaveSim.beam_energy_map_and_transmit_time_map(images)

  # Display results.
  imshowall(sim_params, images, beam_energy_map, transmit_time_map, dbrange)

  return images
end


images = main()
