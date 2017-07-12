include("ultrasim.jl")
include("view.jl")
using UltraSim

# Run the simulator and display results.
function main()
  focus = 0.03  # [m]
  steer = 10.0  # [deg]
  aperture_size = 0.02  # [m]
  dbrange = 40

  # Compute focusing delays for the elements of the phased array.
  trans_delays = UltraSim.delays_from_focus_and_steer(focus, steer, aperture_size)

  # Run the simulation.
  sim_params, images = UltraSim.ultrasim(trans_delays)
  beam_energy_map, transmit_time_map = UltraSim.beam_energy_map_and_transmit_time_map(images)

  # Display results.
  imshowall(sim_params, images, beam_energy_map, transmit_time_map, dbrange)

  return images
end


images = main()
