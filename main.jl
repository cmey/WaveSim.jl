include("ultrasim.jl")
include("view.jl")
using UltraSim

# Run the simulator and display results.
function main()
  focus = 0.03  # [m]
  aperture_size = 0.02  # [m]

  # Compute focusing delays for the elements of the phased array.
  trans_delays = UltraSim.delays_from_focus(focus, aperture_size)

  # Run the simulation.
  images = UltraSim.ultrasim(trans_delays)
  beam_energy_map, transmit_time_map = UltraSim.beam_energy_map_and_transmit_time_map(images)

  # Display results.
  imshowall(images, beam_energy_map, transmit_time_map)

  return images
end


images = main()
