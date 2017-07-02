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

  # Display results.
  imshowall(images)

  return images
end


images = main()
