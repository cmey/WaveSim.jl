include("ultrasim.jl")
include("view.jl")
using UltraSim

# test simulator
function main()
  focus = 0.03
  aperture_size = 0.01
  trans_delays = UltraSim.delays_from_focus(focus, aperture_size)  # elements firing delay

  # run the simulation
  images = UltraSim.ultrasim(trans_delays)

  # display result
  imshowall(images)

  return images
end


images = main()
