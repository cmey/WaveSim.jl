include("UltraSim.jl")
include("view.jl")
using UltraSim

# test simulator
function main()
  # number of elements
  const num_transducers = 50
  # elements firing delay
  trans_delays = zeros(num_transducers)
  # run the simulation
  images = UltraSim.ultrasim(trans_delays)

  return images
end


images = main()

imshowall(images)
