using BenchmarkTools
include("UltraSim.jl")
using UltraSim

function main()
  # number of elements
  const num_transducers = 5
  # elements firing delay
  trans_delays = zeros(num_transducers)
  # run the simulation
  images = UltraSim.ultrasim(trans_delays)
end

@benchmark images = main()
