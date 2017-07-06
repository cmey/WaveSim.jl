using BenchmarkTools
include("ultrasim.jl")
using UltraSim

function main()
  focus = 0.03
  aperture_size = 0.01
  trans_delays = UltraSim.delays_from_focus(focus, aperture_size)  # elements firing delay
  # run the simulation
  sim_params, images = UltraSim.ultrasim(trans_delays)
end

# TODO: FIXME: This doesn't get displayed.
# @benchmark images = main()
