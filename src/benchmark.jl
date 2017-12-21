using BenchmarkTools
include("WaveSim.jl")
using WaveSim

function main()
  focus = 0.03  # [m]
  steer = 0.0  # [deg]
  aperture_size = 0.01  # [m]

  trans_delays = WaveSim.delays_from_focus_and_steer(focus, steer, aperture_size)  # elements firing delay

  # run the simulation
  sim_params, images = WaveSim.wavesim(trans_delays)
end

# TODO: FIXME: This doesn't get displayed.
# @benchmark main()
println("Can now execute: \"@benchmark main()\"")
