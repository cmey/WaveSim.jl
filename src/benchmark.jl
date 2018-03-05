using BenchmarkTools
include("WaveSim.jl")
using WaveSim

function main()
  sim_params = WaveSimParameters(
    focus_depth = 0.03,  # [m]
    steer_angle = 10.0,  # [deg]
    aperture_size = 0.01,  # [m]
    temporal_res = 0.5e-6,  # [s]
    spatial_res = [32, 64]  # [pixels]
  )

  # Compute focusing delays for the elements of the phased array.
  trans_delays = WaveSim.delays_from_focus_and_steer(sim_params)

  # Run the simulation.
  images = WaveSim.wavesim(trans_delays, sim_params)

  return images, sim_params
end

display(@benchmark main())
