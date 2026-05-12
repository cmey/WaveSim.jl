using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

try
  using CUDA
catch
end

using BenchmarkTools
using WaveSim

function main()
  sim_params = WaveSimParameters(
    focus_depth = 0.03,
    steer_angle = 10.0,
    aperture_size = 0.02,
    temporal_res = 0.5e-6,
    spatial_res = [64, 128],
  )

  trans_delays = WaveSim.delays_from_focus_and_steer(sim_params)
  sim_params = WaveSim.autores(sim_params, trans_delays, multiplier = 0.2, min_spatial_res = (96, 192))

  println("CUDA backend available: ", WaveSim.cuda_backend_available())
  if !WaveSim.cuda_backend_available()
    println("No functional CUDA device detected; benchmark skipped.")
    return nothing
  end

  cpu_images = WaveSim.wavesim(trans_delays, sim_params; backend = :cpu)
  gpu_images = WaveSim.wavesim(trans_delays, sim_params; backend = :cuda)

  max_abs_error = maximum(abs.(cpu_images .- gpu_images))
  rel_error = max_abs_error / max(maximum(abs.(cpu_images)), 1f-6)

  cpu_trial = @benchmark WaveSim.wavesim($trans_delays, $sim_params; backend = :cpu)
  gpu_trial = @benchmark WaveSim.wavesim($trans_delays, $sim_params; backend = :cuda)

  println("Output comparison:")
  println("  max abs error = ", max_abs_error)
  println("  relative error = ", rel_error)

  println("CPU benchmark:")
  display(cpu_trial)

  println("GPU benchmark:")
  display(gpu_trial)

  return nothing
end

main()