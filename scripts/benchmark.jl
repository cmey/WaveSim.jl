using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using BenchmarkTools
using Profile
using WaveSim

function run_serial_profile_workload(trans_delays, sim_params)
  transducers, transducer_normals, tvec, pulse_length, apodization_vec, pix_coord_xs, pix_coord_ys, transducers_that_are_firing = WaveSim.init(trans_delays, sim_params)
  images = zeros(Float32, (sim_params.spatial_res[1], sim_params.spatial_res[2], length(tvec)))
  effective_transducer_pitch = WaveSim.slice_directivity_pitch(sim_params.beamplot_axes, sim_params.transducer_pitch, sim_params.transducer_pitch_elevation)

  for i_time in 1:length(tvec)
    t = tvec[i_time]
    image = view(images, :, :, i_time)
    WaveSim.simulate_one_time_step!(image, t, trans_delays, pulse_length, sim_params.tx_frequency, sim_params.c, sim_params.spatial_res, sim_params.pulse_shape_func, apodization_vec, sim_params.directivity_func, effective_transducer_pitch, sim_params.attenuation_coefficient, pix_coord_xs, pix_coord_ys, transducers, transducer_normals, sim_params.beamplot_axes, transducers_that_are_firing)
  end

  return images
end

function main()
  sim_params = WaveSimParameters(
    focus_depth = 0.03,
    steer_angle = 10.0,
    aperture_size = 0.02,
    temporal_res = 0.5e-6,
    spatial_res = [64, 128]
  )

  trans_delays = WaveSim.delays_from_focus_and_steer(sim_params)
  sim_params = WaveSim.autores(sim_params, trans_delays, multiplier=0.2, min_spatial_res=(96, 192))

  println("Benchmark input:")
  println("  spatial_res = ", sim_params.spatial_res)
  println("  transducers = ", length(trans_delays))

  println("Benchmarking steady-state runtime...")
  trial = @benchmark WaveSim.wavesim($trans_delays, $sim_params)
  display(trial)

  println("Warmup run...")
  images = WaveSim.wavesim(trans_delays, sim_params)

  println("Profiling run (serial, no progress output)... this may take a few seconds")
  Profile.clear()
  @profile begin
    for i in 1:2
      images = run_serial_profile_workload(trans_delays, sim_params)
    end
  end

  println("Profile results:")
  println("\nTree view (deeper stack):")
  Profile.print(format=:tree, C=false, groupby=:none, maxdepth=30, noisefloor=1.0)

  println("\nHotspots (flat view):")
  Profile.print(format=:flat, C=false, groupby=:none, sortedby=:count, mincount=1)

  return images, sim_params
end

main()