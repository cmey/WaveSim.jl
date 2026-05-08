using WaveSim, StaticArrays, InteractiveUtils

sim_params = WaveSimParameters(
  focus_depth = 0.03,
  steer_angle = 10.0,
  aperture_size = 0.01,
  temporal_res = 0.5e-6,
  spatial_res = @SVector [64, 128]
)

trans_delays = WaveSim.delays_from_focus_and_steer(sim_params)
transducers, tvec, pulse_length, apodization_vec, pix_xs, pix_ys, firing = WaveSim.init(trans_delays, sim_params)

image = zeros(Float32, (sim_params.spatial_res[1], sim_params.spatial_res[2]))
t = tvec[10]

println("Running @code_warntype for simulate_one_time_step! — output follows:")
InteractiveUtils.@code_warntype WaveSim.simulate_one_time_step!(image, t, trans_delays, pulse_length, sim_params.tx_frequency, sim_params.c, sim_params.spatial_res, sim_params.pulse_shape_func, apodization_vec, sim_params.directivity_func, sim_params.transducer_pitch, pix_xs, pix_ys, transducers, firing)
