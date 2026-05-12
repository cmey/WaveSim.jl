module WaveSimCUDAExt

using CUDA
using Parameters: @unpack

import WaveSim

const CUDA_THREADS = (8, 8, 4)

function cuda_backend_available_impl()
  return CUDA.functional()
end

@inline function _cuda_pixel_coordinates(beamplot_axes::Symbol, lateral::Float32, depth::Float32)
  return WaveSim.slice_pixel_coordinates(beamplot_axes, lateral, depth)
end

function _wavesim_cuda_kernel!(
  images,
  temporal_res,
  trans_delays,
  pulse_length,
  tx_frequency,
  c,
  pulse_shape_func,
  apodization_vec,
  directivity_func,
  transducer_pitch,
  attenuation_coefficient,
  pix_coord_xs,
  pix_coord_ys,
  transducers,
  transducer_normals,
  beamplot_axes,
  transducers_that_are_firing,
)
  x = (blockIdx().x - 1) * blockDim().x + threadIdx().x
  y = (blockIdx().y - 1) * blockDim().y + threadIdx().y
  t = (blockIdx().z - 1) * blockDim().z + threadIdx().z

  nx = size(images, 1)
  ny = size(images, 2)
  nt = size(images, 3)
  if x > nx || y > ny || t > nt
    return
  end

  one_over_c = 1.0f0 / c
  α_factor = - (attenuation_coefficient * tx_frequency / 1f6) / 20.0f0

  pix_coord_x = pix_coord_xs[x]
  pix_coord_y = pix_coord_ys[y]
  pix_x, pix_y, pix_z = _cuda_pixel_coordinates(beamplot_axes, pix_coord_x, pix_coord_y)
  current_time = (t - 1) * temporal_res

  amp = 0.0f0
  @inbounds for idx in 1:length(transducers_that_are_firing)
    i_trans = transducers_that_are_firing[idx]
    t_coord = transducers[i_trans]
    n_coord = transducer_normals[i_trans]
    trans_delay = trans_delays[i_trans]

    dx = pix_x - t_coord[1]
    dy = pix_y - t_coord[2]
    dz = pix_z - t_coord[3]
    dist_to_transducer_not_zero = max(sqrt(dx * dx + dy * dy + dz * dz), 1.0f-6)
    wave_spreading_factor = 1.0f0 / dist_to_transducer_not_zero
    time_to_reach = dist_to_transducer_not_zero * one_over_c + trans_delay

    if time_to_reach <= current_time <= time_to_reach + pulse_length
      cos_theta = clamp((dx * n_coord[1] + dy * n_coord[2] + dz * n_coord[3]) / dist_to_transducer_not_zero, -1.0f0, 1.0f0)
      directivity_factor = 0.0f0
      if cos_theta >= 0.0f0
        θ = acos(cos_theta)
        directivity_factor = directivity_func(θ, tx_frequency, c, transducer_pitch)
      end

      dist_cm = dist_to_transducer_not_zero * 100.0f0
      attenuation_factor = 10.0f0 ^ (α_factor * dist_cm)
      amp += apodization_vec[i_trans] * pulse_shape_func((current_time - time_to_reach) * tx_frequency * 2.0f0) * wave_spreading_factor * directivity_factor * attenuation_factor
    end
  end

  @inbounds images[x, y, t] = amp
  return
end

function wavesim_cuda_impl(trans_delays, sim_params)
  if !CUDA.functional()
    error("CUDA.jl is available but no functional CUDA device was detected.")
  end

  @unpack tx_frequency, c, temporal_res, spatial_res, pulse_shape_func, transducer_pitch, transducer_pitch_elevation, attenuation_coefficient, beamplot_axes, directivity_func = sim_params
  transducers, transducer_normals, tvec, pulse_length, apodization_vec, pix_coord_xs, pix_coord_ys, transducers_that_are_firing = WaveSim.init(trans_delays, sim_params)
  effective_transducer_pitch = WaveSim.slice_directivity_pitch(beamplot_axes, transducer_pitch, transducer_pitch_elevation)

  CUDA.allowscalar(false)

  trans_delays_d = CuArray(trans_delays)
  apodization_vec_d = CuArray(apodization_vec)
  pix_coord_xs_d = CuArray(pix_coord_xs)
  pix_coord_ys_d = CuArray(pix_coord_ys)
  transducers_d = CuArray(transducers)
  transducer_normals_d = CuArray(transducer_normals)
  transducers_that_are_firing_d = CuArray(transducers_that_are_firing)

  images_d = CUDA.zeros(Float32, spatial_res[1], spatial_res[2], length(tvec))
  blocks = (
    cld(spatial_res[1], CUDA_THREADS[1]),
    cld(spatial_res[2], CUDA_THREADS[2]),
    cld(length(tvec), CUDA_THREADS[3]),
  )

  CUDA.@sync @cuda threads=CUDA_THREADS blocks=blocks _wavesim_cuda_kernel!(
    images_d,
    temporal_res,
    trans_delays_d,
    pulse_length,
    tx_frequency,
    c,
    pulse_shape_func,
    apodization_vec_d,
    directivity_func,
    effective_transducer_pitch,
    attenuation_coefficient,
    pix_coord_xs_d,
    pix_coord_ys_d,
    transducers_d,
    transducer_normals_d,
    beamplot_axes,
    transducers_that_are_firing_d,
  )

  return Array(images_d)
end

end