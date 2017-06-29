# UltraSim.jl, ultrasound simulator.
# Wave propagation computation.
# Christophe MEYER, 2016
module UltraSim

# Configuration
# speed of sound
const c = 1540.0  # [m/s]
# transmit center frequency
const F0 = 3_000_000.0  # [Hz]
# resolution at which to divide the simulation time span
const temporal_res = 0.1 * 1e-6  # [s]
# spatial resolution of the simulation
const spatial_res = [512, 512]  # [pixels]
# field of view, x=0 centered on aperture center, z=0 at aperture plane
const fov = [0.04, 0.04]  # [m]
# time span to simulate
const end_simulation_time = 20.0 * 1e-6  # [s] start at 0 s

# physical length of the transducer array
const transducer_array_size = 0.03  # [m]
# shape of the transmit pulse
pulse_shape_func(phase) = cos(phase)
# length of the transmit pulse
const pulse_length = 2 / F0  # [s] 2 complete cycles


# compute dependent parameters given global configuration
function init(trans_delays)
  wavelength = c/F0 # [m]
  n_transducers = length(trans_delays)
  # spacing between elements
  # const transducer_pitch = 0.1 * 1e-3  # [m]
  transducer_pitch = transducer_array_size / n_transducers
  x_transducers = [transducer_pitch * itrans for itrans in 1:n_transducers]
  tvec = linspace(0.0, end_simulation_time, end_simulation_time / temporal_res)
  image_pitch = fov ./ spatial_res
  return wavelength, x_transducers, tvec, image_pitch
end


# simulate one time step of wave propagation
function simulate_one_time_step!(image, t, image_pitch, x_transducers, trans_delays, wavelength)
  # println("simulate_one_time_step")
  transducers_that_are_firing = find(trans_delays .>= 0)
  trans_coord = [0.0, 0.0]  # [m] space
  pix_coord = [0.0, 0.0]
  for x in 1:spatial_res[1]
    for y in 1:spatial_res[2]
      pix_coord[1] = x * image_pitch[1]  # space
      pix_coord[2] = y * image_pitch[2]  # space

      # for transducers that will fire...
      for i_trans in transducers_that_are_firing
        xt = x_transducers[i_trans]  # index
        trans_coord[1] = xt
        trans_coord[2] = 0.0
        # from pixel to transducer
        dist_to_transducer = sqrt((trans_coord[1] - pix_coord[1])^2 +
                                  (trans_coord[2] - pix_coord[2])^2)
        trans_delay = trans_delays[i_trans]
        time_to_reach = dist_to_transducer / c + trans_delay
        # if the transducer wave reached this pixel...
        if time_to_reach <= t <= time_to_reach+pulse_length
          amp = pulse_shape_func((t - time_to_reach) * F0 * 2 * pi)
        else
          amp = 0.0
        end
        # wave interference
        image[x,y] += amp
      end
    end
  end
end


# run the simulation time steps
function ultrasim(trans_delays)
  wavelength, x_transducers, tvec, image_pitch = init(trans_delays)

  images = zeros(Float32, (length(tvec), spatial_res[1], spatial_res[2]))

  for (i_time, t) in enumerate(tvec)
    image = view(images, i_time, :, :)
    simulate_one_time_step!(image, t, image_pitch, x_transducers, trans_delays, wavelength)
  end

  return images
end

# test simulator
function main()
  # number of elements
  const num_transducers = 5

  trans_delays = zeros(num_transducers)
  images = ultrasim(trans_delays)
  n_sim_time_steps = length(images)
  # while true
  #   fig = figure()
  #   for i_image in 1:n_sim_time_steps
  #     imshow(images[i_image])
  #     sleep(0.1)
  #   end
  #   close(fig)
  # end
  return images
end

end  # module UltraSim
