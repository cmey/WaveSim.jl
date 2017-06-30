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
const temporal_res = 0.1e-6  # [s]
# spatial resolution of the simulation
const spatial_res = [512, 512]  # [pixels]
# field of view, x=0 centered on aperture center, z=0 at aperture plane
const fov = [4e-2, 4e-2]  # [m]
# time span to simulate
const end_simulation_time = 20.0e-6  # [s] starts at 0 s

# physical length of the transducer array
const transducer_array_size = 0.03  # [m]
# spacing between physical elements of transducer array
const transducer_pitch = 208e-6  # [m]
# shape of the transmit pulse
pulse_shape_func(phase) = cos(phase)
# length of the transmit pulse
const pulse_length = 2 / F0  # [s] 2 complete cycles


# compute dependent parameters given global configuration
function init(trans_delays)
  wavelength = c/F0 # [m]
  n_transducers = length(trans_delays)
  x_transducers = [transducer_pitch * itrans for itrans in 1:n_transducers]  # [m]
  x_transducers += fov[1] / 2 - mean(extrema(x_transducers))  # center around FOV in x
  time_vec = collect(0.0:temporal_res:end_simulation_time)
  image_pitch = fov ./ spatial_res
  return wavelength, x_transducers, time_vec, image_pitch
end


# simulate one time step of wave propagation
function simulate_one_time_step!(image, t, image_pitch, x_transducers, trans_delays, wavelength)
  # println("simulate_one_time_step")
  transducers_that_are_firing = find(trans_delays .>= 0)
  trans_coord = [0.0, 0.0]  # [m] space
  pix_coord = [0.0, 0.0]
  for y in 1:spatial_res[2]
    pix_coord[2] = y * image_pitch[2]  # [m]
    for x in 1:spatial_res[1]
      pix_coord[1] = x * image_pitch[1]  # [m]

      # for transducers that will fire...
      for i_trans in transducers_that_are_firing
        xt = x_transducers[i_trans]  # index
        trans_delay = trans_delays[i_trans]
        trans_coord[1] = xt
        trans_coord[2] = 0.0
        # from pixel to transducer
        dist_to_transducer = sqrt((trans_coord[1] - pix_coord[1])^2 +
                                  (trans_coord[2] - pix_coord[2])^2)
        time_to_transducer = dist_to_transducer / c
        time_to_reach = time_to_transducer + trans_delay
        # if the transducer wave reached this pixel...
        if time_to_reach <= t <= time_to_reach+pulse_length
          amp = pulse_shape_func((t - time_to_reach) * F0 * 2pi)
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

  images = zeros(Float32, (spatial_res[1], spatial_res[2], length(tvec)))

  for (i_time, t) in enumerate(tvec)
    image = view(images, :, :, i_time)
    simulate_one_time_step!(image, t, image_pitch, x_transducers, trans_delays, wavelength)
  end

  return images
end

# test simulator
function main()
  # number of elements
  const num_transducers = 50

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
