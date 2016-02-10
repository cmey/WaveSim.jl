using PyPlot
#using SIUnits
#using SIUnits.ShortUnits

# Config
const c = 1540# * m / s
const F0 = 3_000_000# * Hz
const end_simulation_time = 0.00001# * ms # start at 0 s

const image_res = [128, 128]# * m^-1
const image_fov = [0.01, 0.01]# * m
const n_time_steps = 30

const n_transducers = 5
pulse_shape_func(phase) = cos(phase)
const pulse_length = 1 / F0 # 1 complete cycle
const pitch_transducers = 0.0002# * μm

function simulate_one_time_step!(image, t, image_pitch, x_transducers, wavelength)

  for x in range(1, image_res[1])
    for y in range(1, image_res[2])
      pix_index = [x, y]
      pix_coord = pix_index .* image_pitch

      for xt in x_transducers
        trans_coord = [xt, 0]# * m
        dist_to_transducer = norm(trans_coord .- pix_coord)
        time_to_reach = dist_to_transducer / c
        if time_to_reach <= t <= time_to_reach+pulse_length
          amp = pulse_shape_func(dist_to_transducer*wavelength)
        else
          amp = 0.0
        end
        image[x,y] += amp
      end
    end
  end
end

function ultrasim()
  wavelength = c/F0
  x_transducers = [pitch_transducers*itrans for itrans in range(1, n_transducers)]
  tvec = linspace(0, float(end_simulation_time), n_time_steps)
  image_pitch = image_fov ./ image_res
  image = zeros(Float32, (image_res[1], image_res[2]))

  for t in tvec
    simulate_one_time_step!(image, t, image_pitch, x_transducers, wavelength)

    figure(1)
    imshow(image)
    show()
  end

end

ultrasim()
