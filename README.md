# Ultrasound simulator

Simulates the propagation of waves from multiple transmitting elements.

## Usage

```
include("ultrasim.jl")
using UltraSim

# Setup parameters.
focus = 0.03  # [m]
steer = 0.0  # [deg]
aperture_size = 0.01  # [m]

# Compute focusing delays for the elements of the phased array.
trans_delays = UltraSim.delays_from_focus_and_steer(focus, steer, aperture_size)

# Run the simulation.
images = UltraSim.ultrasim(trans_delays)

# Display results.
include("view.jl")
imshowall(images)
```

Visualize the wave propagating through space, over time:

[TODO: add animated gif here]

Get a spatial heatmap of where the transmitted energy is sent:

[TODO: add picture here]

Know the transmit time delay everwhere in space:

[TODO: add picture here]

## Tips

### Parallelization

The code supports multi-threading, make use of it by setting:

`export JULIA_NUM_THREADS=4` (or whatever number of cores your machine has), before starting `julia`

or:

    JULIA_NUM_THREADS=`getconf _NPROCESSORS_ONLN` julia

