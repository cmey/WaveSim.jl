include("../src/benchmark.jl")
include("../src/save.jl")

sim_params, images = main();

beam_energy_map, transmit_time_map = WaveSim.beam_energy_map_and_transmit_time_map(images);

saveall(sim_params, images, beam_energy_map, transmit_time_map, ".")
