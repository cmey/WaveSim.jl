include("../src/benchmark.jl")
include("../src/save.jl")

images, sim_params = main();

beam_energy_map, transmit_time_map = WaveSim.beam_energy_map_and_transmit_time_map(images, sim_params);

saveall(images, beam_energy_map, transmit_time_map, sim_params, ".")
