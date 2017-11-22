# WaveSim.jl, wave propagation simulator.
# Space-time model defining when transducers fire.
# Christophe MEYER, 2016
using Lazy

# the space-time model
@defonce immutable TimingBoard
  # 
  clicked::AbstractMatrix
  # how much time does the full time axis represent
  timescale
end

# outer constructor
function newboard(m, n, timescale=30*1e-6)
  TimingBoard(fill(0, (m, n)), # start with nothing firing
              timescale)
end

# compute transducer delays from the model
# note: each one can fire only once
function delays(board)
  m, n = size(board.clicked)
  time_step = board.timescale / n
  # find at which time each transducer must fire
  delays = [time_step*(findfirst(board.clicked[i_trans,:])-1) for i_trans in 1:m]
end
