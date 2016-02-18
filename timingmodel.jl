#### Timing Model ####

using Lazy

@defonce immutable TimingBoard
  clicked::AbstractMatrix
  timescale
end

function newboard(m, n, timescale=30*1e-6)
  TimingBoard(fill(0, (m, n)),
              timescale)
end

function delays(board)
  m, n = size(board.clicked)
  time_step = board.timescale / n
  delays = [time_step*(findfirst(board.clicked[i_trans,:])-1) for i_trans in 1:m]
end
