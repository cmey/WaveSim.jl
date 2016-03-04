# UltraSim.jl, ultrasound simulator.
# Web front-end using Escher.jl
# Christophe MEYER, 2016
using Colors
using Compose
using Escher
using Gadfly

#### Model ####

# the simulator
include("ultrasim.jl")
# the timing model for transducer firings
include("timingmodel.jl")
# size of the space-time user config of firings
board_size = (20, 10)

### Update functions ###

# update the space-time model based on user clicks
function update(board, click_coords)
  i, j = click_coords
  # toggle firings
  clicked = copy(board.clicked)
  should_remove = clicked[i, j] != 0
  fill!(sub(clicked, i, :), 0) # reset column
  if !should_remove
    clicked[i, j] = 1
  end
  new_board = TimingBoard(clicked, board.timescale)
  return new_board
end

# signal indicating user click
click_signal = Signal((0, 0))
# signal indicating initial state of space-time model
initial_board_signal = Signal(TimingBoard,
	                          newboard(board_size[1], board_size[2]))
# signal indicating space-time model was modified, and link to update function
board_signal = flatten(
  map(initial_board_signal) do b
    foldp(update, b, click_signal; typ=TimingBoard)
  end
)

### View (UI elements) ###

# a square box
box(content, color) =
  inset(Escher.middle,
    fillcolor(color, size(1em, 1em, empty)),
    Escher.fontsize(2em, content)) |> paper(1) |> Escher.pad(0.2em)

clicked_icon = box("x", "#111")
unclicked_icon = box("", "#fff")

# a clickable element of the board
block(board, i, j) =
  return intent(constant((i, j)), clickable(
    board.clicked[i, j] != 0 ? clicked_icon : unclicked_icon)) >>> click_signal

# the simulation results, shared
global images = ultrasim(delays(newboard(board_size[1], board_size[2])))

# render the new board, and run a new simulation
function showboard(board::TimingBoard)
    m, n = size(board.clicked)
    b = vbox([hbox([block(board, i, j) for i in 1:m]) for j in 1:n])
    
    trans_delays = delays(board)
    global images = ultrasim(trans_delays)
    
    return b
end

# web page render
function main(window)
  # we use clickable boxes
  push!(window.assets, "widgets")

  # frame rate to display the simulation results
  eventloop = every(1/3)
  it = 0

  # vertical align
  vbox(
    vskip(2em),
    title(3, "UltraSim.jl"),
    vskip(2em),
    # display the board
    map(showboard, board_signal, typ=Tile),
    # then display the results
    map(eventloop) do _
      n_sim_time_steps = length(images)
      it = (it % n_sim_time_steps) + 1

  	  spy(images[it])
  	  # bitmap("test", rand(UInt8,9), 0, 1, 3, 3, :image)
    end
  ) |> packacross(center)
end
