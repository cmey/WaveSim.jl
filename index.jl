using Colors
using Gadfly
using Compose

#### Model ####

include("ultrasim.jl")
include("timingmodel.jl")
board_size = (20, 10)

### Update ###

function update(board, click_coords)
  i, j = click_coords
  clicked = copy(board.clicked)
  should_remove = clicked[i, j] != 0
  fill!(sub(clicked, i, :), 0) # reset column
  if !should_remove
    clicked[i, j] = 1
  end
  new_board = TimingBoard(clicked, board.timescale)
  println(size(new_board.clicked), " ", new_board.clicked)
  return new_board
end

click_signal = Signal((0, 0))
initial_board_signal = Signal(TimingBoard,
	                          newboard(board_size[1], board_size[2]))
board_signal = flatten(
  map(initial_board_signal) do b
    foldp(update, b, click_signal; typ=TimingBoard)
  end
)

### View ###

colors = ["#fff", colormap("reds", 7)]

box(content, color) =
  inset(Escher.middle,
    fillcolor(color, size(1em, 1em, empty)),
    Escher.fontsize(2em, content)) |> paper(1) |> Escher.pad(0.2em)

number(x) = box(x == -1 ? "" : string(x) |> fontweight(800), colors[x+2])
clicked_icon = box("x", "#111")
unclicked_icon = box("", "#fff")
block(board, i, j) =
  return intent(constant((i, j)), clickable(
    board.clicked[i, j] != 0 ? clicked_icon : unclicked_icon)) >>> click_signal

 
global images = ultrasim(delays(newboard(board_size[1], board_size[2])))

function showboard(board::TimingBoard)
    m, n = size(board.clicked)
    b = vbox([hbox([block(board, i, j) for i in 1:m]) for j in 1:n])
    
    trans_delays = delays(board)
    global images = ultrasim(trans_delays)
    
    return b
end

function main(window)
  push!(window.assets, "widgets")

  eventloop = every(1/3)
  it = 0

  vbox(
    vskip(2em),
    title(3, "UltraSim.jl"),
    vskip(2em),
    map(showboard, board_signal, typ=Tile),
    map(eventloop) do _
      n_sim_time_steps = length(images)
      it = (it % n_sim_time_steps) + 1

  	  spy(images[it])
  	  # bitmap("test", rand(UInt8,9), 0, 1, 3, 3, :image)
    end
  ) |> packacross(center)
end
