using Colors
using Gadfly

#### Model ####

include("ultrasim.jl")
include("timingmodel.jl")

### Update ###

function update(board, click_coords)
  i, j = click_coords
  clicked = copy(board.clicked)
  fill!(sub(clicked, i, :), 0) # reset column
  clicked[i, j] = 1
  return TimingBoard(clicked, board.timescale)
end

click_signal = Signal((0, 0))
initial_board_signal = Signal(TimingBoard, newboard(20, 10))
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
    
function showboard(board::TimingBoard)
    m, n = size(board.clicked)
    b = vbox([hbox([block(board, i, j) for i in 1:m]) for j in 1:n])
    
    trans_delays = delays(board)
    images = ultrasim(trans_delays)
    vbox(
      b,
      # spy(rand(32,32))
      spy(images[end])
    )
end

function main(window)
  push!(window.assets, "widgets")

  #images = ultrasim()
  #n_sim_time_steps = length(images)

  # Run at 30 FPS
  # eventloop = every(1/5)

  it = 0
  # map(eventloop) do _
  	#it = (it % n_sim_time_steps) + 1
  	
    vbox(
      vskip(2em),
      title(3, "ultrasim"),
      vskip(2em),
      map(showboard, board_signal, typ=Tile),
    ) |> packacross(center)
  # end
end
