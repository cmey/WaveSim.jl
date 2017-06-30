include("UltraSim.jl")
include("view.jl")

using UltraSim
images = UltraSim.main()

imshow4d(images)

