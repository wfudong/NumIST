module NumISTNext

using LinearAlgebra

include("kdv/spectral_etdrk4.jl")
include("kdv/paper_examples.jl")

include("ist/finite_section.jl")

include("crum/types.jl")
include("crum/utilities.jl")
include("crum/seed_trajectories.jl")
include("crum/one_fold.jl")
include("crum/two_fold.jl")
include("crum/algorithm.jl")

end
