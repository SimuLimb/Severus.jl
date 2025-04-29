module Severus

using Comodo
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Rotations
using Comodo.Statistics
using Geogram
using FEBio 

include("functions.jl")

# Export imported modules for later possible use
export Comodo
export Geogram
export FEBio

# Export functions
export severusdir, smoothboundary, taperfunction, bezierclose, set_z_level_boundary!, getboundaryset

end # module Severus
