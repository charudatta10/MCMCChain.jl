module MCMCChain

import Showoff: showoff
import StatsBase: autocor, autocov, countmap, counts, describe, predict,
       quantile, sample, sem, summarystats
import LinearAlgebra: diag

using RecipesBase
import RecipesBase: plot

using Distributions
using SpecialFunctions
using AxisArrays

export Chain
export plot, traceplot, meanplot, densityplot, histogramplot, mixeddensityplot, autcorplot
export describe

# export diagnostics functions
export discretediag, gelmandiag, gewekediag, heideldiag, rafterydiag

abstract type AbstractChain end

"""
    Chains type

Parameters:

- `value`: `iterations × variables × chains` Data array
"""
struct Chain{T<:Real} <: AbstractChain
    values::AxisArray
    uniquenames::Dict{Symbol,Int}
    logevidence::Vector{T}
    weights::Matrix{T}
end

# imports
#include("utils.jl")

include("chains.jl")
#include("chainsummary.jl")
#include("discretediag.jl")
#include("fileio.jl")
#include("gelmandiag.jl")
#include("gewekediag.jl")
#include("heideldiag.jl")
#include("mcse.jl")
#include("rafterydiag.jl")
#include("stats.jl")
#include("plot.jl")

end # module
