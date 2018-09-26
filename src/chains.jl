#################### Chains ####################


#################### Constructors ####################
"""
    Chain(value::Array{T,3})

Arguments:
- `value`: Samples obtained from a MCMC simulation. 
           Dimensions are: iterations x variables x chains.

Optional Arguments:
- `start = 1`
- `thin = 1`
- `names = String[]`
- `uniquenames = Dict{Symbol,Int}()`
- `chains = Int[]`
- `weights = Matrix{T}()`
- `logevidence = T[]`

"""
function Chain(
                value::Array{T, 3};
                start = 1,
                thin = 1,
                names = AbstractString[],
                uniquenames = Dict{Symbol, Int}(),
                chains = Int[],
                weights = Matrix{T}(undef, 0, 0),
                logevidence = T[]
               ) where {T<:Union{Real, Missing}}

    niters, nvars, nchains = size(value)

    axisnames = [:iteration, :variable, :chain]
    axisvals = [ 
                 range(start, step = thin, length = niters),
                 names,
                 map(i -> "Chain$i", chains)
    ]

    # Set values if necessary.
    if isempty(axisvals[2])
        axisvals[2] = map(i -> "Param$i", 1:nvars)
    end

    if isempty(axisvals[3])
        axisvals[3] = map(i -> "Chain$i", 1:nchains)
    end

    if isempty(logevidence)
        logevidence = zeros(T, nchains)
    end

    if isempty(weights)
        weights = ones(niters, nvars)
    end

    # Several sanity checks.
    @assert length(axisvals[2]) == nvars begin
        throw(DimensionMismatch("size(value,2) != length(names)"))
    end
    @assert length(axisvals[3]) == nchains begin
        throw(DimensionMismatch("size(value,3) != length(chains)"))
    end
    @assert length(logevidence) == nchains begin
        throw(DimensionMismatch("size(value,3) != length(logevidence)"))
    end
    @assert size(weights,1) == niters begin
        throw(DimensionMismatch("size(value,1) != size(weights,1)"))
    end
    @assert size(weights,2) == nchains begin
        throw(DimensionMismatch("size(value,2) != size(weights,2)"))
    end

    # Construction of Axis tuples.
    axs = ntuple(i -> Axis{axisnames[i]}(axisvals[i]), 3)
    
    # Set and prepare unique names fiels.
    if isempty(uniquenames)
        uniquenames = Dict(gensym() => i for i in 1:nvars)
    elseif length(uniquenames) != p
        throw(DimensionMismatch("size(value,2) != length(uniquenames)"))
    end

    A = AxisArray(convert(Array{Union{Missing, T}, 3}, value), axs...)
    return Chain{T}(A, uniquenames, logevidence, weights)
end

function Chain(
                iters::Int,
                params::Int;
                start = 1,
                thin = 1,
                chains = 1,
                names = AbstractString[],
                uniquenames = Dict{Symbol, Int}()
               )
    value = ones(length(start:thin:iters), params, chains)
    fill!(value, NaN)

    Chain(value, start=start, thin=thin, names=names, uniquenames = uniquenames)
end

function Chain(
                value::Matrix{T};
                start = 1,
                thin = 1,
                names = AbstractString[],
                uniquenames = Dict{Symbol, Int}(),
                chains = 1
               ) where {T<:Union{Real, Missing}}

    Chain(
        reshape(value, size(value, 1), size(value, 2), 1), 
        start=start,
        thin=thin,
        names=names,
        uniquenames = uniquenames,
        chains=Int[chains]
    )
end

function Chain(
                value::Vector{T};
                start = 1,
                thin = 1,
                names = "Param#1",
                uniquenames = gensym(),
                chains = 1
               ) where {T<:Union{Real, Missing}}
    Chain(
         reshape(value, length(value), 1, 1),
         start=start,
         thin=thin,
         names=AbstractString[names],
         uniquenames = Dict{Symbol, Int}(uniquenames => 1),
         chains=Int[chains]
    )
end

#################### Concatenation ####################
function Base.cat(chn::Chain{T}, chns::Chain{T}...; dims) where {T}
    
end


#################### Base Methods ####################
Base.getindex(c::Chain, i...) = getindex(c.values, i...)
Base.setindex!(c::Chain, v, i...) = setindex!(c.values, v, i...)
Base.keys(c::Chain) = keys(c.uniquenames)


function Base.show(io::IO, c::AbstractChain)
  print(io, "Object of type \"$(summary(c))\"\n\n")
  #println(io, header(c))
  show(io, c.values)
end

Base.size(c::Chain) = size(c.values)
Base.size(c::Chain, ind) = size(c)[ind]

Base.first(c::AbstractChain) = first(c.range)
Base.step(c::AbstractChain) = step(c.range)
Base.last(c::AbstractChain) = last(c.range)

#################### Auxilliary Functions ####################

function combine(c::AbstractChain)
  n, p, m = size(c.value)
  value = Array{Float64}(n * m, p)
  for j in 1:p
    idx = 1
    for i in 1:n, k in 1:m
      value[idx, j] = c.value[i, j, k]
      idx += 1
    end
  end
  value
end

function header(c::AbstractChain)
  string(
    "Iterations = $(first(c)):$(last(c))\n",
    "Thinning interval = $(step(c))\n",
    "Chains = $(join(map(string, c.chains), ","))\n",
    "Samples per chain = $(length(c.range))\n"
  )
end

function indiscretesupport(c::AbstractChain,
                           bounds::Tuple{Real, Real}=(0, Inf))
  nrows, nvars, nchains = size(c.value)
  result = Array{Bool}(undef, nvars * (nrows > 0))
  for i in 1:nvars
    result[i] = true
    for j in 1:nrows, k in 1:nchains
      x = c.value[j, i, k]
      if !isinteger(x) || x < bounds[1] || x > bounds[2]
        result[i] = false
        break
      end
    end
  end
  result
end

function link(c::AbstractChain)
  cc = copy(c.value)
  for j in 1:length(c.names)
    x = cc[:, j, :]
    if minimum(x) > 0.0
      cc[:, j, :] = maximum(x) < 1.0 ? logit.(x) : log.(x)
    end
  end
  cc
end
