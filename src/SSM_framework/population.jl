abstract type AbstractPopulation{L} end

struct Population{L} <: AbstractPopulation{L}
    id::UUID
    lineages::Vector{L}
end

function size(p::Population)
    return length(p.lineages),sum([length(lineage) for lineage in p.lineages])
end