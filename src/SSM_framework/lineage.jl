abstract type AbstractLineage{C} end

mutable struct Lineage{C} <: AbstractLineage{C}
    id::UUID
    population::UUID
    cells::Vector{C}
end

function Base.size(l::Lineage)
    return length(l.cells)
end

function Base.size(l::Lineage, t::AbstractFloat)
    return count([c.birthtime <= t for c in l.cells])
end
