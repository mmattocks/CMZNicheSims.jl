abstract type AbstractCell end

mutable struct LabelCell <: AbstractCell 
    time::Float64
    tδ′::Float64
    g1′::Float64
    s′::Float64
    label::Float64
end
