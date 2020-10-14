mutable struct Thymidine_Population{L} <: AbstractPopulation{L}
    id::UUID
    lineages::Vector{L}
    refractory::Float64
    Tc_μ::Float64
    Tc_σ::Float64
    s_frac::Float64
end