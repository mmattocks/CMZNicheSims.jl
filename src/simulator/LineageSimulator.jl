abstract type AbstractSSM_MC_run{P<:AbstractPopulation} end

mutable struct SSM_MC_run{P} <: AbstractSSM_MC_run{P}
    popsets::Vector{Vector{P}}
    run::Bool
end