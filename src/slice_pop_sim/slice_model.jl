struct Slice_Record <: GMC_NS_Model_Record
    trajectory::Int64
    i::Int64
    pos::Vector{Float64}
    path::String
    log_Li::Float64
end

struct Slice_Model <: GMC_NS_Model
    trajectory::Int64
    i::Int64

    θ::Vector{Float64}
    log_Li::Float64

    pos::Vector{Float64}
    v::Vector{Float64}

    disp_mat::Array{Float64} #matrix for mean & 95%CI plot of model output

    function Slice_Model(trajectory::Int64, i::Int64, θ::Vector{Float64}, pos::Vector{Float64}, v::Vector{Float64}, obs::Vector{Tuple{Vector{Float64},Vector{Float64}}}, T::Vector{Float64}, popdist::Distribution, lens_model::Lens_Model, mc_its::Int64, phs::Int64; v_init=false)
        pparams=θ[1:2*phs];phase_ends=θ[2*phs+1:(2*phs+1)+(phs-2)]

        log_lh,disp_mat=slice_mc_llh(popdist, lens_model, phase_ends, pparams, mc_its, T, obs)
    
        v_init && (v=rand(MvNormal(length(θ),1.)))

        new(trajectory, i, θ, log_lh, pos, v, disp_mat)
    end
end

construct_slice(args...;kwargs...) = Slice_Model(args...;kwargs...)
