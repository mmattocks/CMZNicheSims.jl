mutable struct Thymidine_Record <: GMC_NS_Model_Record
    trajectory::Int64
    i::Int64
    pos::Vector{Float64}
    path::String
    log_Li::Float64
end

struct Thymidine_Model <: GMC_NS_Model
    trajectory::Int64
    i::Int64

    θ::Vector{Float64}
    log_Li::Float64

    pos::Vector{Float64}
    v::Vector{Float64}

    disp_mat::Matrix{Float64} #matrix for mean & 95%CI plot of model output
end

function thymidine_constructor(trajectory, i, θ, pos, v, obs, T, pulse, mc_its, end_time; v_init=false)
    mod(length(θ),7)!=0 && throw(ArgumentError("θ must contain 7 values per lineage population!"))
    pulse<0 && throw(ArgumentError("pulse length must be >=0!"))

    n_pops=length(θ)/7
    pparams=Vector{Tuple{LogNormal, LogNormal, Float64, Float64, Float64}}()
    for pop in 1:n_pops
        lpμ, lpσ², tcμ, tcσ², g1_frac, s_frac, sis_frac= θ[Int64(1+((pop-1)*7)):Int64(7*pop)]        
        lpσ=sqrt(lpσ²); tcσ=sqrt(tcσ²);
        push!(pparams,(LogNormal(lpμ,lpσ),LogNormal(tcμ,tcσ),g1_frac,s_frac,sis_frac))
    end

    log_lh,disp_mat=thymidine_mc_llh(pparams, mc_its, end_time, pulse, T, obs)
    
    v_init && (v=rand(MvNormal(length(θ),1.)))

    Thymidine_Model(trajectory, i, θ, log_lh, pos, v, disp_mat)
end