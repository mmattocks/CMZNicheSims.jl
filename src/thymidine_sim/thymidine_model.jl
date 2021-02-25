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

function thymidine_constructor(trajectory, i, θ, pos, v, obs, popdist, T, pulse, mc_its, end_time; v_init=false)
    lt=length(θ)
    
    lt!=5 && mod(lt-5,6)!=0 && throw(ArgumentError("θ must contain tcμ, tcσ², g1_frac, s_frac, sis_frac values per lineage population, and one popfrac per population beyond the first!"))
    pulse<0 && throw(ArgumentError("pulse length must be >=0!"))

    n_pops=1
    lt>5 && (n_pops+=(lt-5)/6)

    pparams=Vector{Tuple{LogNormal, Float64, Float64, Float64}}()
    for pop in 1:n_pops
        tcμ, tcσ², g1_frac, s_frac, sis_frac= θ[Int64(1+((pop-1)*5)):Int64(5*pop)]
        tcσ=sqrt(tcσ²);
        push!(pparams,(LogNormal(tcμ,tcσ),g1_frac,s_frac,sis_frac))
    end

    n_pops > 1 ? (pop_fracs=θ[Int64(5*n_pops)+(pop-1):end]) : (    pop_fracs=[])

    log_lh,disp_mat=thymidine_mc_llh(popdist, pop_fracs, pparams, mc_its, end_time, pulse, T, obs)
    
    v_init && (v=rand(MvNormal(length(θ),1.)))

    Thymidine_Model(trajectory, i, θ, log_lh, pos, v, disp_mat)
end