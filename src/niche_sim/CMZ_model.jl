struct CMZ_Record <: GMC_NS_Model_Record
    trajectory::Int64
    i::Int64
    pos::Vector{Float64}
    path::String
    log_Li::Float64
end

mutable struct CMZ_Model <: GMC_NS_Model
    trajectory::Int64
    i::Int64

    θ::Vector{Float64}
    log_Li::Float64

    pos::Vector{Float64}
    v::Vector{Float64}

    disp_mat::Array{Float64} #matrix for mean & 95%CI plot of model output

    function CMZ_Model(trajectory, i, θ, pos, v, obs, box, T, popdist, voldist, noise_dist, mc_its, phs; v_init=false)
        bind_phases!(θ,phs)
        pparams=θ[1:2*phs];phase_ends=θ[2*phs+1:(2*phs+1)+(phs-2)];volconst=θ[end]

        log_lh,disp_mat=CMZ_mc_llh(popdist, voldist, volconst, phase_ends, pparams, noise_dist, mc_its, T, obs)
    
        v_init && (v=rand(MvNormal(length(θ),1.)))

        new(trajectory, i, θ, log_lh, pos, v, disp_mat)
    end
end

construct_CMZ(args...;kwargs...) = CMZ_Model(args...;kwargs...)

function bind_phases!(θ,phs)
    if phs>2
        phase_ends=θ[2*phs+1:(2*phs+1)+(phs-2)]
        for pe in 1:phs-2
            phase_ends[pe]>phase_ends[pe+1]&& (phase_ends[pe+1]=phase_ends[pe]+1)
        end
    end
end
