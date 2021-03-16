struct MultiSlice_Model <: GMC_NS_Model
    trajectory::Int64
    i::Int64

    θ::Vector{Float64}
    log_Li::Float64

    pos::Vector{Float64}
    v::Vector{Float64}

    slices::Vector{Slice_Model}
end

function construct_multislice(trajectory::Int64, i::Int64, θ::Vector{Float64}, pos::Vector{Float64}, v::Vector{Float64}, obs::Vector{Vector{Vector{Float64}}}, T::Vector{Float64}, popdists::Vector{<:Distribution}, lens_model::Lens_Model, mc_its::Int64, phs::Int64; v_init=false)
    pparams=θ[1:2*phs];phase_ends=θ[2*phs+1:(2*phs+1)+(phs-2)]

    log_lh=0.

    slices=Vector{Slice_Model}()
    for (popdist,pobs) in zip(popdists, obs)
        llh,disp_mat=slice_mc_llh(popdist, lens_model, phase_ends, pparams, mc_its, T, pobs)
        log_lh+=llh

        push!(slices,Slice_Model(trajectory, i, θ, log_lh, pos, v, disp_mat))
    end

    v_init && (v=rand(MvNormal(length(θ),1.)))

    return MultiSlice_Model(trajectory, i, θ, log_lh, pos, v, slices)
end

function construct_decay_multislice(trajectory::Int64, i::Int64, θ::Vector{Float64}, pos::Vector{Float64}, v::Vector{Float64}, obs::Vector{Vector{Vector{Float64}}}, T::Vector{Float64}, popdists::Vector{<:Distribution}, lens_model::Lens_Model, mc_its::Int64; v_init=false)
    log_lh=0.

    slices=Vector{Slice_Model}()
    for (popdist,pobs) in zip(popdists, obs)
        llh,disp_mat=decay_slice_mc_llh(popdist, lens_model, θ, mc_its, T, pobs)
        log_lh+=llh

        push!(slices,Slice_Model(trajectory, i, θ, log_lh, pos, v, disp_mat))
    end

    v_init && (v=rand(MvNormal(length(θ),1.)))

    return MultiSlice_Model(trajectory, i, θ, log_lh, pos, v, slices)
end