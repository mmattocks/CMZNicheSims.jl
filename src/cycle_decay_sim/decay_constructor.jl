function construct_decay_slice(trajectory::Int64, i::Int64, θ::Vector{Float64}, pos::Vector{Float64}, v::Vector{Float64}, obs::Vector{Vector{Float64}}, T::Vector{Float64}, popdist::Distribution, lens_model::Lens_Model, mc_its::Int64; v_init=false)
    log_lh,disp_mat=decay_slice_mc_llh(popdist, lens_model, θ, mc_its, T, obs)

    v_init && (v=rand(MvNormal(length(θ),1.)))

    return Slice_Model(trajectory, i, θ, log_lh, pos, v, disp_mat)
end