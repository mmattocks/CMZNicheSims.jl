mutable struct Thymidine_Ensemble <: GMC_NS_Ensemble
    path::String

    model_initλ::Function
    models::Vector{Thymidine_Record}

    contour::Float64
    log_Li::Vector{Float64}
    log_Xi::Vector{Float64}
    log_wi::Vector{Float64}
    log_Liwi::Vector{Float64}
    log_Zi::Vector{Float64}
    Hi::Vector{Float64}

    obs::Vector{Vector{<:Integer}}
    priors::Vector{<:Distribution}
    constants::Vector{<:Any}
    #T, pulse, mc_its, end_time
    box::Matrix{Float64}

    sample_posterior::Bool
    posterior_samples::Vector{Thymidine_Record}

    GMC_Nmin::Int64

    GMC_τ_death::Float64
    GMC_init_τ::Float64
    GMC_tune_μ::Int64
    GMC_tune_α::Float64
    GMC_tune_PID::NTuple{3,Float64}

    GMC_timestep_η::Float64
    GMC_reflect_η::Float64
    GMC_exhaust_σ::Float64
    GMC_chain_κ::Int64

    t_counter::Int64
end

Thymidine_Ensemble(path::String, no_models::Integer, obs::AbstractVector{<:AbstractVector{<:Integer}}, priors::AbstractVector{<:Distribution}, constants, box, GMC_settings; sample_posterior::Bool=true) =
Thymidine_Ensemble(
    path,
    thymidine_constructor,
    assemble_TMs(path, no_models, obs, priors, constants, box)...,
    [-Inf], #L0 = 0
	[0], #ie exp(0) = all of the prior is covered
	[-Inf], #w0 = 0
	[-Inf], #Liwi0 = 0
	[-1e300], #Z0 = 0
	[0], #H0 = 0,
    obs,
    priors,
    constants, #T, pulse, mc_its, end_time
    box,
    sample_posterior,
    Vector{Thymidine_Record}(),
    GMC_settings...,
	no_models+1)

function assemble_TMs(path::String, no_trajectories::Integer, obs, priors, constants, box)
	ensemble_records = Vector{Thymidine_Record}()
	!isdir(path) && mkpath(path)
    @showprogress 1 "Assembling Thymidine_Model ensemble..." for trajectory_no in 1:no_trajectories
		model_path = string(path,'/',trajectory_no,'.',1)
        if !isfile(model_path)
            proposal=rand.(priors)

            pos=to_unit_ball.(proposal,priors)
            box_bound!(pos,box)
            θvec=to_prior.(pos,priors)

            model = thymidine_constructor(trajectory_no, 1, θvec, pos, [0.], obs, constants...; v_init=true)

			serialize(model_path, model) #save the model to the ensemble directory
			push!(ensemble_records, Thymidine_Record(trajectory_no, 1, pos,model_path,model.log_Li))
		else #interrupted assembly pick up from where we left off
			model = deserialize(model_path)
			push!(ensemble_records, Thymidine_Record(trajectory_no, 1, model.pos,model_path,model.log_Li))
		end
	end

	return ensemble_records, minimum([record.log_Li for record in ensemble_records])
end

function Base.show(io::IO, m::Thymidine_Model, e::Thymidine_Ensemble; progress=false)
    T=e.constants[1]

    catobs=vcat(e.obs...)
    ymax=max(maximum(m.disp_mat),maximum(catobs))
    ymin=min(minimum(m.disp_mat),minimum(catobs))

    plt=lineplot(T,m.disp_mat[:,2],title="Thymidine_Model $(m.trajectory).$(m.i), log_Li $(m.log_Li)",color=:green,name="μ count", ylim=[ymin,ymax])
    lineplot!(plt,T,m.disp_mat[:,1],color=:magenta,name="95% CI")
    lineplot!(plt,T,m.disp_mat[:,3],color=:magenta)

    Ts=vcat([[T[n] for i in 1:length(e.obs[n])] for n in 1:length(T)]...)
    scatterplot!(plt,Ts, catobs, color=:yellow, name="Obs")

    show(io, plt)
    println()
    println("θ: $(m.θ)")

    (progress && return nrows(plt.graphics)+7);
end
