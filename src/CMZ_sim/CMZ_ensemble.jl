mutable struct CMZ_Ensemble <: GMC_NS_Ensemble
    path::String

    model_initλ::Function
    models::Vector{CMZ_Record}

    contour::Float64
    log_Li::Vector{Float64}
    log_Xi::Vector{Float64}
    log_wi::Vector{Float64}
    log_Liwi::Vector{Float64}
    log_Zi::Vector{Float64}
    Hi::Vector{Float64}

    obs::AbstractVector{<:Tuple{<:AbstractVector{<:Float64},<:AbstractVector{<:Float64}}}
    priors::Vector{<:Distribution}
    constants::Vector{<:Any}
    #T, popdist, voldist, mc_its, phs
    box::Matrix{Float64}

    sample_posterior::Bool
    posterior_samples::Vector{CMZ_Record}

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

CMZ_Ensemble(path::String, no_models::Integer, obs::AbstractVector{<:Tuple{<:AbstractVector{<:Float64},<:AbstractVector{<:Float64}}}, priors::AbstractVector{<:Distribution}, constants, box, GS_settings; sample_posterior::Bool=true) =
CMZ_Ensemble(
    path,
    construct_CMZ,
    assemble_CMs(path, no_models, obs, priors, constants, box)...,
    [-Inf], #L0 = 0
	[0], #X0 = 1
	[-Inf], #w0 = 0
	[-Inf], #Liwi0 = 0
	[-Inf], #Z0 = 0
	[0], #H0 = 0,
    obs,
    priors,
    constants, #T, popdist, voldist, mc_its, phs
    box,
    sample_posterior,
    Vector{CMZ_Record}(),
    GS_settings...,
    no_models+1)

function assemble_CMs(path::String, no_trajectories::Integer, obs, priors, constants, box)
    ensemble_records = Vector{CMZ_Record}()
    !isdir(path) && mkpath(path)
    phs=constants[6]
    @showprogress 1 "Assembling CMZ_Model ensemble..." for trajectory_no in 1:no_trajectories
        model_path = string(path,'/',trajectory_no,'.',1)
        if !isfile(model_path)
            proposal=rand.(priors)
            proposal[2*phs+1:(2*phs+1)+(phs-2)]=sort(proposal[2*phs+1:(2*phs+1)+(phs-2)])

            pos=to_unit_ball.(proposal,priors)
            box_bound!(pos,box)
            θvec=to_prior.(pos,priors)
            
            model = construct_CMZ(trajectory_no, 1, θvec, pos, [0.], obs, constants...; v_init=true)
    
            serialize(model_path, model) #save the model to the ensemble directory
            push!(ensemble_records, CMZ_Record(trajectory_no, 1, pos,model_path,model.log_Li))
        else #interrupted assembly pick up from where we left off
            model = deserialize(model_path)
            push!(ensemble_records, CMZ_Record(trajectory_no, 1, model.pos,model_path,model.log_Li))
        end
    end

    return ensemble_records, minimum([record.log_Li for record in ensemble_records])
end


function Base.show(io::IO, m::CMZ_Model, e::CMZ_Ensemble; progress=false)
    T=e.constants[1]

    catpobs=vcat([e.obs[t][1] for t in 1:length(T)]...)
    ymax=max(maximum(m.disp_mat[:,:,1]),maximum(catpobs))
    ymin=min(minimum(m.disp_mat[:,:,1]),minimum(catpobs))

    plt=lineplot(T,m.disp_mat[:,2,1],title="CMZ_Model $(m.trajectory).$(m.i), log_Li $(m.log_Li)",color=:green,name="μ pop", ylim=[ymin,ymax])
    lineplot!(plt,T,m.disp_mat[:,1,1],color=:magenta,name="95% CI")
    lineplot!(plt,T,m.disp_mat[:,3,1],color=:magenta)

    Ts=vcat([[T[n] for i in 1:length(e.obs[n][1])] for n in 1:length(T)]...)
    scatterplot!(plt,Ts, catpobs, color=:yellow, name="Obs")

    show(io, plt)
    println()

    catvobs=vcat([e.obs[t][2] for t in 1:length(T)]...)
    ymax=max(maximum(m.disp_mat[:,:,1]),maximum(catvobs))
    ymin=min(minimum(m.disp_mat[:,:,1]),minimum(catvobs))

    plt=lineplot(T,m.disp_mat[:,2,2],color=:blue,name="μ vol", ylim=[ymin,ymax])
    lineplot!(plt,T,m.disp_mat[:,1,2],color=:yellow,name="95% CI")
    lineplot!(plt,T,m.disp_mat[:,3,2],color=:yellow)

    Ts=vcat([[T[n] for i in 1:length(e.obs[n][2])] for n in 1:length(T)]...)
    scatterplot!(plt,Ts, catvobs, color=:yellow, name="Obs")

    show(io, plt)

    println()
    println("θ: $(m.θ)")
    println("v: $(m.v)")

    (progress && return nrows(plt.graphics)+26);
end