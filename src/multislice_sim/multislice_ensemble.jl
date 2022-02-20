
mutable struct MultiSlice_Ensemble <: GMC_NS_Ensemble
    path::String

    model_initλ::Function
    models::Vector{Slice_Record}

    contour::Float64
    log_Li::Vector{Float64}
    log_Xi::Vector{Float64}
    log_wi::Vector{Float64}
    log_Liwi::Vector{Float64}
    log_Zi::Vector{Float64}
    Hi::Vector{Float64}

    obs::AbstractVector{<:AbstractVector{<:AbstractVector{<:Float64}}}
    priors::Vector{<:Distribution}
    constants::Vector{<:Any}
    #T, popdist, lens_model, mc_its, phs
    box::Matrix{Float64}

    sample_posterior::Bool
    posterior_samples::Vector{Slice_Record}

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

MultiSlice_Ensemble(path::String, no_models::Integer, obs::AbstractVector{<:AbstractVector{<:AbstractVector{<:Float64}}}, priors::AbstractVector{<:Distribution}, constants, box, GS_settings; sample_posterior::Bool=true) =
MultiSlice_Ensemble(
    path,
    construct_multislice,
    assemble_MSMs(path, no_models, obs, priors, constants, box)...,
    [-Inf], #L0 = 0
	[0], #X0 = 1
	[-Inf], #w0 = 0
	[-Inf], #Liwi0 = 0
	[-Inf], #Z0 = 0
	[0], #H0 = 0,
    obs,
    priors,
    constants, #T, popdist, lens_model, mc_its, phs
    box,
    sample_posterior,
    Vector{Slice_Record}(),
    GS_settings...,
    no_models+1)

MultiSlice_Decay_Ensemble(path::String, no_models::Integer, obs::AbstractVector{<:AbstractVector{<:AbstractVector{<:Float64}}}, priors::AbstractVector{<:Distribution}, constants, box, GS_settings; sample_posterior::Bool=true) =
MultiSlice_Ensemble(
        path,
        construct_decay_multislice,
        assemble_MSDMs(path, no_models, obs, priors, constants, box)...,
        [-Inf], #L0 = 0
        [0], #X0 = 1
        [-Inf], #w0 = 0
        [-Inf], #Liwi0 = 0
        [-Inf], #Z0 = 0
        [0], #H0 = 0,
        obs,
        priors,
        constants, #T, popdist, lens_model, mc_its
        box,
        sample_posterior,
        Vector{Slice_Record}(),
        GS_settings...,
        no_models+1)

function assemble_MSMs(path::String, no_trajectories::Integer, obs, priors, constants, box)
    ensemble_records = Vector{Slice_Record}()
    !isdir(path) && mkpath(path)
    phs=constants[5]
    @showprogress 1 "Assembling Slice_Model ensemble..." for trajectory_no in 1:no_trajectories
        model_path = string(path,'/',trajectory_no,'.',1)
        if !isfile(model_path)
            proposal=rand.(priors)
            proposal[2*phs+1:(2*phs+1)+(phs-2)]=sort(proposal[2*phs+1:(2*phs+1)+(phs-2)])

            pos=to_unit_ball.(proposal,priors)
            box_bound!(pos,box)
            θvec=to_prior.(pos,priors)
            
            model = construct_multislice(trajectory_no, 1, θvec, pos, [0.], obs, constants...; v_init=true)
    
            serialize(model_path, model) #save the model to the ensemble directory
            push!(ensemble_records, Slice_Record(trajectory_no, 1, pos,model_path,model.log_Li))
        else #interrupted assembly pick up from where we left off
            model = deserialize(model_path)
            push!(ensemble_records, Slice_Record(trajectory_no, 1, model.pos,model_path,model.log_Li))
        end
    end

    return ensemble_records, minimum([record.log_Li for record in ensemble_records])
end

function assemble_MSDMs(path::String, no_trajectories::Integer, obs, priors, constants, box)
    ensemble_records = Vector{Slice_Record}()
    !isdir(path) && mkpath(path)
    @showprogress 1 "Assembling Slice_Model ensemble..." for trajectory_no in 1:no_trajectories
        model_path = string(path,'/',trajectory_no,'.',1)
        if !isfile(model_path)
            proposal=rand.(priors)
            pos=to_unit_ball.(proposal,priors)
            box_bound!(pos,box)
            θvec=to_prior.(pos,priors)
            
            model = construct_decay_multislice(trajectory_no, 1, θvec, pos, [0.], obs, constants...; v_init=true)
    
            serialize(model_path, model) #save the model to the ensemble directory
            push!(ensemble_records, Slice_Record(trajectory_no, 1, pos,model_path,model.log_Li))
        else #interrupted assembly pick up from where we left off
            model = deserialize(model_path)
            push!(ensemble_records, Slice_Record(trajectory_no, 1, model.pos,model_path,model.log_Li))
        end
    end

    return ensemble_records, minimum([record.log_Li for record in ensemble_records])
end

function Base.show(io::IO, m::MultiSlice_Model, e::MultiSlice_Ensemble; progress=false)
    T=e.constants[1]

    println("MultiSlice_Model $(m.trajectory).$(m.i)")

    rows=0
    for (n,(obs, slice)) in enumerate(zip(e.obs,m.slices))
        catpobs=vcat([obs[t] for t in 1:length(T)]...)
        ymax=max(maximum(slice.disp_mat[:,:]),maximum(catpobs))
        ymin=min(minimum(slice.disp_mat[:,:]),minimum(catpobs))

        plt=lineplot(T,slice.disp_mat[:,2],title="Slice $n",color=:green,name="μ pop", ylim=[ymin,ymax])
        lineplot!(plt,T,slice.disp_mat[:,1],color=:magenta,name="95% CI")
        lineplot!(plt,T,slice.disp_mat[:,3],color=:magenta)

        Ts=vcat([[T[n] for i in 1:length(obs[n])] for n in 1:length(T)]...)
        scatterplot!(plt,Ts, catpobs, color=:yellow, name="Obs")

        rows+=nrows(plt.graphics)
        show(io, plt)
        println()
    end

    println()
    println("log_Li: $(m.log_Li)")
    println("θ: $(m.θ)")
    println("v: $(m.v)")

    progress ? (return rows+14) : (return)
end
