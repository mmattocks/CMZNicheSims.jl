Decay_Ensemble(path::String, no_models::Integer, obs::AbstractVector{<:AbstractVector{<:Float64}}, priors::AbstractVector{<:Distribution}, constants, box, GMC_settings; sample_posterior::Bool=true) =
Slice_Ensemble(
    path,
    construct_decay_slice,
    assemble_DMs(path, no_models, obs, priors, constants, box)...,
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
    GMC_settings...,
    no_models+1)

function assemble_DMs(path::String, no_trajectories::Integer, obs, priors, constants, box)
    ensemble_records = Vector{Slice_Record}()
    !isdir(path) && mkpath(path)
    @showprogress 1 "Assembling Slice_Model ensemble..." for trajectory_no in 1:no_trajectories
        model_path = string(path,'/',trajectory_no,'.',1)
        if !isfile(model_path)
            proposal=rand.(priors)
            pos=to_unit_ball.(proposal,priors)
            box_bound!(pos,box)
            θvec=to_prior.(pos,priors)
            
            model = construct_decay_slice(trajectory_no, 1, θvec, pos, [0.], obs, constants...; v_init=true)
    
            serialize(model_path, model) #save the model to the ensemble directory
            push!(ensemble_records, Slice_Record(trajectory_no, 1, pos,model_path,model.log_Li))
        else #interrupted assembly pick up from where we left off
            model = deserialize(model_path)
            push!(ensemble_records, Slice_Record(trajectory_no, 1, model.pos,model_path,model.log_Li))
        end
    end

    return ensemble_records, minimum([record.log_Li for record in ensemble_records])
end