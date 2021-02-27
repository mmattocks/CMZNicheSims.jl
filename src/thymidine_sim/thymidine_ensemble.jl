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
    constants, #popdist, T, pulse, mc_its, end_time
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
    T=e.constants[2]

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

function print_MAP_output(e::Thymidine_Ensemble, path::String=e.path)
    MLEmod=deserialize(e.models[findmax([m.log_Li for m in e.models])[2]].path)
    θ=MLEmod.θ

    n_pops=1
    lt>5 && (n_pops+=(lt-5)/6)

    pparams=Vector{Tuple{LogNormal, Float64, Float64, Float64}}()
    for pop in 1:n_pops
        tcμ, tcσ², g1_frac, s_frac, sis_frac= θ[Int64(1+((pop-1)*5)):Int64(5*pop)]        
        tcσ=sqrt(tcσ²);
        push!(pparams,(LogNormal(tcμ,tcσ),g1_frac,s_frac,sis_frac))
    end

    n_pops > 1 ? (pop_fracs=θ[Int64(5*n_pops)+(pop-1):end]) : (    pop_fracs=[])

    popdist, T, pulse, mc_its, end_time = e.constants

    joint_DNPs=thymidine_mc_llh(popdist, pop_fracs, pparams, Int(1e6), end_time, pulse, T, obs)
    
    max_ct=0
    for DNP in joint_DNPs
        max_ct<maximum(DNP.support) && (max_ct=maximum(DNP.support))
    end
    cts=LinRange(0,max_ct,max_ct+1)
    probs=zeros(length(T),length(cts))
    for (t,DNP) in enumerate(joint_DNPs)
        for (ct, p) in zip(DNP.support,DNP.p)
            ct_idx=findfirst(isequal(ct),cts)
            ct_idx !== nothing && (probs[t,ct_idx]=p)
        end
    end
    catobs=vcat(e.obs...)
    ymax=Int64(round(maximum(catobs)*1.1))
    ys=[0,ymax]
    Ts=vcat([[T[n] for i in 1:length(e.obs[n])] for n in 1:length(T)]...)

    MAPout=Plots.heatmap(T,cts,transpose(probs),ylims=ys)
    Plots.scatter!(MAPout,Ts,catobs,marker=:cross,markercolor=:green)

    savefig(MAPout, path)
end

function print_marginals(e::Thymidine_Ensemble, path::String;  param_names=["LogNormal Tc μ", "LogNormal Tc σ²", "G1 Fraction", "S Fraction", "Sister Shift Fraction"])
    wtvec=zeros(0)
    lt=length(e.priors)
    θvecs=[zeros(0) for param in 1:lt]
    for (n,rec) in enumerate(e.posterior_samples)
        m=deserialize(rec.path)
        for i in 1:length(θvecs)
            push!(θvecs[i],m.θ[i])
        end
        push!(wtvec,e.log_Liwi[n+1])
    end

    for rec in e.models
        m=deserialize(rec.path)
        for i in 1:length(θvecs)
            push!(θvecs[i],m.θ[i])
        end
        push!(wtvec,(lps(m.log_Li,lps(e.log_Xi[end],-log(length(e.models))))))
    end

    n_pops=1
    lt>5 && (n_pops+=(lt-5)/6)
    n_pops>1 ? (plotrows=6; og_pn=copy(param_names)) : (plotrows=5)


    for pop in 2:n_pops
        vcat(param_names,vcat(og_pn,"Population $pop Fraction"))
    end

    plots=Vector{Plots.Plot}()
    for (n,θvec) in enumerate(θvecs)
        p_label=param_names[n]
        occursin("Fraction",p_label) ? (xls=[0,1]) : (xls=[quantile(e.priors[n],.01),quantile(e.priors[n],.99)])


        θkde=kde(θvec,weights=exp.(wtvec))
        plt=plot(e.priors[n], color=:darkmagenta, fill=true, fillalpha=.5, label="Prior", xlabel=p_label, ylabel="Probability density", xlims=xls)
        plot!(θkde.x,θkde.density, color=:green, fill=true, fillalpha=.75, label="Posterior")
        push!(plots,plt)
        n_pops>1 && n==5 && push!(plots,plot())
    end

    combined=plot(plots..., layout=grid(plotrows,n_pops),size=(600*n_pops,1500))

    savefig(combined, path)
end