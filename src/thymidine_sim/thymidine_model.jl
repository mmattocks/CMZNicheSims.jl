struct Thymidine_Record <: GMC_NS_Model_Record
    path::String
    log_Li::Float64
end

mutable struct Thymidine_Model <: GMC_NS_Model
    id::Int64
    origin::Int64

    θ::Vector{Float64}
    v::Vector{Float64}
    log_Li::Float64
    fate_dist::Categorical

    retain_run::Bool #determines whether simulated MC run data is serialized or discarded
    run::SSM_MC_run

    disp_mat::Matrix{Float64} #matrix for mean & 95%CI plot of model output

    function Thymidine_Model(id, origin, θ, v, obs, T, pulse, mc_its, end_time, retain_run; v_init=false)
        mod(length(θ),6)!=0 && throw(ArgumentError("θ must contain 6 values per lineage population!"))
        pulse<0 && throw(ArgumentError("pulse length must be >=0!"))
        bound_θ!(θ)

        n_pops=length(θ)/6

        (pd, rt, tc, sf) = (Vector{Normal}(),Vector{Float64}(),Vector{Normal}(),Vector{Float64}())
        pvecs=[pd,rt,tc,sf]

        for pop in 1:n_pops
            pμ, pσ², r, tcμ, tcσ², s = θ[Int64(1+((pop-1)*6)):Int64(6*pop)]
            pσ=sqrt(pσ²); tcσ=sqrt(tcσ²)

            push!(pd,Normal(pμ,pσ)); push!(rt,r)
            push!(tc,Normal(tcμ,tcσ)); push!(sf,s)
        end

        fate_dist=Categorical([1.,0.,0.])

        smr=init_thymidine_MC_run(pvecs...,mc_its)
        exec_t_MC_run!(smr, [:refractory, :Tc_μ, :Tc_σ, :s_frac], end_time, pulse, fate_dist)
        log_lh,disp_mat=thymidine_ssm_mc_llh(smr,T,obs)

        !retain_run && (smr=SSM_MC_run(Vector{Vector{Population}}(),true))

        v_init && (v=rand(MvNormal(length(θ),1.)))

        new(id, origin, θ, v, log_lh, fate_dist, retain_run, smr, disp_mat)
    end
end
                function bound_θ!(θ)
                    npops=length(θ)/6
                    θ[findall(θi->θi<0., θ)].=nextfloat(0.)
                    for p in 1:npops
                        θ[Int64(3*p)]<1. && (θ[Int64(3*p)]=1.)
                        θ[Int64(4*p)]<1. && (θ[Int64(4*p)]=1.)
                        θ[Int64(6*p)]>1. && (θ[Int64(6*p)]=1.)
                    end
                    return θ
                end

thymidine_constructor(params...) = Thymidine_Model(params...)
