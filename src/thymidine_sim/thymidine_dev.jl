function d_Thymidine_Model(trajectory, i, θ, v, obs, T, pulse, mc_its, end_time, retain_run; v_init=false)
    mod(length(θ),8)!=0 && throw(ArgumentError("θ must contain 8 values per lineage population!"))
    pulse<0 && throw(ArgumentError("pulse length must be >=0!"))
    bound_θ!(θ)

    n_pops=length(θ)/8
    pparams=Vector{Tuple{LogNormal, LogNormal, Float64, Normal, Float64}}()
    for pop in 1:n_pops
        pμ, pσ², r, tcμ, tcσ², sμ, sσ², sis_frac= θ[Int64(1+((pop-1)*8)):Int64(8*pop)]
        pσ=sqrt(pσ²); tcσ=sqrt(tcσ²); sσ=sqrt(sσ²)

        push!(pparams,(LogNormal(pμ,pσ),LogNormal(tcμ,tcσ),r,Normal(sμ,sσ),sis_frac))
    end

    fate_dist=Categorical([1.,0.,0.])

    log_lh,disp_mat=d_thymidine_ssm_mc_llh(pparams, mc_its, end_time, pulse, T, obs)
    
    !retain_run && (smr=SSM_MC_run(Vector{Vector{Population}}(),true))

    v_init && (v=rand(MvNormal(length(θ),1.)))

    Thymidine_Model(trajectory, i, θ, v, log_lh, fate_dist, retain_run, smr, disp_mat)
end

function d_thymidine_ssm_mc_llh(pparams, mc_its, end_time, pulse, T::Vector{<:AbstractFloat}, obs::Vector{<:AbstractVector{<:Integer}})
    n_pops=length(pparams); times=length(T)
    popset_labelled=zeros(Int64,mc_its,n_pops,times)
    
    Threads.@threads for it in 1:mc_its
        for p in 1:n_pops
            pop_dist,tc_dist,r,s,sis_frac = pparams[p]
            n_lineages=Int64(max(1,round(rand(pop_dist))))
            plv=view(popset_labelled,it,p,:)

            for l in 1:Int64(floor(n_lineages/2))
                d_sim_lin_pair!(plv, T, end_time, pulse, tc_dist, r, s, sis_frac)
            end

            for l in 1:n_lineages%2
                d_sim_lineage!(plv, T, end_time, pulse, tc_dist, r, s, sis_frac)
            end
        end
    end

    pop_DNPs=Matrix{DiscreteNonParametric}(undef,n_pops,times)
    for p in 1:n_pops
        Threads.@threads for t in 1:times
            pop_DNPs[p,t]=fit(DiscreteNonParametric,popset_labelled[:,p,t])
        end
    end

    joints=Vector{DiscreteNonParametric}(undef,times)    
    Threads.@threads for t in 1:times
        joints[t]=joint_DNP_sum(pop_DNPs[:,t])
    end

    log_lhs=Vector{Float64}(undef,times)
    Threads.@threads for t in 1:times
        log_lhs[t]=lps(logpdf(joints[t],obs[t]))
    end

    log_lh=lps(log_lhs)

    disp_mat=zeros(times,3)
    Threads.@threads for t in 1:times
        disp_mat[t,:]=[quantile(joints[t],.025),mean(joints[t]),quantile(joints[t],.975)]
    end
    return log_lh, disp_mat
end

function d_sim_lin_pair!(plv, T, end_time, pulse, Tc, r, s, sis_frac)
    tδ1, s1 = d_refractory_cycle_model(r, Tc, s)
    tδ2, s2 = (1+rand([-1,1])*sis_frac).*(tδ1, s1)

    birth_time=rand(Uniform(-min(tδ1,tδ2),0.))

    
    cells1=[LabelCell(birth_time,tδ1,s1,0.)]
    cells2=[LabelCell(birth_time,tδ2,s2,0.)]

    for cellvec in [cells1,cells2]
        ci=1
        while ci<=length(cellvec)
            cycle_sim!(plv, ci, cellvec, T, end_time, pulse, Tc,r,s, sis_frac)
            ci+=1
        end
    end    
end


function d_sim_lineage!(plv, T, end_time, pulse, Tc, r, s, sis_frac)
    tδ1, s1 = d_refractory_cycle_model(r, Tc, s)
    birth_time=rand(Uniform(-tδ1,0.))
    cells=[LabelCell(birth_time,tδ1,s1,0.)]

    ci=1
    while ci<= length(cells)
        cycle_sim!(plv, ci, cells, T, end_time, pulse, Tc, r, s, sis_frac)
        ci+=1
    end
end

function cycle_sim!(plv, ci::Int64, cells::Vector{LabelCell}, T::Vector{Float64}, end_time::Float64, pulse::Float64, Tc::LogNormal, r::Float64, sd::Normal, sis_frac::Float64)
    n_times=length(T);

    cell=cells[ci]
    t=cell.time; tδ′=cell.tδ′; s′=cell.s′; l=cell.label

    first_cycle=true
    while t<end_time
        tδ,s=0.,0.
        first_cycle ? (tδ=tδ′;s=s′;first_cycle=false) : ((tδ, s)=d_refractory_cycle_model(r, Tc, sd))

        s_start=t+r
        s_end=s_start+s
        cycle_mitosis=t+tδ
        cycle_events=[s_start,s_end,cycle_mitosis]

        oldl=l; cycle_label=false
        if (s_start <= pulse && s_end >= 0) #some s-phase overlaps with the pulse in this case
            cycle_label=true
            l=update_label(l, s, s_start, s_end, pulse)
        end

        tidxs=t.<T.<cycle_mitosis
        n_idxs=sum(tidxs)
        (length(tidxs)>0 && (oldl > DETECTION_THRESHOLD || l > DETECTION_THRESHOLD)) 
        
        if n_idxs>0 && (oldl > DETECTION_THRESHOLD || l > DETECTION_THRESHOLD)
            cycle_label ? (label_outcomes=[oldl,-Inf,l]) : (label_outcomes=[oldl,oldl,oldl])
            cubuf=falses(n_idxs)
            for (n,timept) in enumerate(T[tidxs])
                count_labelled_at_T!(n,cubuf,timept, cycle_events, label_outcomes, oldl, s_start, s, pulse)
            end
            plv[tidxs]+=cubuf
        end

        l/=2 #mitosis halves label for the daughters

        tδ2, s2 = (1+rand([-1,1])*sis_frac).*(tδ, s)

        cycle_mitosis < end_time && push!(cells,LabelCell(cycle_mitosis,tδ2,s2,l))
        t+=tδ
    end
end

function d_refractory_cycle_model(refractory, Tc, s_dist)
    cycle_time=rand(Tc)
    s=max(eps(),rand(s_dist))
    tδ=max(refractory + s, cycle_time)
    return tδ, s
end

function update_label(label, s, s_start, s_end, pulse)
        s_overlap=min(pulse,s_start+s_end)-max(0,s_start)
        nuc_label_frac=s_overlap/s
        return min(nuc_label_frac+label,1.) #label before mitosis time will be what's been accumulated from past + this cycle's s phase
end

function count_labelled_at_T!(n, cubuf, timept, cycle_events, label_outcomes, oldl, s_start, s, pulse)
    outcome_idx=timept.<cycle_events
    timept_label=label_outcomes[outcome_idx][1]
    timept_label==-Inf && (timept_label=partial_label(timept,oldl,s_start,s,pulse))
    timept_label >= DETECTION_THRESHOLD && (cubuf[n]=true)
end

function partial_label(timept, label, s_start, s, pulse)
    s_overlap=min(pulse,timept)-max(0.,s_start)
    nuc_label_frac=s_overlap/s
    return min(nuc_label_frac+label,1.)
end
