function thymidine_mc_llh(pop_dist, pop_fracs, pparams, mc_its, end_time, pulse, T::Vector{<:AbstractFloat}, obs::Vector{<:AbstractVector{<:Integer}}, display=false)
    n_pops=length(pparams); times=length(T)
    popset_labelled=zeros(Int64,mc_its,n_pops,times)
    max_labelled=maximum.(obs)

    Threads.@threads for it in 1:mc_its
        pop_lins=zeros(Int64,n_pops)
        pop_lins[1]=Int64(max(1,round(rand(pop_dist))))
        for p in 2:n_pops
            pop_lins[p]=Int64(round(pop_fracs[p-1]*pop_lins[p-1]))
            pop_lins[p-1]-=pop_lins[p]
        end
        
        for p in 1:n_pops
            tc_dist,g1_frac,s_frac,sis_frac = pparams[p]
            plv=view(popset_labelled,it,p,:)

            for l in 1:Int64(floor(pop_lins[p]/2))
                !all(max_labelled .< sum(popset_labelled[it,:,:],dims=1)[1,:]) && sim_lin_pair!(plv, T, end_time, pulse, tc_dist, g1_frac, s_frac, sis_frac)
            end

            for l in 1:pop_lins[p]%2
                !all(max_labelled .< sum(popset_labelled[it,:,:],dims=1)[1,:]) && sim_lineage!(plv, T, end_time, pulse, tc_dist, g1_frac, s_frac, sis_frac)
            end
        end
    end

    pop_DNPs=Matrix{DiscreteNonParametric}(undef,n_pops,times)
    for p in 1:n_pops
        Threads.@threads for t in 1:times
            pop_DNPs[p,t]=fit(DiscreteNonParametric,popset_labelled[:,p,t])
        end
    end

    if n_pops > 1
        joints=Vector{DiscreteNonParametric}(undef,times)    
        Threads.@threads for t in 1:times
            joints[t]=joint_DNP_sum(pop_DNPs[:,t])
        end
    else
        joints=pop_DNPs[1,:]
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
    display ? (return log_lh, disp_mat, joints) : (return log_lh, disp_mat)
end

function joint_DNP_sum(cs)
    npops=length(cs)
    ct_supports=Vector{Vector{Int64}}()
    for c in cs
        push!(ct_supports,c.support[[c.p...].>0.])
    end

    possible_sums=collect(Iterators.product(ct_supports...))
    possible_totals=unique(sum.(possible_sums))
    total_mat=sum.(possible_sums)
    
    probs=zeros(length(possible_totals))

    Threads.@threads for (n,pt) in collect(enumerate(possible_totals))
        psums=possible_sums[findall(t->t==pt,total_mat)]
        probs[n]=logsumexp([lps([logpdf(cs[p],ps[p]) for p in 1:npops]) for ps in psums])
    end

    return DiscreteNonParametric(possible_totals, exp.(probs))
end


    
