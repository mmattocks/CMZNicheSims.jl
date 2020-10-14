function thymidine_ssm_mc_llh(mcr::SSM_MC_run, T::Vector{<:AbstractFloat}, obs::Vector{<:AbstractVector{<:Integer}})
    times=length(T); its=length(mcr.popsets)
    pops=length(mcr.popsets[1])

    popset_labelled=zeros(its,pops,times)
    Threads.@threads for (it,popset) in collect(enumerate(mcr.popsets))
        for (n,pop) in enumerate(popset)
            for lineage in pop.lineages
                for (tn, t) in enumerate(T)
                    popset_labelled[it,n,tn]+=label_ct_at_time(lineage,t)
                end
            end
        end
    end

    pop_DNPs=[fit(DiscreteNonParametric,popset_labelled[:,p,t]) for p in 1:pops, t in 1:times]

    joints=[joint_DNP_sum(pop_DNPs[:,t]) for t in 1:times]

    log_lh=sum(sum.([logpdf(joints[t],obs[t]) for t in 1:times]))

    disp_mat=zeros(times,3)
    Threads.@threads for t in 1:times
        disp_mat[t,:]=[quantile(joints[t],.025),mean(joints[t]),quantile(joints[t],.975)]
    end

    return log_lh, disp_mat
end

                function label_ct_at_time(l::AbstractLineage,t::AbstractFloat)
                    ct=0
                    for c in l.cells
                        idx=findlast(tm->tm<t,c.events)
                        if !(isnothing(idx))
                            c.label[idx] > 0. && (ct+=1)
                        end
                    end
                    return ct
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
    
