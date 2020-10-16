function g_thymidine_ssm_mc_llh(mcr::SSM_MC_run, T::Vector{<:AbstractFloat}, obs::Vector{<:AbstractVector{<:Integer}})
    times=length(T); its=length(mcr.popsets)
    pops=length(mcr.popsets[1])

    popset_labelled=zeros(its,pops,times)
    #Threads.@threads 
    for (it,popset) in collect(enumerate(mcr.popsets))
        for (n,pop) in enumerate(popset)
            for lineage in pop.lineages
                popset_labelled[it,n,:]+=g_label_ct_at_time(lineage,T)
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

function g_label_ct_at_time(l::GraphLineage, times::Vector{Float64})
    g=l.graph
    cts=zeros(Int64,length(times)) #potentially no cells are labelled

    for (nt,t) in enumerate(times)
        t_mitoses=neighborhood(g,1,t)
        for m in t_mitoses
            if get_prop(g,m,:label) > 0.
                for nb in outneighbors(g,m)
                    !(nb in t_mitoses) && (cts[nt]+=1)
                end
            end
        end
    end

    return cts
end