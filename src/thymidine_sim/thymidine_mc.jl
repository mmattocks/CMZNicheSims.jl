
function init_thymidine_MC_run(pop_dists::Vector{Normal},refrac_times::Vector{Float64}, Tcdists::Vector{Normal}, s_fracs::Vector{Float64}, mc_its::Int64)
    n_pops=length(pop_dists)
    pops=Vector{Vector{Thymidine_Population}}(undef, mc_its)
    
    Threads.@threads for it in 1:mc_its
        pvec=Vector{Thymidine_Population}()
        for p in 1:n_pops
            pop_id=uuid4()

            n_lineages=Int64(max(1,round(rand(pop_dists[p]))))
            lineages=Vector{Lineage}(undef, n_lineages)
            progenitor_refrac=refrac_times[p]
            for l in 1:n_lineages
                lineage_id=uuid4()
                progenitor_id=uuid4()
                tδ=max(0.,rand(Tcdists[p]))+progenitor_refrac
                birth_time=rand(Uniform(-tδ,0.))

                cells=[SPhaseCell(progenitor_id,progenitor_id,true,[birth_time],[0.])]

                lineages[l]=Lineage(lineage_id,pop_id,cells)
            end
            push!(pvec,Thymidine_Population(pop_id,lineages,progenitor_refrac,params(Tcdists[p])...,s_fracs[p]))
        end
        pops[it]=pvec
    end
    return SSM_MC_run(pops,false)
end

function exec_t_MC_run!(smr, pop_params::Vector{Symbol}, end_time, m_params...)
    Threads.@threads for (n,popset) in collect(enumerate(smr.popsets))
        for (p,pop) in enumerate(popset)
            p_params=[getproperty(pop,param) for param in pop_params]
            for (i,l) in enumerate(pop.lineages)
                t_lineage_sim!(l, end_time, m_params..., p_params...)
            end
        end
    end
    smr.run=true
end

function t_lineage_sim!(l::AbstractLineage, end_time::AbstractFloat, params...)
    c_idx=1
    start=time()
    while c_idx <= length(l.cells)
        c=l.cells[c_idx]
        t=c.events[end]
        while c.mitotic && t <= end_time
            tδ,d=thymidine_cycle!(c,t,end_time, params...)
            t+=tδ
            t<= end_time && push!(l.cells,d)
        end
        c_idx+=1
    end
end