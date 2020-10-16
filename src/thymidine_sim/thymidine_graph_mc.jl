function G_Thymidine_Model(id, origin, θ, v, obs, T, pulse, mc_its, end_time, retain_run; v_init=false)
    mod(length(θ),6)!=0 && throw(ArgumentError("θ must contain 6 values per lineage population!"))
    pulse<0 && throw(ArgumentError("pulse length must be >=0!"))
    bound_θ!(θ)

    n_pops=length(θ)/6

    (pd, pp) = (Vector{Normal}(),Vector{Vector{Float64}}())

    for pop in 1:n_pops
        pμ, pσ², r, tcμ, tcσ², s = θ[Int64(1+((pop-1)*6)):Int64(6*pop)]
        pσ=sqrt(pσ²); tcσ=sqrt(tcσ²)

        push!(pd, Normal(pμ,pσ))
        push!(pp, [r, tcμ, tcσ, s])
    end


    smr=graphrun(pulse, end_time, pd, pp, mc_its)
    log_lh,disp_mat=g_thymidine_ssm_mc_llh(smr,T,obs)

    !retain_run && (smr=SSM_MC_run(Vector{Vector{Population}}(),true))

    v_init && (v=rand(MvNormal(length(θ),1.)))

    Thymidine_Model(id, origin, θ, v, log_lh, Categorical([1.]), retain_run, smr, disp_mat)
end

function graphrun(pulse::Float64, end_time::Float64, pop_dists::Vector{Normal}, pop_params::Vector{Vector{Float64}}, mc_its::Int64)
    n_pops=length(pop_dists)
    pops=Vector{Vector{Thymidine_Population}}(undef, mc_its)
    
    Threads.@threads for it in 1:mc_its
        pvec=Vector{Thymidine_Population}()
        for p in 1:n_pops
            pop_id=uuid4(); pparams=pop_params[p]
                        
            n_lineages=Int64(max(1,round(rand(pop_dists[p]))))
            lineages=Vector{GraphLineage}(undef, n_lineages)
            for l in 1:n_lineages
                lineages[l]=graphsim(pulse, end_time, pparams...)
            end
            push!(pvec,Thymidine_Population(pop_id,lineages,pparams...))
        end
        pops[it]=pvec
    end
    return SSM_MC_run(pops,true)
end

#[refractory][S][Tc---M]
function graphsim(pulse::Float64, end_time::Float64, refractory, Tcμ, Tcσ, s_frac)
    l,mts=initialize_lineage(pulse, end_time, refractory, Tcμ, Tcσ, s_frac)
    g=l.graph
    
    while length(mts) > 0
        #determine which mitosis vertex to use and what the time at that vertex is
        m=popfirst!(mts)
        mt=time_to_m(l,m)
        #determine the mitosis vertex indices to be added and add them
        v=length(vertices(g)); d1idx=v+1; d2idx=v+2
        add_vertices!(g,2); add_edge!(g,m,d1idx); add_edge!(g,m,d2idx)
        #calculate cycle times and assign the appropriate attributes to edges and vertices
        m_label=get_prop(g,m,:label)
        for idx in [d1idx,d2idx]
            tδ,s=refractory_cycle_model(refractory, Tcμ, Tcσ, s_frac)
            set_prop!(g,m,idx, :weight, tδ)
            lbl=check_cell_label(mt, refractory, s, pulse, m_label)
            set_prop!(g,idx,:label,lbl)
            mt+tδ<end_time && push!(mts,idx) #push idx to mitoses to sim if mitotic event would occur before sim end
        end
    end
    
    return l
end

function initialize_lineage(pulse, end_time, refractory, Tcμ, Tcσ, s_frac)
    g=MetaDiGraph(SimpleDiGraph(2)) #instantiate lineage with 1 founder; edge is its first cycle
    add_edge!(g,1,2)
    tδ,s=refractory_cycle_model(refractory, Tcμ, Tcσ, s_frac)
    set_prop!(g,1,2, :weight, tδ)
    start_time=rand(Uniform(0. - tδ,0.)) #select random time at or before pulse for cycle to begin
    label=check_cell_label(start_time, refractory, s, pulse) #fraction of nucleus labelled in first cycle
    set_prop!(g,1, :label, label*2.)
    set_prop!(g,2, :label, label) #first cycle terminates in mitotic event, children of which have half of the label from first cycle
    start_time+tδ>=end_time ? (mts=[]) : (mts=[2]) #mitoses to sim
    return GraphLineage(start_time,g),mts
end


function check_cell_label(t, refractory, s, pulse_length, curr_label=0.)
    lbl=curr_label
    if (t+refractory <= pulse_length && t+refractory+s >= 0) #some s-phase overlaps with the pulse in this case, track if detectable
        s_overlap=min(pulse_length,t+refractory+s)-max(0,t)
        nuc_label_frac=s_overlap/s
        lbl=min(nuc_label_frac+curr_label,1.)
    end

    if lbl >= DETECTION_THRESHOLD
        return lbl
    else
        return 0.
    end
end

function time_to_m(l,m)
    g=l.graph
    return l.start_time + sum([get_prop(g,e.src,e.dst,:weight) for e in a_star(g,1,m)])
end