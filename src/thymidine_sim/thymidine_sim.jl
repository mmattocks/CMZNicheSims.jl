DETECTION_THRESHOLD=.15
MIN_CT=4.

function sim_lin_pair!(plv, T, end_time, pulse, Tc, g1_frac, s_frac, sis_frac)
    tδ1, g1_1, s1 = cycle_model(Tc, g1_frac, s_frac)
    tδ2, g1_2, s2 = (1+rand([-1,1])*sis_frac).*(tδ1, g1_1, s1)

    mint=-min(tδ1,tδ2)
    mint==0. ? (birth_time=mint) : (birth_time=rand(Uniform(mint,0.)))

    cells1=[LabelCell(birth_time,tδ1,g1_1,s1,0.)]
    cells2=[LabelCell(birth_time,tδ2,g1_2,s2,0.)]

    for cellvec in [cells1,cells2]
        ci=1
        while ci<=length(cellvec)    
            cycle_sim!(plv, ci, cellvec, T, end_time, pulse, Tc, g1_frac,s_frac, sis_frac)
            ci+=1
        end
    end    
end


function sim_lineage!(plv, T, end_time, pulse, Tc, g1_frac, s_frac, sis_frac)
    tδ1, g1, s1 = cycle_model(Tc, g1_frac, s_frac)
    birth_time=rand(Uniform(-tδ1,0.))
    cells=[LabelCell(birth_time,tδ1,g1,s1,0.)]

    ci=1
    while ci<= length(cells)
        cycle_sim!(plv, ci, cells, T, end_time, pulse, Tc, g1_frac, s_frac, sis_frac)
        ci+=1
    end
end

function cycle_sim!(plv, ci::Int64, cells::Vector{LabelCell}, T::Vector{Float64}, end_time::Float64, pulse::Float64, Tc::LogNormal, g1_frac::Float64, s_frac::Float64, sis_frac::Float64)
    n_times=length(T);

    cell=cells[ci]
    t=cell.time; tδ′=cell.tδ′; g1′=cell.g1′; s′=cell.s′; l=cell.label

    first_cycle=true

    while t<end_time
        tδ,s=0.,0.
        first_cycle ? (tδ=tδ′;g1=g1′;s=s′;first_cycle=false) : ((tδ, g1, s)=cycle_model(Tc, g1_frac, s_frac))

        s_start=t+g1
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

        tδ2, g1_2, s2 = (1+rand([-1,1])*sis_frac).*(tδ, g1, s)
        tδ2 < MIN_CT && (tδ2=tδ; g1_2=g1; s2 = s;) #prevent sister shift from resulting in unrealistically fast proliferation

        cycle_mitosis < end_time && push!(cells,LabelCell(cycle_mitosis,tδ2,g1_2,s2,l))
        t+=tδ
    end
end

function cycle_model(Tc, g1_frac, s_frac)
    tδ=max(MIN_CT,rand(Tc)) #floor of 1 hr for any distribution of cycle times
    s=tδ * s_frac
    g1=g1_frac * (tδ - s)
    return tδ, g1, s
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
