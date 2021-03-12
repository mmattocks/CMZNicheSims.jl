function decay_slice_mc_llh(popdist::Distribution, lm::Lens_Model, dparams::Vector{Float64}, mc_its::Int64, T::Vector{Float64}, obs::Vector{Vector{Float64}})
    times=length(T)
    events=Int64.(sort(T))
    results=zeros(mc_its,length(events))

    nt=Threads.nthreads()

    init_day=events[1]
    init_ct, decay_const, exit_rate=dparams
    cts=[init_ct*exp(decay_const*elapsed) for elapsed in collect(events[1]:1:events[end]).-init_day]

    Threads.@threads for t in 1:nt
        t==nt ? (idxs=1+(t-1)*Int(floor(mc_its/nt)):mc_its) : (idxs=(1+(t-1)*Int(floor(mc_its/nt))):t*Int(floor(mc_its/nt)))

        results[idxs,1].=rand.(popdist)

        day=init_day; next_event=2
        curr_pops=copy(results[idxs,1,1])
        while day <= events[end]
            days_elapsed = day - init_day
            cycle_time=cts[day-init_day+1]
            pop_factor=max(eps(),(2^(24/cycle_time)-exit_rate))
            circexit=circumferential_exit(lm,day-1,day,curr_pops)

            curr_pops=max.(min.(curr_pops.*pop_factor,MAXVAL).-circexit,eps())

            day in events[2:end] && (
                results[idxs,next_event].=curr_pops; 
                next_event+=1)
            day+=1
        end
    end

    tidxs=[findfirst(t->t==time,events) for time in T]

    pop_lns=Vector{LogNormal}(undef,times)
    pop_lns[1]=popdist
    Threads.@threads for t in 2:times
            pop_lns[t]=fit(LogNormal,results[:,tidxs[t]])
    end

    pop_lhs=Vector{Float64}(undef,times-1)

    Threads.@threads for t in 1:times-1
        pop_lhs[t]=lps(logpdf(pop_lns[t+1],obs[t+1]))
    end
 
    log_lh=lps(pop_lhs)

    disp_mat=zeros(times,3)
    Threads.@threads for t in 1:times
        disp_mat[t,:]=[quantile(pop_lns[t],.025),mean(pop_lns[t]),quantile(pop_lns[t],.975)]
    end

    return log_lh, disp_mat
end