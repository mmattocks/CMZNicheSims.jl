const MAXVAL=prevfloat(Inf)

function slice_mc_llh(popdist::Distribution, lm::Lens_Model, phase_ends::Vector{Float64}, pparams::Vector{Float64}, mc_its::Int64, T::Vector{Float64}, obs::Vector{Vector{Float64}})
    times=length(T)
    pends=floor.(phase_ends)
    events=Int64.(sort(unique(vcat(T,pends))))
    results=zeros(mc_its,length(events))

    nt=Threads.nthreads()
    Threads.@threads for t in 1:nt
        t==nt ? (idxs=1+(t-1)*Int(floor(mc_its/nt)):mc_its) : (idxs=(1+(t-1)*Int(floor(mc_its/nt))):t*Int(floor(mc_its/nt)))
     
        phase=1
        cycle_time,exit_rate=pparams[1+((phase-1)*2):2+((phase-1)*2)]

        results[idxs,1,1].=rand.(popdist)

        last_ph_ch=1;
        next_event=2;
        while next_event<=length(events)
            n=events[next_event]-events[last_ph_ch]

            pop_factor=max(eps(),(2^(24/cycle_time)-exit_rate))
            pop_factor==1. && (pop_factor=pop_factor+eps())

            lastpops=view(results,idxs,last_ph_ch)
            circ_exit=circumferential_exit(lm,events[last_ph_ch],events[next_event],lastpops)

            results[idxs,next_event].=max.(min.(lastpops.*pop_factor^n,MAXVAL).-circ_exit,eps())

            if events[next_event] in pends #after updating pop and vol to t, update phase parameters for next event, make t₀ vol₀ and pop₀ the phase change vals
                phase+=1; last_ph_ch=next_event
                phase <= length(pparams)/2 && ((cycle_time,exit_rate)=pparams[1+((phase-1)*2):2+((phase-1)*2)])
            end #advance phase if necessary

            next_event+=1
        end
    end

    tidxs=[findfirst(t->t==time,events) for time in T]

    pop_lns=Vector{LogNormal}(undef,times)
    pop_lns[1]=popdist
    try
        Threads.@threads for t in 2:times
                pop_lns[t]=fit(LogNormal,results[:,tidxs[t]])
        end
    catch
        return -Inf, zeros(0,0,0)
    end

    pop_lhs=Vector{Float64}(undef,times-1)

    Threads.@threads for t in 1:times-1
        pop_lhs[t]=lps(logpdf(pop_lns[t],obs[t+1]))
    end
    
    log_lh=lps(pop_lhs)

    disp_mat=zeros(times,3)
    Threads.@threads for t in 1:times
        disp_mat[t,:]=[quantile(pop_lns[t],.025),mean(pop_lns[t]),quantile(pop_lns[t],.975)]
    end

    return log_lh, disp_mat
end