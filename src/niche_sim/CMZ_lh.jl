const MAXVAL=prevfloat(Inf)

function CMZ_mc_llh(popdist::Distribution, voldist::Distribution, volconst::Float64, phase_ends::Vector{Float64}, pparams::Vector{Float64}, mc_its::Int64, T::Vector{<:AbstractFloat}, obs::Vector{<:Tuple{<:AbstractVector{<:Float64},<:AbstractVector{<:Float64}}})
    times=length(T)
    pends=floor.(phase_ends)
    events=Int64.(sort(unique(vcat(T,pends))))
    results=zeros(mc_its,length(events),2)

    nt=Threads.nthreads()
    Threads.@threads for t in 1:nt
        t==nt ? (idxs=1+(t-1)*Int(floor(mc_its/nt)):mc_its) : (idxs=(1+(t-1)*Int(floor(mc_its/nt))):t*Int(floor(mc_its/nt)))
     
        phase=1
        cycle_time,exit_rate=pparams[1+((phase-1)*2):2+((phase-1)*2)]

        results[idxs,1,1].=rand.(popdist)
        results[idxs,1,2].=rand.(voldist)

        last_ph_ch=1;
        next_event=2;
        while next_event<=length(events)
            n=events[next_event]-events[last_ph_ch]

            pop_factor=max(eps(),(2^(24/cycle_time)-exit_rate))
            pop_factor==1. && (pop_factor=pop_factor+eps())
            vol_factor=volconst*exit_rate*(1-pop_factor^(n-1))/(1-pop_factor)

            lastpops=view(results,idxs,last_ph_ch,1)
            lastvols=view(results,idxs,last_ph_ch,2)

            results[idxs,next_event,1].=max.(min.(lastpops.*pop_factor^n,MAXVAL),eps())
            results[idxs,next_event,2].=min.(lastvols.+(lastpops.*vol_factor),MAXVAL)

            if events[next_event] in pends #after updating pop and vol to t, update phase parameters for next event, make t₀ vol₀ and pop₀ the phase change vals
                phase+=1; last_ph_ch=next_event
                phase <= length(pparams)/2 && ((cycle_time,exit_rate)=pparams[1+((phase-1)*2):2+((phase-1)*2)])
            end #advance phase if necessary

            next_event+=1
        end
    end

    tidxs=[findfirst(t->t==time,events) for time in T]

    pop_lns=Vector{LogNormal}(undef,times)
    vol_lns=Vector{LogNormal}(undef,times)
    pop_lns[1]=popdist
    vol_lns[1]=voldist

    try
        Threads.@threads for t in 2:times
                pop_lns[t]=fit(LogNormal,results[:,tidxs[t],1])
                vol_lns[t]=fit(LogNormal,results[:,tidxs[t],2])
        end
    catch
        return -Inf, zeros(0,0,0)
    end

    pop_lhs=Vector{Float64}(undef,times-1)
    vol_lhs=Vector{Float64}(undef,times-1)

    Threads.@threads for t in 1:times-1
        pop_lhs[t]=lps(logpdf(pop_lns[t+1],obs[t+1][1]))
        vol_lhs[t]=lps(logpdf(vol_lns[t+1],obs[t+1][2]))
    end
    
    log_lh=lps(lps(pop_lhs),lps(vol_lhs))

    disp_mat=zeros(times,3,2)
    Threads.@threads for t in 1:times
        disp_mat[t,:,1]=[quantile(pop_lns[t],.025),mean(pop_lns[t]),quantile(pop_lns[t],.975)]
        disp_mat[t,:,2]=[quantile(vol_lns[t],.025),mean(vol_lns[t]),quantile(vol_lns[t],.975)]
    end

    return log_lh, disp_mat
end