const MAXVAL=prevfloat(Inf)

function CMZ_mc_llh(popdist, voldist, volconst, phase_ends, pparams, noise_dist, mc_its, T::Vector{<:AbstractFloat}, obs::Vector{<:Tuple{<:AbstractVector{<:Float64},<:AbstractVector{<:Float64}}})
    times=length(T)
    results=zeros(mc_its,times,2)
    pends=floor.(phase_ends)
    events=sort(unique(vcat(T,pends)))

    Threads.@threads for it in 1:mc_its
        phase=1
        cycle_time,exit_rate=pparams[1+((phase-1)*2):2+((phase-1)*2)]

        pop₀=rand(popdist); vol₀=rand(voldist)
        popₜ=pop₀;volₜ=vol₀
        results[it,1,:]=[popₜ,volₜ]
        
        t₀=Int64(events[1]); t=t₀; next_event=2; next_tidx=2
        while next_event<=length(events)
            popₜ==0 && (results[it,next_tidx:end,:].=[0. vol];break) 
            popₜ==MAXVAL && (results[it,next_tidx:end,:].=[MAXVAL volₜ];break)
            volₜ==MAXVAL && (results[it,next_tidx:end,:].=[popₜ MAXVAL];break)
          
            n=Int64(events[next_event])-t₀

            pop_factor=(2^(24/cycle_time)-exit_rate)
            pop_factor==1. && (pop_factor=pop_factor+eps())
            volₜ=vol₀+((volconst*pop₀*exit_rate*(1-pop_factor^(n-1)))/(1-pop_factor))
            popₜ=pop₀*pop_factor^n

            volₜ=ifelse(volₜ<MAXVAL, volₜ,MAXVAL)
            popₜ=ifelse(popₜ>0., ifelse(popₜ<MAXVAL,popₜ,MAXVAL),0.)

            if t in pends #after updating pop and vol to t, update phase parameters for next event, make t₀ vol₀ and pop₀ the phase change vals
                phase+=1
                phase <= length(pparams)/2 && ((cycle_time,exit_rate)=pparams[1+((phase-1)*2):2+((phase-1)*2)])
                phase>length(phase_ends) ? (phase_end=maximum(T)+1) : (phase_end=phase_ends[phase])
                t₀=t; vol₀=volₜ; pop₀=popₜ
            end #advance phase if necessary

            t in T && (results[it,next_tidx,:]=[popₜ,volₜ]; next_tidx+=1)#write pop,vol

            t=events[next_event];  next_event+=1
        end
    end

    pop_lns=Vector{LogNormal}(undef,times)
    vol_lns=Vector{LogNormal}(undef,times)
    try
        Threads.@threads for t in 1:times
                pop_lns[t]=fit(LogNormal,results[:,t,1])
                vol_lns[t]=fit(LogNormal,results[:,t,2])
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