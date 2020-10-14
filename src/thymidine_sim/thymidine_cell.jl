const DETECTION_THRESHOLD=.01

mutable struct SPhaseCell <: AbstractCell 
    id::UUID
    parent::UUID

    mitotic::Bool

    events::Vector{Float64}
    label::Vector{Float64}
end


#[refractory][S][Tc---M]
function thymidine_cycle!(c::SPhaseCell,t::Float64, end_time::Float64, pulse_length::Float64, fate_dist::Categorical, refractory::Float64, Tcμ::Float64, Tcσ::Float64, s_frac::Float64)
    tδ,s=refractory_cycle_model(refractory, Tcμ, Tcσ, s_frac)

    push!(c.events,t+refractory+s) #s-phase exit is an event

    if (t+refractory <= pulse_length && t+refractory+s >= 0) #some s-phase overlaps with the pulse in this case, track if detectable
        s_overlap=min(pulse_length,t+refractory+s)-max(0,t)
        nuc_label_frac=s_overlap/s
        if nuc_label_frac+c.label[end] >= DETECTION_THRESHOLD
            push!(c.label,min(nuc_label_frac+c.label[end],1.))
        else
            push!(c.label,c.label[end])
        end
    else
        push!(c.label, c.label[end])
    end

    if t+tδ <= end_time #if the mitosis is to occur before the simulation end time, each daughter cell should inherit half the label fraction
        push!(c.events,t+tδ) #mitosis is an event
        c.label[end]/2 > DETECTION_THRESHOLD ? push!(c.label, c.label[end]/2) : push!(c.label, 0.)#zero out label if it goes below detection threshold
    end

    fate=rand(fate_dist) #get the fate of the daughter cells
    if fate==1 #PP mitosis
        d_mitotic=true
    elseif fate==2 #PD mitosis
        d_mitotic=false
    elseif fate==3 #DD mitosis
        d_mitotic=false
        c.mitotic=false
    end

    d=SPhaseCell(uuid4(),c.id, d_mitotic, [c.events[end]], [c.label[end]])
    return tδ,d
end

function refractory_cycle_model(refractory, Tcμ, Tcσ, s_frac)
    Tc=max(0,rand(Normal(Tcμ,Tcσ)))
    tδ=refractory + Tc
    s=tδ*s_frac
    return tδ, s
end