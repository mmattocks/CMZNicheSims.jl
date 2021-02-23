struct Lens_Model
    factor::Float64
    power::Float64
    section::Float64
end

function circumferential_exit(lm,t1,t2,pops)
    t1c=lm.factor*(t1^lm.power)
    sliceno=t1c/lm.section

    t2c=lm.factor*(t2^lm.power)
    nc=t2c-t1c
    psg=nc/sliceno
    return (psg/lm.section).*pops
end

function Base.length(lm::Lens_Model)
    return 1
end