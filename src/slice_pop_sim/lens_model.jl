struct Lens_Model
    factor::Float64
    power::Float64
    section::Float64
end

function circumferential_exit(lm,t1,t2,pop)
    t1c=lm.factor*(t1^lm.power)
    sliceno=t1c/14

    t2c=lm.factor*(t2^lm.power)
    nc=t2c-t1c
    psg=nc/sliceno
    return (psg/lm.section)*pop
end