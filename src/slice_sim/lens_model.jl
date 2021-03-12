struct Lens_Model
    factor::Float64
    power::Float64
    section::Float64
end

function circumferential_exit(lm,t1,t2,pops)
    t1c=lm.factor*(t1^lm.power) #lens circumference at time 1, according to lens model
    t2c=lm.factor*(t2^lm.power) #lens circumference at time 2
    nc=t2c-t1c #new circumference grown btw t1 and t2


    sliceno=t1c/lm.section #number of slices implied by lm's section size for t1 lens diameter
    psg=nc/sliceno #amount of circumferential growth each CMZ slice is responsible for (per slice growth)
    return (psg/lm.section).*pops #per slice growth divided by slice section thickness supplies the implied number of RPCs leaving CMZ to supply growing retina
end

function Base.length(lm::Lens_Model)
    return 1
end