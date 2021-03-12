module CMZNicheSims
    using Distributions, GMC_NS, UnicodePlots, StatsPlots
    import ProgressMeter: @showprogress,move_cursor_up_while_clearing_lines
    import BioBackgroundModels:lps
    import StatsFuns:logsumexp
    import Serialization:serialize,deserialize
    import KernelDensity:kde,AbstractKDE

    const MAXVAL=prevfloat(Inf) #don't permit infinite CMZ populations or retinal volumes

    include("thymidine_sim/thymidine_cell.jl")
    include("thymidine_sim/thymidine_model.jl")
    include("thymidine_sim/thymidine_ensemble.jl")
    export Thymidine_Ensemble
    include("thymidine_sim/thymidine_sim.jl")
    include("thymidine_sim/thymidine_lh.jl")
    include("CMZ_sim/CMZ_lh.jl")
    include("CMZ_sim/CMZ_model.jl")
    include("CMZ_sim/CMZ_ensemble.jl")
    export CMZ_Ensemble
    include("slice_sim/lens_model.jl")
    export Lens_Model
    include("slice_sim/slice_lh.jl")
    include("slice_sim/slice_model.jl")
    include("slice_sim/slice_ensemble.jl")
    export Slice_Ensemble
    include("cycle_decay_sim/decay_lh.jl")
    include("cycle_decay_sim/decay_constructor.jl")
    include("cycle_decay_sim/decay_ensemble.jl")
    export Decay_Ensemble
    include("multislice_sim/multislice_model.jl")
    include("multislice_sim/multislice_ensemble.jl")
    export MultiSlice_Ensemble, MultiSlice_Decay_Ensemble
end # module
