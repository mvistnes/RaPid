using Test
import RaPidSCOPF as SCOPF

@testset "RaPidSCOPF" begin
    include("test_dc_opf.jl")
    
    include("test_dc_power_flow.jl")
    
    include("test_imml.jl")
    
    include("test_ieee_rts_benders.jl")
end
