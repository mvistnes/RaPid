using Test
import .RaPidSCOPF as SCOPF

@testset "Test dc power flow" begin
    include("test_dc_power_flow.jl")
end

@testset "Test IMML" begin
    include("test_imml.jl")
end

@testset "Test Benders" begin
    include("test_ieee_rts_benders.jl")
end
