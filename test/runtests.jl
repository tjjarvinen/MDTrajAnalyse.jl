using Test
using MDTrajAnalyse


@testset "Cell" begin
    c = CubicCell(1.1)
    a = parsecell(1.1)
    aa = parsecell(1.1,1.1,1.1)
    @test c == a
    @test c == aa

end
