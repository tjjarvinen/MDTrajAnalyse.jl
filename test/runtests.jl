using Test
using MDTrajAnalyse


@testset "Cell" begin
    c = CubicCell(1.1)
    a = parsecell(1.1)
    aa = parsecell(1.1,1.1,1.1)
    naa = parsecell(1.2,1.1,1.1)
    @test c == a
    @test c == aa
    @test typeof(aa) != typeof(naa)

end


@testset "Trajectorys" begin
    t = Trajectory(rand(3,4,5))
    tt = PeriodicConstCellTrajectory(t.xyz, CubicCell(1.0))
    @test view(t,1) == view(tt,1)
    distances(t, 1,2)
    compute_rdf(tt, 1:1, 2:3)
end
