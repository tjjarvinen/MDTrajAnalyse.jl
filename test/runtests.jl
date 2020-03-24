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

    v = VariableCell([CubicCell(2.3), CubicCell(3.4), OrthorombicCell(2.3,2.4,2.5)] )
    @test length(v) == 3
    @test lastindex(v) == 3

    q = iterate(v)

    @test q[1] == v[1]

    @test length(volume(v)) == 3
end


@testset "Trajectorys" begin
    t = Trajectory(rand(3,4,5))
    tt = Trajectory(t.xyz; cell=CubicCell(1.0))
    @test view(t,1) == view(tt,1)
    distances(t, 1,2)
    compute_rdf(tt, 1:1, 2:3)
end
