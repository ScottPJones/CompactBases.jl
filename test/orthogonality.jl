@testset "Orthogonality" begin
    @test orthogonality(BSpline(LinearKnotSet(7, 0, 1, 11))) ==
        CompactBases.NonOrthogonal()

    @test orthogonality(FEDVR(range(0, stop=1, length=11), 7)) ==
        CompactBases.Orthonormal()

    SFD = StaggeredFiniteDifferences(0.1, 0.3, 0.1, 10.0)
    @test orthogonality(SFD) == CompactBases.Orthonormal()
    @test orthogonality(SFD[:,2:end-1]) == CompactBases.Orthonormal()

    @test orthogonality(StaggeredFiniteDifferences(11, 0.1)) ==
        CompactBases.Orthogonal()

    @test orthogonality(FiniteDifferences(range(0, stop=1, length=11))) ==
        CompactBases.Orthogonal()
end
