@testset "Basis transforms" begin
    A = BSpline(LinearKnotSet(7, 0, 2, 11))
    A′ = A[:,2:end-1]
    B = FEDVR(range(0.0, stop=2.5, length=5), 6)[:,2:end-1]
    C = StaggeredFiniteDifferences(0.03, 0.3, 0.01, 3.0)
    D = FiniteDifferences(range(0.25, stop=1.75, length=31))
    E = BSpline(LinearKnotSet(9, -1, 3, 15))[:,2:end]
    F = FiniteDifferences(range(1, stop=4, step=0.1))

    ex(B,f) = B \ f.(axes(B,1))
    bt(A,B,f) = basis_transform(A,B,ex(B,f))
    bt2(A,B,f) = basis_transform(A,B)*ex(B,f)

    gauss = x -> exp(-(x-1.0).^2/(2*0.2^2))
    # u = x -> (0.5 < x < 1.5)*one(x)

    for (A,B,f,rtol) in [(A,A,gauss,1e-13),
                         (A′,A,gauss,1e-3),
                         (B,A,gauss,5e-3),
                         (C,A,gauss,5e-3),
                         (D,A,gauss,5e-3),
                         (E,A,gauss,1e-1),
                         (F,A,gauss,5e-3),

                         (A,A′,gauss,1e-3),
                         (A′,A′,gauss,1e-13),
                         (B,A′,gauss,5e-3),
                         (C,A′,gauss,5e-3),
                         (D,A′,gauss,5e-3),
                         (E,A′,gauss,1e-1),
                         (F,A′,gauss,5e-3),

                         (A,B,gauss,1e-1),
                         (A′,B,gauss,1e-1),
                         (B,B,gauss,1e-13),
                         (C,B,gauss,1e-2),
                         (D,B,gauss,1e-2),
                         (E,B,gauss,1e-1),
                         (F,B,gauss,1e-2),

                         (A,C,gauss,1e-2),
                         (A′,C,gauss,5e-3),
                         (B,C,gauss,1e-2),
                         (C,C,gauss,1e-13),
                         (D,C,gauss,1e-2),
                         (F,C,gauss,1e-2)]
        c = ex(A, f)
        d = bt(A, B, f)
        d2 = bt2(A, B, f)

        (norm(c-d)/norm(c) ≥ rtol || norm(c-d2)/norm(c) ≥ rtol) &&
            @error "Basis transform accuracy test failed" A B rtol norm(c-d)/norm(c) norm(c-d)/norm(c) < rtol norm(c-d2)/norm(c) norm(c-d2)/norm(c) < rtol

        @test d ≈ c rtol=rtol
        @test d2 ≈ c rtol=rtol
    end
end
