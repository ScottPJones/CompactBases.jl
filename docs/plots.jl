using PythonPlot
using Jagot.plotting
plot_style("ggplot")
using PythonCall
PathEffects = pyimport("matplotlib.patheffects")
using Statistics
using Random
using Colors

using LinearAlgebra
using ArnoldiMethod

using CompactBases
import CompactBases: applied
using IntervalSets

function mean_color(color::String)
    c = parse(Colorant, color)
    mean([c.r,c.g,c.b])
end

mean_color(color::Py) = mean_color(pyconvert(String, color))

lerp(a,b,t) = (1-t)*a + t*b

mean_position(x, ϕ) = ϕ'*Diagonal(x)*ϕ/(ϕ'ϕ)

function savedocfig(name,dir="figures")
    fig = gcf()
    filename = joinpath(@__DIR__, "src", dir, "$(name).svg")
    savefig(filename,
            transparent=true,
            facecolor=fig.get_facecolor())
    PythonPlot.close("all")
    if isfile(filename)
        println("Saved $(name) to $(filename)")
    else
        @warn "Saving $(name) to $(filename) failed"
    end
end

function logo()
    t = ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3)
    r = range(first(t), stop=last(t), length=301)
    R = BSpline(t)
    χ = R[r, :]

    cfigure("logo",figsize=(7,3)) do
        plot(r, χ, length(r) < 30 ? ".-" : "-", linewidth=4)
        margins(0.1, 0.1)
        axis("off")
    end
    savedocfig("logo", "assets")
end

function simple_example()
    a,b = 0.0,1.0 # Extents
    N = 13 # Number of nodes
    k = 5 # Order of FE-DVR/B-splines
    x = range(a, stop=b, length=1000)

    Δx = (b-a)/(N+1) # Grid spacing
    # Standard, uniform finite-differences
    fd = FiniteDifferences(N, Δx)
    # Staggered, uniform finite-differences
    sfd_uni = StaggeredFiniteDifferences(N, Δx)
    # Staggered, non-uniform finite-differences
    sfd_nonuni = StaggeredFiniteDifferences(0.01, 0.5, 0.1, b)

    # Finite-element boundaries
    tf = range(a, stop=b, length=N+2)
    # We can vary the polynomial order in each element
    forder = vcat(7, fill(4,length(tf)-2))
    # By indexing the second dimension, we can implement Dirichlet0
    # boundary conditions.
    fem = FEDVR(tf, forder)[:,2:end-1]

    tb = ExpKnotSet(k, -2.0, log10(b), N+1)
    splines = BSpline(tb)[:,2:end-1]

    bases = [fd, sfd_uni, sfd_nonuni, fem, splines]
    m = length(bases)

    cfigure("simple example", figsize=(6,9)) do
        for (i,basis) in enumerate(bases)
            csubplot(m,1,i,nox=i<m) do
                plot(x, basis[x,:])
            end
        end
    end
    savedocfig("simple_example")
end

function restricted_bases()
    B = FEDVR(range(0, stop=10.0, length=5), 7)
    B̃ = B[:, 3:17]

    x = range(0, stop=10.0, length=1000)
    χ = B[x, :]
    χ̃ = B̃[x, :]

    cfigure("restricted bases") do
        csubplot(211,nox=true) do
            plot(x, χ)
            ylabel("Original basis")
        end
        csubplot(212) do
            plot(x, χ̃)
            ylabel("Restricted basis")
        end
    end
    savedocfig("restricted_bases")
end

function densities()
    f = x -> sin(2π*x)
    g = x -> x*exp(-x)
    h = x -> f(x)*g(x)

    rmin = 0.01
    rmax = 10.0
    ρ = 0.05
    N = ceil(Int, rmax/ρ)
    k = 7
    Nn = 71

    ρmin=rmin
    ρmax=0.5
    α=0.01

    coeff_centers = Vector{Vector{Float64}}()
    density_coeffs = Vector{Vector{Float64}}()
    errors = Vector{Vector{Float64}}()
    ϕ = Vector{Vector{Float64}}()

    rr = 10.0 .^ range(-2, stop=log10(rmax), length=3000)

    for R in [FiniteDifferences(N, ρ),
              StaggeredFiniteDifferences(ρmin, ρmax, α, rmax),
              FEDVR(range(0, stop=rmax, length=Nn), k),
              BSpline(LinearKnotSet(k, 0, rmax, Nn))]
        r = axes(R,1)

        cf = R \ f.(r)
        cg = R \ g.(r)
        ch = R \ h.(r)

        ρ = Density(applied(*,R,cf), applied(*,R,cg))

        push!(coeff_centers, centers(R))
        push!(density_coeffs, ρ.ρ)
        push!(errors, ρ.ρ-ch)
        push!(ϕ, R[rr,:]*ρ.ρ)
    end

    cfigure("mutual density", figsize=(7,10)) do
        csubplot(411,nox=true) do
            semilogx(rr, f.(rr), linewidth=1.0)
            semilogx(rr, g.(rr), linewidth=1.0)
            semilogx(rr, h.(rr), linewidth=2.0)
            legend([L"f(x)", L"g(x)", L"h(x)"])
        end
        csubplot(412, nox=true) do
            for ϕ in ϕ
                loglog(rr, abs.(ϕ-h.(rr)))
            end
            yl = pyconvert(Tuple, ylim())
            ylim(max(yl[1],1e-15), yl[2])
            legend(["Finite-differences",
                    "Log–linear finite-differences",
                    "FE-DVR",
                    "B-splines"])
            ylabel("Reconstruction error")
        end
        csubplot(413, nox=true) do
            for (r,c) in zip(coeff_centers, density_coeffs)
                semilogx(r, c, ".-", linewidth=1.0)
            end
            legend(["Finite-differences",
                    "Log–linear finite-differences",
                    "FE-DVR",
                    "B-splines"])
            ylabel("Expansion coefficients")
        end
        csubplot(414) do
            for (r,e) in zip(coeff_centers, errors)
                loglog(r, abs.(e), ".-", linewidth=1.0)
            end
            legend(["Finite-differences",
                    "Log–linear finite-differences",
                    "FE-DVR",
                    "B-splines"])
            ylabel("Coefficient error")
            xlabel(L"r")
        end
    end

    savedocfig("mutual_densities")
end

function diagonal_operators()
    Z = 1.0
    rmax = 20*(2^1.05)/Z

    ρ = 0.1/Z
    N = ceil(Int, rmax/ρ)
    ufd = StaggeredFiniteDifferences(N, ρ)

    ρmin=0.1/Z
    ρmax=0.6
    α=0.002
    nufd = StaggeredFiniteDifferences(ρmin, ρmax, α, rmax)

    fedvr = FEDVR(range(0, stop=rmax, length=20),
                  vcat(10,fill(7,18)))[:,2:end-1]

    bsplines = BSpline(ExpKnotSet(5, -2.0, log10(rmax), 100))[:,2:end-1]

    R = bsplines

    r = axes(R,1)
    rr = clamp.(10.0 .^ range(min(-2,log10(r[1])), stop=log10(rmax), length=3000),
                0, rightendpoint(R))
    χ = R[rr,:]

    r = axes(R,1)
    V = r -> (0.314269680527354*r^2 + 0.20951312035157*r)*exp(-3*r/2);
    Vc = R \ V.(r);
    rc = R \ identity.(r);

    L = DiagonalOperator(applied(*, R, Vc));

    Vr = L*rc

    Vr_exact = rr .* V.(rr)

    cfigure("diagonal operators") do
        csubplot(211,nox=true) do
            semilogx(rr, χ*Vr, label=L"V(r)r")
            semilogx(rr, Vr_exact, "--", label=L"$V(r)r$ exact")
            semilogx(rr, V.(rr), label=L"V(r)")
            yl = ylim()
            semilogx(rr, rr, label=L"r")
            ylim(yl)
            legend()
        end
        csubplot(212) do
            for (R,label) in [
                (ufd, "Uniform finite-differences"),
                (nufd, "Non-uniform finite-differences"),
                (fedvr, "FE-DVR"),
                (bsplines, "B-splines")
            ]
                r = axes(R,1)
                Vc = R \ V.(r)
                rc = R \ identity.(r)

                L = DiagonalOperator(applied(*, R, Vc))

                Vr = L*rc
                χ = R[rr,:]
                loglog(rr, abs.(χ*Vr - Vr_exact), label=label)
            end
            xlabel(L"r")
            ylabel("Error")
            legend()
        end
    end

    savedocfig("diagonal_operators")
end

function basis_transforms()
    function compare_bases(x, bases, functions)
        m = length(bases)
        n = length(functions)

        cs = [B \ f.(axes(B,1))
              for B in bases, f in functions]
        es = extrema.(cs)
        ymins = [minimum(first.(es[i,:])) for i = 1:m]
        ymaxs = [maximum(last.(es[i,:])) for i = 1:m]

        cfigure("basis comparison", figsize=(10,16)) do
            for i = 1:m
                B = bases[i]
                a = axes(B,1)
                χ = B[x,:]
                xl = [mean_position(x, view(χ, :, j)) for j in 1:size(B,2)]
                csubplot(m,n+1,(i,1), nox=i<m) do
                    plot(x, χ)
                    i == m && xlabel(L"x")
                    reltext(0.1, 0.8, string('A'+i-1))
                end
                for j = 1:n
                    f = functions[j]
                    c = cs[i,j]
                    csubplot(2m,n+1,(2i-1,j+1), nox=true, noy=j<n) do
                        plot(x, χ*c)
                        plot(x, f.(x), "k--")
                        margins(0,0.1)
                        ylim(-0.2,1.25)
                        axes_labels_opposite(:y)
                        i == m && xlabel(L"x")
                    end
                    csubplot(2m,n+1,(2i,j+1), nox=i<m, noy=j<n) do
                        plot(xl, c, ".-")
                        for i′ = 1:m
                            i′ == i && continue
                            A = bases[i′]
                            d = basis_transform(B, A, cs[i′,j])
                            plot(xl, d, ".--", label=string('A'+i′-1))
                        end
                        xlim(x[1],x[end])
                        dy = (ymaxs[i]-ymins[i])
                        ylim(ymins[i]-0.1dy, ymaxs[i]+0.1dy)
                        legend(loc=2)
                        axes_labels_opposite(:y)
                        i == m && xlabel(L"x")
                    end
                end
            end
        end
    end

    A = BSpline(LinearKnotSet(7, 0, 2, 11))
    B = FEDVR(range(0.0, 2.5, 5), 6)[:,2:end-1]
    C = StaggeredFiniteDifferences(0.03, 0.3, 0.01, 3.0)
    D = FiniteDifferences(range(0.25,1.75,length=31))
    E = BSpline(LinearKnotSet(9, -1, 3, 15))[:,2:end]

    gauss = x -> exp(-(x-1.0).^2/(2*0.2^2))
    u = x -> (0.5 < x < 1.5)*one(x)

    x = range(-1.5, stop=3.5, length=1001)

    compare_bases(x, (A,B,C,D,E), (gauss,u))

    savedocfig("basis_transforms")
end

macro echo(expr)
    println(expr)
    :(@time $expr)
end

@info "Documentation plots"
fig_dir = joinpath(@__DIR__, "src", "figures")
mkpath(fig_dir)
@assert isdir(fig_dir)
@echo logo()
@echo simple_example()
@echo restricted_bases()
include("bspline_plots.jl")
include("fd_plots.jl")
@echo densities()
@echo diagonal_operators()
@echo basis_transforms()
