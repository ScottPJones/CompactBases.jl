using Documenter, CompactBases, LinearAlgebra

isdefined(Main, :NOPLOTS) && NOPLOTS || include("plots.jl")

makedocs(;
    modules=[CompactBases],
    format=Documenter.HTML(assets = ["assets/latex.js"],
                           mathengine = Documenter.MathJax()),
    pages = [
        "Home" => "index.md",
        "Overview" => "overview.md",
        "Inner products & norms" => "inner_products.md",
        "Operators" => [
            "Diagonal operators" => "diagonal_operators.md",
        ],
        "Densities" => "densities.md",
        "Finite-differences" => [
            "Overview" => "fd_overview.md",
            "Non-uniform grids" => "fd_non_uniform.md",
        ],
        "B-splines" => [
            "Theory" => "bsplines/theory.md",
            "Usage" => [
                "Basis creation" => "bsplines/usage.md",
                "Knot sets" => "bsplines/knot_sets.md",
                "Splines" => [
                    "Spline creation & evaluation" => "bsplines/splines.md",
                    "Function approximation" => "bsplines/function_approximation.md",
                ],
                "Approximating operators" => "bsplines/operators.md",
                "Examples" => [
                    "Differentiating functions" => "bsplines/differentiation.md",
                    "Ordinary differential equations" => "bsplines/odes.md",
                    "Eigenproblems" => "bsplines/eigenproblems.md"
                ]
            ],
        ],
        "Basis transforms" => "basis_transforms.md"
    ],
    repo="https://github.com/JuliaApproximation/CompactBases.jl/blob/{commit}{path}#L{line}",
    sitename="CompactBases.jl",
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com>",
    doctest=false
)

deploydocs(;
    repo="github.com/JuliaApproximation/CompactBases.jl",
    push_preview = true,
)
