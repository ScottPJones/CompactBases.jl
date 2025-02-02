using CompactBases
import CompactBases: locs,
    unrestricted_basis, restriction_extents

using IntervalSets
using QuasiArrays
import QuasiArrays: AbstractQuasiArray,  AbstractQuasiMatrix, MulQuasiArray
using ContinuumArrays
import ContinuumArrays: ℵ₁, Inclusion

using LinearAlgebra
using BandedMatrices
using BlockBandedMatrices
using SparseArrays

using LazyArrays
import LazyArrays: materialize, Dot
using FillArrays

using ArnoldiMethod
using Random

using Test

function vecdist(a::AbstractVector, b::AbstractVector,
                 ϵ = eps(eltype(a)))
    δ = √(sum(abs2, a-b))
    δ, δ/√(sum(abs2, a .+ ϵ))
end

include("derivative_accuracy_utils.jl")

@testset "CompactBases" begin
    include("restricted_bases.jl")
    include("fd/runtests.jl")
    include("fedvr/runtests.jl")
    include("bsplines/runtests.jl")
    include("interpolation.jl")
    include("inner_products.jl")
    include("densities.jl")
    include("linear_operators.jl")
    include("orthogonality.jl")
    include("basis_transforms.jl")
end
