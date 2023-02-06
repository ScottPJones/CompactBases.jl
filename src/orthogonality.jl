"""
    BasisOrthogonality

Base type for the orthogonality trait of a basis.
"""
abstract type BasisOrthogonality end

"""
    Ortho <: BasisOrthogonality

Intermediate type for the trait for all orthogonal bases.
"""
abstract type Ortho <: BasisOrthogonality end

"""
    Orthogonal <: Ortho

Concrete type for the trait of orthogonal bases, i.e. with diagonal metric.
"""
struct Orthogonal <: Ortho end

"""
    Orthonormal <: Ortho

Concrete type for the trait of orthonormal bases, i.e. with metric `I`.
"""
struct Orthonormal <: Ortho end

"""
    NonOrthogonal <: BasisOrthogonality

Concrete type for the trait of non-orthogonal bases, i.e. with non-diagonal metric.
"""
struct NonOrthogonal <: BasisOrthogonality end

"""
    orthogonality(B)

Returns the [`BasisOrthogonality`](@ref) trait of the basis `B`.

# Examples

```julia
julia> orthogonality(BSpline(LinearKnotSet(7, 0, 1, 11)))
CompactBases.NonOrthogonal()

julia> orthogonality(FEDVR(range(0, stop=1, length=11), 7))
CompactBases.Orthonormal()

julia> orthogonality(StaggeredFiniteDifferences(0.1, 0.3, 0.1, 10.0))
CompactBases.Orthonormal()

julia> orthogonality(FiniteDifferences(range(0, stop=1, length=11)))
CompactBases.Orthogonal()
```
"""
function orthogonality end

export orthogonality
