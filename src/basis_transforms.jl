# # * General case

@doc raw"""
    basis_transform(A, B)

Compute the (possibly dense) basis transform matrix to go from `B` to
`A`, via [`basis_overlap`](@ref) for orthogonal basis `A`:

```math
B^HA
```
"""
basis_transform(::Ortho,         A, B) = basis_overlap(A, B)

@doc raw"""
    basis_transform(A, B)

Compute the (possibly dense) basis transform matrix to go from `B` to
`A`, via [`basis_overlap`](@ref) for non-orthogonal basis `A`:

```math
(B^HB)^{-1}B^HA
```
"""
basis_transform(::NonOrthogonal, A, B) = (A'A) \ basis_overlap(A, B)

"""
    basis_transform(A, B)

Compute the (possibly dense) basis transform matrix to go from `B` to
`A`, via [`basis_overlap`](@ref), dispatching on the
[`orthogonality`](@ref) trait of `A`.
"""
basis_transform(A, B) = basis_transform(orthogonality(A), A, B)

# * To finite-differences

"""
    basis_overlap(A::FiniteDifferencesOrRestricted, B)

Transforming to finite-differences is simple: just evaluate the basis
functions of the source basis on the grid points, possibly weighting
them.
"""
function basis_overlap(A::FiniteDifferencesOrRestricted, B)
    xa = locs(A)
    wa = weights(A)

    B[xa,:] ./ wa
end

function basis_transform(A::FiniteDifferencesOrRestricted, B, c::AbstractVector)
    xa = locs(A)
    wa = weights(A)

    (B[xa,:]*c) ./ wa
end

# * To FEDVR

"""
    basis_overlap(A::FEDVROrRestricted, B)

Transforming to FEDVR is much the same as for finite-differences;
evaluate at the quadrature nodes and scale by the weights. However,
this has to be done per-element, so we feed it through the
interpolation routines already setup for this.
"""
function basis_overlap(A::FEDVROrRestricted, B)
    T = promote_type(eltype(A), eltype(B))
    S = Matrix{T}(undef, size(A,2), size(B,2))

    xa = axes(A, 1)
    for j ∈ axes(B,2)
        # This is ugly; it is necessary since B[x,j] checks that x ⊆
        # axes(B, 1), whereas B[[x],j:j] does not, and for the purposes of the
        # transform, we wish to evaluate the basis functions as vanishing
        # outside the axis, which the latter does.
        f = x -> first(B[[x],j:j])
        S[:,j] .= A \ f.(xa)
    end

    S
end

function basis_transform(A::FEDVROrRestricted, B, c::AbstractVector)
    xa = axes(A, 1)
    # This is ugly; it is necessary since B[x,:] checks that x ⊆
    # axes(B, 1), whereas B[[x],:] does not, and for the purposes of the
    # transform, we wish to evaluate the basis functions as vanishing
    # outside the axis, which the latter does.
    f = x -> dot(B[[x],:],c)
    A \ f.(xa)
end

# * To B-splines
# ** From polynomial bases

function _basis_transform_2_bspline_common(A::BSplineOrRestricted,
                                          B::Union{BSplineOrRestricted,FEDVROrRestricted})
    k = order(A)
    k′ = max(k, maximum(order(B)))
    N = num_quadrature_points(k, k′-k+2)

    x,w = lgwt(knotset(A), N)
    W = Diagonal(w)

    x, A[x,:]'W
end

"""
    basis_overlap(A::BSplineOrRestricted, B::Union{BSplineOrRestricted,FEDVROrRestricted})

From a piecewise polynomial basis, i.e. B-splines or FEDVR, we can
use Gaussian quadrature to compute the overlap matrix elements
exactly. It is possible that we could derive better results, if we
used the fact that FEDVR basis functions are orthogonal in the sense
of the Gauss–Lobatto quadrature, but we leave that for later.
"""
function basis_overlap(A::BSplineOrRestricted, B::Union{BSplineOrRestricted,FEDVROrRestricted})
    T = promote_type(eltype(A), eltype(B))
    S = Matrix{T}(undef, size(A,2), size(B,2))

    x,χA = _basis_transform_2_bspline_common(A, B)

    for j ∈ axes(B,2)
        S[:,j] .= χA*B[x,j]
    end

    S
end

function basis_transform(A::BSplineOrRestricted, B::Union{BSplineOrRestricted,FEDVROrRestricted},
                         c::AbstractVector)
    T = promote_type(eltype(A), eltype(B), eltype(c))
    d = zeros(T, size(A, 2))

    x,χA = _basis_transform_2_bspline_common(A, B)

    for j ∈ axes(B,2)
        d .+= c[j] .* (χA*B[x,j])
    end

    (A'A) \ d
end

# ** From finite-differences (or any other basis)

"""
    basis_transform(A::BSplineOrRestricted, B::AbstractFiniteDifferences)

This works via Vandermonde interpolation, which is potentially
unstable. In this case, we do not provide [`basis_overlap`](@ref).
"""
function basis_transform(A::BSplineOrRestricted, B::AbstractFiniteDifferences)
    T = promote_type(eltype(A), eltype(B))
    S = Matrix{T}(undef, size(A,2), size(B,2))

    x = locs(B)
    χA = A[x,:]
    for j ∈ axes(B,2)
        S[:,j] .= χA \ Vector(B[x,j])
    end

    S
end

function basis_transform(A::BSplineOrRestricted, B::AbstractFiniteDifferences, c::AbstractVector)
    T = promote_type(eltype(A), eltype(B))
    N = size(B,2)

    x = locs(B)
    χA = A[x,:]
    χB = B[x,:]

    χA \ (χB*c)
end

# * Exports

export basis_transform, basis_overlap
