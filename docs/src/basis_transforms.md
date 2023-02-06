# Transforms between different bases

Sometimes, we have a function ``f`` expanded over the basis functions
of one basis ``A``,
```math
f_A = A\vec{c}
```
and wish to re-express as an expansion over the
basis functions of another basis ``B``,
```math
f_B = B\vec{d}.
```
We thus want to solve ``f_A = f_B`` for ``\vec{d}``; this is very
easy, and follows the same reasoning as in [Solving
equations](@ref). We begin by projecting into the space ``B`` by using
the projector ``BB^H``:
```math
BB^HA\vec{c} =
BB^HB\vec{d}.
```
This has to hold for any ``B``, as long as the basis functions of
``B`` are linearly independent, and is then equivalent to
```math
B^HA\vec{c} =
B^HB\vec{d},
```
the solution of which is
```math
\vec{d} =
(B^HB)^{-1}
B^HA
\vec{c}
\defd
S_{BB}^{-1}
S_{BA}
\vec{c}.
```

[`basis_overlap`](@ref) computes ``S_{BA}=B^HA`` and
[`basis_transform`](@ref) the full transformation matrix
``(B^HB)^{-1}B^HA``; sometimes it may be desirable to perform the
transformation manually in two steps, since the combined matrix may be
dense whereas the constituents may be relatively sparse. For
orthogonal bases, ``B^HB`` is naturally diagonal.

There is no requirement that the two bases span the same intervals,
i.e. `axes(A,1)` and `axes(B,1)` need not be fully overlapping. Of
course, if they are fully disjoint, the transform matrix will be
zero. Additionally, some functions may not be well represented in a
particular basis, e.g. the step function is not well approximated by
B-splines (of order ``>1``), but poses no problem for
finite-differences, except at the discontinuity. To judge if a
transform is faithful thus amount to comparing the transformed
expansion coefficients with those obtained by expanding the function
in the target basis directly, instead of comparing reconstruction
errors.

## Examples

We will now illustrate transformation between B-splines, FEDVR, and
various finite-differences:
```julia
julia> A = BSpline(LinearKnotSet(7, 0, 2, 11))
BSpline{Float64} basis with typename(LinearKnotSet)(Float64) of order k = 7 on 0.0..2.0 (11 intervals)

julia> B = FEDVR(range(0.0, 2.5, 5), 6)[:,2:end-1]
FEDVR{Float64} basis with 4 elements on 0.0..2.5, restricted to elements 1:4, basis functions 2..20 ⊂ 1..21

julia> C = StaggeredFiniteDifferences(0.03, 0.3, 0.01, 3.0)
Staggered finite differences basis {Float64} on 0.0..3.066973760326056 with 90 points with spacing varying from 0.030040496962651875 to 0.037956283408020486

julia> D = FiniteDifferences(range(0.25,1.75,length=31))
Finite differences basis {Float64} on 0.2..1.8 with 31 points spaced by Δx = 0.05

julia> E = BSpline(LinearKnotSet(9, -1, 3, 15))[:,2:end]
BSpline{Float64} basis with typename(LinearKnotSet)(Float64) of order k = 9 on -1.0..3.0 (15 intervals), restricted to basis functions 2..23 ⊂ 1..23
```

We use these bases to expand two test functions:
```math
f_1(x) = \exp\left[
-\frac{(x-1)^2}{2\times0.2^2}
\right]
```
and
```math
f_2(x) = \theta(x-0.5) - \theta(x-1.5),
```
where ``\theta(x)`` is the [Heaviside step
function](https://en.wikipedia.org/wiki/Heaviside_step_function).

```julia
julia> gauss(x) = exp(-(x-1.0).^2/(2*0.2^2))
gauss (generic function with 1 method)

julia> u(x) = (0.5 < x < 1.5)*one(x)
u (generic function with 1 method)
```

If for example, we first expand ``f_1`` in B-splines (``A``), and then wish to
project onto the staggered finite-differences (``C``), we do the
following:

```julia
julia> c = A \ gauss.(axes(A, 1))
17-element Vector{Float64}:
 -0.0003757237742430358
  0.001655433792947186
 -0.003950561899347059
  0.00714774775995981
 -0.011051108016208545
  0.0173792655408227
  0.04090160980771093
  0.6336420285181446
  1.3884690430556483
  0.6336420285181449
  0.04090160980771065
  0.01737926554082302
 -0.011051108016208311
  0.007147747759958865
 -0.00395056189934631
  0.0016554337929467365
 -0.00037572377424269145

julia> S_CA = basis_transform(C, A)
90×17 SparseArrays.SparseMatrixCSC{Float64, Int64} with 1530 stored entries:
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⡇
⠛⠛⠛⠛⠛⠛⠛⠛⠃

julia> S_CA*c
90-element Vector{Float64}:
  3.805531055898412e-5
  1.6890846360940058e-5
 -4.5403814886370545e-5
 -4.121425310424971e-5
  1.2245471085512564e-5
  6.482348697494439e-5
  8.896755347867799e-5
  9.703404969825733e-5
  0.0001290654548866452
  0.00023543152355538763
  0.0004663267687150117
  ⋮
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0

julia> basis_transform(C, A, c) # Or transform without forming matrix
90-element Vector{Float64}:
  3.80553105589841e-5
  1.6890846360940065e-5
 -4.540381488637052e-5
 -4.121425310424972e-5
  1.2245471085512601e-5
  6.482348697494439e-5
  8.896755347867799e-5
  9.703404969825729e-5
  0.00012906545488664526
  0.0002354315235553877
  0.00046632676871501176
  ⋮
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
```

In the plots below, we show each basis (labelled ``A-E``), and their
respective reconstructions of ``f_1`` and ``f_2`` as solid red lines,
with the true functions shown as dashed black lines. For each
basis–function combination, shown are also the expansion coefficients
as red points joined by solid lines, as well as the expansion
coefficients transformed from the other bases, as labelled in the
legends. As an example, we see that first expanding in
finite-differences ``D`` and then transforming to B-splines ``A`` is
numerically unstable, leading to very large coefficients at the edges
of ``D``. Going _to_ finite-differences is always stable, even though
the results may not be accurate.

![Transforms between bases](../figures/basis_transforms.svg)

## Reference

```@docs
basis_overlap
basis_transform
```
