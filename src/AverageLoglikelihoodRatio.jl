module AverageLoglikelihoodRatio

using Random
using Statistics

export allr_column, allr_matrix, compare_pfms, ALLRResult

# ──────────────────────────────────────────────────────────────────────────── #
#  Types
# ──────────────────────────────────────────────────────────────────────────── #

"""
    ALLRResult

Result of comparing two PFMs via the Average Log-Likelihood Ratio.

# Fields
- `score::Float64`  – best ALLR score across all valid offsets
- `offset::Int`     – offset (0-indexed) of the shorter matrix on the longer one at the best alignment
- `pvalue::Float64` – p-value from the permutation test (NaN when not computed)
"""
struct ALLRResult
    score::Float64
    offset::Int
    pvalue::Float64
end

function Base.show(io::IO, r::ALLRResult)
    pstr = isnan(r.pvalue) ? "not computed" : string(round(r.pvalue; sigdigits=4))
    print(io, "ALLRResult(score=$(round(r.score; sigdigits=6)), offset=$(r.offset), pvalue=$pstr)")
end

# ──────────────────────────────────────────────────────────────────────────── #
#  Small epsilon to avoid log(0)
# ──────────────────────────────────────────────────────────────────────────── #

const PSEUDO_EPS = 1e-10

# ──────────────────────────────────────────────────────────────────────────── #
#  Column-level ALLR
# ──────────────────────────────────────────────────────────────────────────── #

"""
    allr_column(c1, c2; background=fill(0.25, 4))

Compute the Average Log-Likelihood Ratio between two column frequency
vectors `c1` and `c2` (each length-4, summing to ≈1) relative to
`background`.

```math
ALLR(C_1,C_2) = \\frac{
    \\sum_i f_{i,1} \\ln\\!\\left(\\frac{f_{i,2}}{b_i}\\right)
  + \\sum_i f_{i,2} \\ln\\!\\left(\\frac{f_{i,1}}{b_i}\\right)
}{2}
```
"""
function allr_column(c1::AbstractVector{<:Real},
                     c2::AbstractVector{<:Real};
                     background::AbstractVector{<:Real} = fill(0.25, 4))
    length(c1) == length(c2) == length(background) == 4 ||
        throw(ArgumentError("All vectors must have length 4 (A, C, G, T)."))

    s = 0.0
    @inbounds for i in 1:4
        f1 = max(c1[i], PSEUDO_EPS)
        f2 = max(c2[i], PSEUDO_EPS)
        bi = max(background[i], PSEUDO_EPS)
        s += f1 * log(f2 / bi) + f2 * log(f1 / bi)
    end
    return s / 2.0
end

# ──────────────────────────────────────────────────────────────────────────── #
#  Matrix-level ALLR  (aligned – same width)
# ──────────────────────────────────────────────────────────────────────────── #

"""
    allr_matrix(pfm1, pfm2; background=fill(0.25, 4))

Compute the average ALLR across aligned columns of two PFMs of **equal
width**.  Each PFM is a `4 × W` matrix whose columns are frequency vectors.
"""
function allr_matrix(pfm1::AbstractMatrix{<:Real},
                     pfm2::AbstractMatrix{<:Real};
                     background::AbstractVector{<:Real} = fill(0.25, 4))
    size(pfm1, 1) == size(pfm2, 1) == 4 ||
        throw(ArgumentError("PFMs must have 4 rows (A, C, G, T)."))
    size(pfm1, 2) == size(pfm2, 2) ||
        throw(ArgumentError("PFMs must have the same number of columns for allr_matrix. " *
                            "Use compare_pfms for matrices of different widths."))

    W = size(pfm1, 2)
    s = 0.0
    @inbounds for j in 1:W
        s += allr_column(view(pfm1, :, j), view(pfm2, :, j); background=background)
    end
    return s / W
end

# ──────────────────────────────────────────────────────────────────────────── #
#  Sliding-window best ALLR  (different widths allowed)
# ──────────────────────────────────────────────────────────────────────────── #

"""
    _best_allr(pfm_short, pfm_long; background)

Slide `pfm_short` (width `w`) over `pfm_long` (width `W ≥ w`) and return
`(best_score, best_offset)` where offset is 0-indexed.
"""
function _best_allr(pfm_short::AbstractMatrix{<:Real},
                    pfm_long::AbstractMatrix{<:Real};
                    background::AbstractVector{<:Real} = fill(0.25, 4))
    w = size(pfm_short, 2)  # shorter width
    W = size(pfm_long, 2)   # longer width

    best_score  = -Inf
    best_offset = 0

    for offset in 0:(W - w)
        sub_long = view(pfm_long, :, (offset + 1):(offset + w))
        sc = allr_matrix(pfm_short, sub_long; background=background)
        if sc > best_score
            best_score  = sc
            best_offset = offset
        end
    end
    return best_score, best_offset
end

# ──────────────────────────────────────────────────────────────────────────── #
#  Random PFM generation (for permutation null)
# ──────────────────────────────────────────────────────────────────────────── #

"""
    _random_pfm(width, background, rng)

Generate a random PFM of given `width` by sampling columns from a Dirichlet
distribution parameterised so that the expected frequencies equal `background`.
Uses the Gamma-variable trick: draw `g_i ∼ Gamma(α_i, 1)` and normalise.
A concentration parameter `α₀ = 4` is used (relatively flat).
"""
function _random_pfm(width::Int, background::AbstractVector{<:Real}, rng::AbstractRNG)
    α0 = 4.0                                     # concentration
    pfm = Matrix{Float64}(undef, 4, width)
    @inbounds for j in 1:width
        s = 0.0
        for i in 1:4
            # Gamma(α, 1) via the Marsaglia–Tsang method is hidden inside
            # Julia's randexp / randn; we use a simple shape > 1 approach.
            α_i = α0 * background[i]
            g = _rand_gamma(α_i, rng)
            pfm[i, j] = g
            s += g
        end
        # normalise to frequencies
        for i in 1:4
            pfm[i, j] /= s
        end
    end
    return pfm
end

"""
Simple Gamma(α, 1) sampler for α > 0 using the rejection method
(Ahrens–Dieter for α < 1, Marsaglia–Tsang for α ≥ 1).
"""
function _rand_gamma(α::Float64, rng::AbstractRNG)
    if α < 1.0
        # Boost: Gamma(α) = Gamma(α+1) * U^(1/α)
        return _rand_gamma(α + 1.0, rng) * rand(rng)^(1.0 / α)
    end
    # Marsaglia–Tsang method for α ≥ 1
    d = α - 1.0 / 3.0
    c = 1.0 / sqrt(9.0 * d)
    while true
        local x::Float64
        local v::Float64
        while true
            x = randn(rng)
            v = (1.0 + c * x)
            v > 0.0 && break
        end
        v = v^3
        u = rand(rng)
        if u < 1.0 - 0.0331 * (x^2)^2
            return d * v
        end
        if log(u) < 0.5 * x^2 + d * (1.0 - v + log(v))
            return d * v
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────── #
#  P-value by permutation
# ──────────────────────────────────────────────────────────────────────────── #

"""
    _permutation_pvalue(observed, pfm_input, width_target, n_perm, background, rng)

Estimate a p-value by generating `n_perm` random PFMs (of width
`width_target`), computing the best sliding ALLR of `pfm_input` against each,
and counting how often the random score ≥ `observed`.

Uses the standard pseudo-count formula:

```math
P = \\frac{\\#\\{S_{\\text{rand}} \\ge S_{\\text{obs}}\\} + 1}{N + 1}
```
"""
function _permutation_pvalue(observed::Float64,
                             pfm_input::AbstractMatrix{<:Real},
                             width_target::Int,
                             n_perm::Int,
                             background::AbstractVector{<:Real},
                             rng::AbstractRNG)
    count_ge = 0
    w_input = size(pfm_input, 2)
    for _ in 1:n_perm
        rand_pfm = _random_pfm(width_target, background, rng)
        # slide the shorter one over the longer one, same logic as compare_pfms
        if w_input <= width_target
            sc, _ = _best_allr(pfm_input, rand_pfm; background=background)
        else
            sc, _ = _best_allr(rand_pfm, pfm_input; background=background)
        end
        if sc >= observed
            count_ge += 1
        end
    end
    return (count_ge + 1) / (n_perm + 1)
end

# ──────────────────────────────────────────────────────────────────────────── #
#  Main entry point
# ──────────────────────────────────────────────────────────────────────────── #

"""
    compare_pfms(input, target;
                 background = fill(0.25, 4),
                 n_perm     = 1000,
                 compute_pvalue = true,
                 rng        = Random.default_rng())

Compare two Position Frequency Matrices using the **Average Log-Likelihood
Ratio** (ALLR, Ting Wang).

# Arguments
- `input`  : `4 × W₁` PFM (rows = A, C, G, T; columns = positions).
- `target` : `4 × W₂` PFM.

# Keyword arguments
- `background`     : length-4 vector of background nucleotide frequencies
                     (default: uniform `[0.25, 0.25, 0.25, 0.25]`).
- `n_perm`         : number of random PFMs used for the permutation test
                     (default `1000`).
- `compute_pvalue` : set to `false` to skip the permutation test and return
                     only the raw ALLR score (much faster).
- `rng`            : random number generator for reproducibility.

# Returns
An [`ALLRResult`](@ref) with fields `score`, `offset`, and `pvalue`.

# Example
```julia
using AverageLoglikelihoodRatio

# two toy PFMs (4 × W)
pfm1 = [0.9 0.1 0.1;
         0.05 0.8 0.1;
         0.025 0.05 0.7;
         0.025 0.05 0.1]

pfm2 = [0.85 0.1 0.1 0.25;
         0.05 0.75 0.1 0.25;
         0.05 0.05 0.7 0.25;
         0.05 0.1 0.1 0.25]

result = compare_pfms(pfm1, pfm2)
println(result.score)   # best ALLR
println(result.offset)  # where the short one aligns on the long one
println(result.pvalue)  # permutation p-value
```
"""
function compare_pfms(input::AbstractMatrix{<:Real},
                      target::AbstractMatrix{<:Real};
                      background::AbstractVector{<:Real} = fill(0.25, 4),
                      n_perm::Int = 1000,
                      compute_pvalue::Bool = true,
                      rng::AbstractRNG = Random.default_rng())
    size(input, 1) == 4 || throw(ArgumentError("`input` must have 4 rows (A, C, G, T)."))
    size(target, 1) == 4 || throw(ArgumentError("`target` must have 4 rows (A, C, G, T)."))
    length(background) == 4 || throw(ArgumentError("`background` must have length 4."))

    w_in = size(input, 2)
    w_tg = size(target, 2)

    # Slide the shorter PFM over the longer one
    if w_in <= w_tg
        best_score, best_offset = _best_allr(input, target; background=background)
    else
        best_score, best_offset = _best_allr(target, input; background=background)
    end

    # Permutation p-value
    pval = NaN
    if compute_pvalue
        pval = _permutation_pvalue(best_score, input, w_tg, n_perm, background, rng)
    end

    return ALLRResult(best_score, best_offset, pval)
end

end
