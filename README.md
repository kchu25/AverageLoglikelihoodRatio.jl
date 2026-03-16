# AverageLoglikelihoodRatio.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/AverageLoglikelihoodRatio.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/AverageLoglikelihoodRatio.jl/dev/)
[![Build Status](https://github.com/kchu25/AverageLoglikelihoodRatio.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/AverageLoglikelihoodRatio.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/AverageLoglikelihoodRatio.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/AverageLoglikelihoodRatio.jl)

Compare position frequency matrices (PFMs) using the Average Log-Likelihood Ratio (ALLR) from Wang and Stormo (2003). Computes a similarity score between two PFMs of possibly different widths, with optional p-values via permutation test.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/kchu25/AverageLoglikelihoodRatio.jl")
```

## Usage

### Single pair

```julia
using AverageLoglikelihoodRatio

# Each PFM is 4 x W (rows = A, C, G, T; columns = positions)
input = [0.90 0.05 0.05;
         0.05 0.85 0.05;
         0.025 0.05 0.85;
         0.025 0.05 0.05]

target = [0.25 0.88 0.06 0.06 0.25;
          0.25 0.04 0.84 0.09 0.25;
          0.25 0.04 0.05 0.78 0.25;
          0.25 0.04 0.05 0.07 0.25]

result = compare_pfms(input, target)
result.score   # best ALLR across all alignments
result.offset  # 0-indexed position of shorter PFM on longer PFM
result.pvalue  # permutation p-value (1000 permutations by default)
```

### Multiple inputs vs. multiple targets

```julia
inputs  = [pfm_a, pfm_b, pfm_c]
targets = [pfm_x, pfm_y]

# Returns a 3 x 2 Matrix{ALLRResult}
results = compare_pfms(inputs, targets)
results[2, 1]  # comparison of pfm_b vs pfm_x
```

### Options

```julia
compare_pfms(input, target;
    background     = [0.25, 0.25, 0.25, 0.25],  # nucleotide background frequencies
    n_perm         = 1000,                        # number of permutations for p-value
    compute_pvalue = true,                        # set false to skip permutation test
    rng            = Random.default_rng()         # RNG for reproducibility
)
```

## API

| Function | Description |
|---|---|
| `compare_pfms(input, target; ...)` | Compare two PFMs. Returns `ALLRResult`. |
| `compare_pfms(inputs, targets; ...)` | Compare all pairs. Returns `Matrix{ALLRResult}`. |
| `allr_column(c1, c2; background)` | ALLR between two length-4 frequency vectors. |
| `allr_matrix(pfm1, pfm2; background)` | Mean column-wise ALLR for equal-width PFMs. |

## How it works

The shorter PFM is slid across every valid offset of the longer PFM. At each offset, the mean ALLR across aligned columns is computed. The offset with the highest score is reported.

For the p-value, `n_perm` random PFMs are drawn from a Dirichlet distribution parameterized by the background frequencies. The observed score is compared against this null distribution:

```
P = (count(S_random >= S_observed) + 1) / (N_permutations + 1)
```

See `latex/allr.tex` for the full formulas.

## Reference

Wang T, Stormo GD. "Combining phylogenetic data with co-regulated genes to identify regulatory motifs." *Bioinformatics*, 19(18):2369-2380, 2003. https://academic.oup.com/bioinformatics/article/19/18/2369/194379
