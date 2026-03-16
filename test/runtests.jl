using AverageLoglikelihoodRatio
using Test
using Random

@testset "AverageLoglikelihoodRatio.jl" begin

    # ── Column-level ALLR ──────────────────────────────────────────────────
    @testset "allr_column" begin
        # Identical columns → maximum similarity (positive score)
        c = [0.7, 0.1, 0.1, 0.1]
        sc = allr_column(c, c)
        @test sc > 0.0

        # Uniform columns → score ≈ 0 (no information above background)
        u = [0.25, 0.25, 0.25, 0.25]
        sc_uni = allr_column(u, u)
        @test abs(sc_uni) < 1e-8

        # Opposite columns → low / negative score
        c1 = [0.9, 0.033, 0.034, 0.033]
        c2 = [0.033, 0.034, 0.033, 0.9]
        sc_opp = allr_column(c1, c2)
        @test sc_opp < sc   # much worse than identical

        # Symmetry: ALLR(c1, c2) == ALLR(c2, c1)
        @test allr_column(c1, c2) ≈ allr_column(c2, c1)

        # Custom background
        bg = [0.3, 0.2, 0.2, 0.3]
        sc_bg = allr_column(c, c; background=bg)
        @test isfinite(sc_bg)

        # Wrong length → error
        @test_throws ArgumentError allr_column([0.5, 0.5], c)
    end

    # ── Matrix-level ALLR (aligned) ────────────────────────────────────────
    @testset "allr_matrix" begin
        pfm = [0.7 0.1 0.1;
               0.1 0.7 0.1;
               0.1 0.1 0.7;
               0.1 0.1 0.1]

        # Same matrix → positive score
        sc = allr_matrix(pfm, pfm)
        @test sc > 0.0

        # Dimension mismatch → error
        pfm_short = pfm[:, 1:2]
        @test_throws ArgumentError allr_matrix(pfm, pfm_short)

        # 3-row matrix → error
        @test_throws ArgumentError allr_matrix(pfm[1:3, :], pfm[1:3, :])
    end

    # ── Sliding alignment via compare_pfms ─────────────────────────────────
    @testset "compare_pfms – alignment" begin
        # Build a target that contains the input as a sub-motif at offset 2
        motif = [0.9 0.05 0.05;
                 0.05 0.9 0.05;
                 0.025 0.025 0.85;
                 0.025 0.025 0.05]

        bg_col = [0.25, 0.25, 0.25, 0.25]
        # Pad with background columns on each side
        bg_mat = repeat(reshape(bg_col, 4, 1), 1, 2)
        target = hcat(bg_mat, motif, bg_mat)   # 4 × 7, motif at cols 3-5 (offset 2)

        res = compare_pfms(motif, target; compute_pvalue=false)
        @test res.offset == 2
        @test res.score > 0.0
        @test isnan(res.pvalue)
    end

    # ── compare_pfms: input longer than target ─────────────────────────────
    @testset "compare_pfms – input longer" begin
        short = [0.8 0.1;
                 0.1 0.8;
                 0.05 0.05;
                 0.05 0.05]

        long = [0.25 0.8 0.1 0.25;
                0.25 0.1 0.8 0.25;
                0.25 0.05 0.05 0.25;
                0.25 0.05 0.05 0.25]

        res = compare_pfms(long, short; compute_pvalue=false)
        @test res.score > 0.0
        @test 0 ≤ res.offset ≤ 2   # offset within valid range
    end

    # ── Equal-width matrices ───────────────────────────────────────────────
    @testset "compare_pfms – equal width" begin
        pfm = [0.85 0.15;
               0.05 0.75;
               0.05 0.05;
               0.05 0.05]
        res = compare_pfms(pfm, pfm; compute_pvalue=false)
        @test res.score > 0.0
        @test res.offset == 0
    end

    # ── Permutation p-value ────────────────────────────────────────────────
    @testset "compare_pfms – p-value" begin
        rng = MersenneTwister(42)

        # Two very similar motifs → small p-value
        pfm1 = [0.9 0.05 0.05;
                0.05 0.85 0.1;
                0.025 0.05 0.8;
                0.025 0.05 0.05]

        pfm2 = [0.88 0.06 0.06;
                0.04 0.84 0.09;
                0.04 0.05 0.78;
                0.04 0.05 0.07]

        res = compare_pfms(pfm1, pfm2; n_perm=500, rng=rng)
        @test 0.0 < res.pvalue ≤ 1.0
        @test res.pvalue < 0.1   # should be quite significant

        # Two dissimilar motifs → large p-value (not significant)
        rng2 = MersenneTwister(99)
        # One motif is A-rich, the other is T-rich → poor match
        pfm_a = [0.85 0.85 0.85;
                 0.05 0.05 0.05;
                 0.05 0.05 0.05;
                 0.05 0.05 0.05]
        pfm_t = [0.05 0.05 0.05;
                 0.05 0.05 0.05;
                 0.05 0.05 0.05;
                 0.85 0.85 0.85]
        res2 = compare_pfms(pfm_a, pfm_t; n_perm=500, rng=rng2)
        @test res2.pvalue > 0.3   # should not be significant
    end

    # ── Input validation ───────────────────────────────────────────────────
    @testset "input validation" begin
        @test_throws ArgumentError compare_pfms(rand(3, 4), rand(4, 4))
        @test_throws ArgumentError compare_pfms(rand(4, 4), rand(5, 4))
        @test_throws ArgumentError compare_pfms(rand(4, 4), rand(4, 4); background=[0.5, 0.5])
    end

    # ── ALLRResult display ─────────────────────────────────────────────────────────
    @testset "ALLRResult show" begin
        r = ALLRResult(1.234, 3, 0.005, false)
        buf = IOBuffer()
        show(buf, r)
        s = String(take!(buf))
        @test occursin("1.234", s)
        @test occursin("offset=3", s)
        @test occursin("0.005", s)
        @test occursin("forward", s)

        r2 = ALLRResult(0.5, 0, NaN, true)
        buf2 = IOBuffer()
        show(buf2, r2)
        s2 = String(take!(buf2))
        @test occursin("not computed", s2)
        @test occursin("reverse", s2)
    end

    # ── Batch comparison: multiple inputs × multiple targets ───────────────
    @testset "compare_pfms – batch" begin
        pfm1 = [0.9 0.05 0.05;
                0.05 0.85 0.05;
                0.025 0.05 0.85;
                0.025 0.05 0.05]

        pfm2 = [0.88 0.06;
                0.04 0.84;
                0.04 0.05;
                0.04 0.05]

        pfm3 = [0.05 0.05 0.85 0.25;
                0.85 0.05 0.05 0.25;
                0.05 0.85 0.05 0.25;
                0.05 0.05 0.05 0.25]

        inputs  = [pfm1, pfm2]
        targets = [pfm2, pfm3, pfm1]

        # Without p-value (fast)
        results = compare_pfms(inputs, targets; compute_pvalue=false)
        @test size(results) == (2, 3)
        @test all(r -> isa(r, ALLRResult), results)
        @test all(r -> isnan(r.pvalue), results)

        # Each entry should match the single-pair call
        for i in 1:2, j in 1:3
            single = compare_pfms(inputs[i], targets[j]; compute_pvalue=false)
            @test results[i, j].score ≈ single.score
            @test results[i, j].offset == single.offset
        end

        # Self-comparison: pfm vs itself should have high score
        self_results = compare_pfms([pfm1], [pfm1]; compute_pvalue=false)
        @test size(self_results) == (1, 1)
        @test self_results[1, 1].score > 0.0
        @test self_results[1, 1].offset == 0

        # With p-value
        rng = MersenneTwister(123)
        results_pv = compare_pfms([pfm1], [pfm1]; n_perm=200, rng=rng)
        @test size(results_pv) == (1, 1)
        @test 0.0 < results_pv[1, 1].pvalue ≤ 1.0
    end

    # ── Reverse complement ──────────────────────────────────────────────────────
    @testset "compare_pfms – reverse complement" begin
        # A motif and its reverse complement should be a perfect match
        # when reverse_comp=true
        motif = [0.9 0.05 0.05;
                 0.05 0.85 0.05;
                 0.025 0.05 0.85;
                 0.025 0.05 0.05]

        # Manually build the reverse complement:
        # swap rows A↔T (1↔4), C↔G (2↔3), then reverse columns
        motif_rc = motif[[4, 3, 2, 1], end:-1:1]

        # Without reverse_comp: comparing motif to its RC on forward strand
        res_fwd = compare_pfms(motif, motif_rc; compute_pvalue=false)
        @test res_fwd.reverse_complement == false

        # With reverse_comp: should find the RC match and score higher or equal
        res_rc = compare_pfms(motif, motif_rc; compute_pvalue=false, reverse_comp=true)
        @test res_rc.score ≥ res_fwd.score
        # The RC of the RC is the original, so it should be a perfect match
        @test res_rc.reverse_complement == true

        # Score with reverse_comp=true on the RC target should equal
        # the self-comparison score
        res_self = compare_pfms(motif, motif; compute_pvalue=false)
        @test res_rc.score ≈ res_self.score

        # Forward-only: default reverse_complement field is false
        res_default = compare_pfms(motif, motif; compute_pvalue=false)
        @test res_default.reverse_complement == false

        # Reverse comp with p-value
        rng = MersenneTwister(77)
        res_pv = compare_pfms(motif, motif_rc; n_perm=200, reverse_comp=true, rng=rng)
        @test 0.0 < res_pv.pvalue ≤ 1.0
        @test res_pv.reverse_complement == true

        # Batch with reverse_comp
        results = compare_pfms([motif], [motif_rc]; compute_pvalue=false, reverse_comp=true)
        @test size(results) == (1, 1)
        @test results[1, 1].reverse_complement == true
        @test results[1, 1].score ≈ res_self.score
    end
end
