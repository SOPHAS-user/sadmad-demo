# programs/50-base-pk/pumas/d100.jl
#
# d100 — 1-compartment, first-order absorption, constant bioavailability.
# Anchor / failure case: the publication's apparent t1/2z (22–30 h) and
# Vd (>1000 L) demand a 2-cmt structure. d100 is fit so its OFV / GOF can
# document the structural deficiency before moving on.
#
# Random-effect parameterization: SCALAR Normal etas. Each η is its own
# named univariate random effect with its own SD parameter (`ω_cl`, `ω_vc`,
# `ω_ka`). This is mathematically identical to the diagonal-Ω form
# (`PDiagDomain(3)` + `MvNormal(Ω)`) used in d200 onward — same likelihood,
# same optimum — but it makes each η's source explicit and is the cleanest
# starting point for learners. See the "Random effects — three styles"
# callout in 50-base-pk.qmd for the full comparison with the
# `PDiagDomain` (d200/d210/d211/d500) and `PSDDomain` (d600) forms.

d100 = @model begin
    @metadata begin
        desc = "1-cmt FO absorption + constant F (Phase A anchor)"
        timeu = u"hr"
    end

    @param begin
        tvcl ∈ RealDomain(lower = 0)
        tvvc ∈ RealDomain(lower = 0)
        tvka ∈ RealDomain(lower = 0)

        # Scalar IIV SDs — each named, each its own RealDomain. Contrast
        # with d200's `Ω ∈ PDiagDomain(5)` which stores the same numbers
        # in a diagonal matrix indexed `Ω[i,i]`.
        ω_cl ∈ RealDomain(lower = 0)
        ω_vc ∈ RealDomain(lower = 0)
        ω_ka ∈ RealDomain(lower = 0)

        σ_prop ∈ RealDomain(lower = 0)
        σ_add ∈ RealDomain(lower = 0)
    end

    @random begin
        # Scalar form: one Normal per parameter, addressed by name in @pre.
        # `Normal(0, ω)` takes the SD directly (not the variance).
        η_cl ~ Normal(0.0, ω_cl)
        η_vc ~ Normal(0.0, ω_vc)
        η_ka ~ Normal(0.0, ω_ka)
    end

    @pre begin
        CL = tvcl * exp(η_cl)
        Vc = tvvc * exp(η_vc)
        Ka = tvka * exp(η_ka)
    end

    @vars begin
        # Plasma concentration in ng/mL. Shared across @derived (and used
        # by @dynamics in models with MM elimination — see d500/d600).
        cp = 1000 * Central / Vc
    end

    @dynamics Depots1Central1     # 1-cmt + first-order absorption from Depot

    @derived begin
        # Pumas 2.8 canonical combined error: σ_add (ng/mL) + proportional σ_prop.
        DV ~ @. CombinedNormal(cp, σ_add, σ_prop)
    end
end

# Initial estimates: SDs (not variances). Square-rooted from the prior
# diagonal-Ω form to preserve the same effective IIV magnitude.
d100_init = (
    tvcl = 35.0,
    tvvc = 200.0,
    tvka = 0.6,
    ω_cl = sqrt(0.1),
    ω_vc = sqrt(0.1),
    ω_ka = sqrt(0.2),
    σ_prop = 0.20,
    σ_add = 0.5,
)
