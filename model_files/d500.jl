# programs/50-base-pk/pumas/d500.jl
#
# d500 — Phase C base model. Locked structural form for the StarterSADMAD
# popPK base model handed off to Stage 60.
#
# Builds on d210 (2-cmt FO + dose-saturable F) and adds:
#   - Michaelis-Menten elimination on the central compartment
#       dCentral/dt -= Vmax · Cp / (Km + Cp)        Cp = 1000·Central/Vc
#     Required to capture MAD Day-14 BID super-proportionality the
#     publication reports (Cmax,ss 100→300 BID ratio ≈ 5.5×).
#   - IOV on Ka (ηka_oc[OCC]) across two occasions: SAD/FE Period 1 vs
#     MAD Day-14 / FE Period 2.  Encoded the canonical Pumas way — a
#     single shared variance `Ωocc` broadcast across an
#     `MvNormal(Diagonal(fill(Ωocc, 2)))` and indexed by `OCC` in `@pre`.
#     Mirror of the simulation truth (`programs/00-dataprep/pumas/true-model.jl`).
#     IOV on F is deferred (literature suggests it's small for this
#     compound class; can be re-added in covariate stage if residuals
#     demand it).
#   - Food covariate slots are NOT included — Phase B Δ-2LL did not
#     justify them (see Phase B section of Stage 50 QMD).
#
# Population: 70 unique post-SAD-005 participants. FE crossover is now
# encoded as one Subject per participant carrying both periods, so the
# Subject count equals the participant count. See
# `programs/00-dataprep/drivers/load-poppk-population.jl` for derivation.

d500 = @model begin
    @metadata begin
        desc = "2-cmt FO absorption + saturable F + MM elimination + IOV on Ka — Phase C base model"
        timeu = u"hr"
    end

    @param begin
        # --- Disposition (linear-equivalent CL is now Vmax/Km) ----------------
        tvvmax ∈ RealDomain(lower = 0)               # max metabolic rate (mg/h)
        tvkm ∈ RealDomain(lower = 0)               # half-saturation (ng/mL)
        tvvc ∈ RealDomain(lower = 0)               # central V/F (L)
        tvq ∈ RealDomain(lower = 0)               # inter-compartmental Q/F (L/h)
        tvvp ∈ RealDomain(lower = 0)               # peripheral V/F (L)

        # --- Absorption ------------------------------------------------------
        tvka ∈ RealDomain(lower = 0)               # 1/h

        # --- Saturable presystemic F -----------------------------------------
        tvfmax ∈ RealDomain(lower = 0, upper = 1)
        tvkf50 ∈ RealDomain(lower = 0)

        # --- IIV on Vmax / Vc / Q / Vp / Ka / Fmax ---------------------------
        Ω ∈ PDiagDomain(6)

        # --- IOV on Ka -------------------------------------------------------
        Ωocc ∈ RealDomain(lower = 0)

        # --- RUV -------------------------------------------------------------
        σ_prop ∈ RealDomain(lower = 0)
        σ_add ∈ RealDomain(lower = 0)
    end

    @random begin
        η ~ MvNormal(Ω)
        # IOV on Ka — canonical Pumas pattern.  A single shared variance
        # `Ωocc` is broadcast across a 2-component MvNormal so each
        # subject draws one ηka_oc per occasion; `OCC` (covariate) picks
        # the active component in `@pre`.  Two occasions cover the SAD /
        # FE Period 1 vs MAD Day-14 / FE Period 2 design.  Mirrors the
        # simulation truth (`programs/00-dataprep/pumas/true-model.jl`).
        ηka_oc ~ MvNormal(Diagonal(fill(Ωocc, 2)))
    end

    @covariates DOSE_MG OCC

    @pre begin
        Vmax = tvvmax * exp(η[1])
        Km = tvkm
        Vc = tvvc * exp(η[2])
        Q = tvq * exp(η[3])
        Vp = tvvp * exp(η[4])
        Ka = tvka * exp(η[5] + ηka_oc[OCC])
    end

    @dosecontrol begin
        bioav = (Depot = tvfmax * DOSE_MG / (tvkf50 + DOSE_MG) * exp(η[6]),)
    end

    @vars begin
        # Plasma concentration in ng/mL. Single named alias used by both
        # the @dynamics MM term (Vmax·cp / (Km + cp)) and @derived. Km
        # is in ng/mL; the unit conversion lives here so it appears once.
        cp = 1000 * Central / Vc
    end

    @dynamics begin
        Depot' = -Ka * Depot
        Central' =
            Ka * Depot - (Vmax * cp / (Km + cp)) - (Q / Vc) * Central +
            (Q / Vp) * Peripheral
        Peripheral' = (Q / Vc) * Central - (Q / Vp) * Peripheral
    end

    @derived begin
        # Pumas 2.8 canonical combined error: σ_add (ng/mL) + proportional σ_prop.
        DV ~ @. CombinedNormal(cp, σ_add, σ_prop)
    end
end

# Initial estimates: structural Vmax/Km seeded from simulation truth
# (truth Vmax = 33.6 mg/h, Km = 800 ng/mL); other parameters from the
# Phase B fit (`fit_d210_sadfe`).
d500_init = (
    tvvmax = 35.0,
    tvkm = 800.0,
    tvvc = 350.0,
    tvq = 18.0,
    tvvp = 420.0,
    tvka = 0.52,
    tvfmax = 0.71,
    tvkf50 = 149.0,
    Ω = Diagonal([0.07, 0.07, 0.11, 0.05, 0.07, 0.10]),
    Ωocc = 0.05,
    σ_prop = 0.10,
    σ_add = 0.42,
)
