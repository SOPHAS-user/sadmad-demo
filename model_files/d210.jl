# programs/50-base-pk/pumas/d210.jl
#
# d210 — 2-compartment, first-order absorption, dose-saturable
# bioavailability F = Fmax · Dose / (KF50 + Dose). The publication
# explicitly attributes SAD supra-proportionality to "saturable presystemic
# restriction to bioavailability"; d210 is the smallest model that
# encodes that mechanism.

d210 = @model begin
    @metadata begin
        desc = "2-cmt FO absorption + saturable presystemic F (Phase A non-linear-F candidate)"
        timeu = u"hr"
    end

    @param begin
        tvcl ∈ RealDomain(lower = 0)
        tvvc ∈ RealDomain(lower = 0)
        tvq ∈ RealDomain(lower = 0)
        tvvp ∈ RealDomain(lower = 0)
        tvka ∈ RealDomain(lower = 0)
        tvfmax ∈ RealDomain(lower = 0, upper = 1)
        tvkf50 ∈ RealDomain(lower = 0)

        Ω ∈ PDiagDomain(6)        # CL, Vc, Q, Vp, Ka, Fmax
        σ_prop ∈ RealDomain(lower = 0)
        σ_add ∈ RealDomain(lower = 0)
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates DOSE_MG

    @pre begin
        CL = tvcl * exp(η[1])
        Vc = tvvc * exp(η[2])
        Q = tvq * exp(η[3])
        Vp = tvvp * exp(η[4])
        Ka = tvka * exp(η[5])
    end

    @dosecontrol begin
        # Saturable presystemic F = Fmax · Dose / (KF50 + Dose). η[6] is
        # subject-level deviation in Fmax (log-normal).
        bioav = (Depot = tvfmax * DOSE_MG / (tvkf50 + DOSE_MG) * exp(η[6]),)
    end

    @vars begin
        # Plasma concentration in ng/mL. Single named alias used by @derived.
        cp = 1000 * Central / Vc
    end

    @dynamics Depots1Central1Periph1

    @derived begin
        # Pumas 2.8 canonical combined error: σ_add (ng/mL) + proportional σ_prop.
        DV ~ @. CombinedNormal(cp, σ_add, σ_prop)
    end
end

d210_init = (
    tvcl = 35.0,
    tvvc = 200.0,
    tvq = 25.0,
    tvvp = 600.0,
    tvka = 0.6,
    tvfmax = 0.7,
    tvkf50 = 100.0,
    Ω = Diagonal([0.1, 0.1, 0.1, 0.1, 0.2, 0.1]),
    σ_prop = 0.20,
    σ_add = 0.5,
)
