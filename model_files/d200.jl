# programs/50-base-pk/pumas/d200.jl
#
# d200 — 2-compartment, first-order absorption, constant bioavailability.
# Industry-standard starting point for a small-molecule oral PK program.
# Should improve markedly on d100's terminal-phase mis-fit but cannot
# capture the SAD supra-proportionality (constant F).

d200 = @model begin
    @metadata begin
        desc = "2-cmt FO absorption + constant F (Phase A baseline)"
        timeu = u"hr"
    end

    @param begin
        tvcl ∈ RealDomain(lower = 0)
        tvvc ∈ RealDomain(lower = 0)
        tvq ∈ RealDomain(lower = 0)
        tvvp ∈ RealDomain(lower = 0)
        tvka ∈ RealDomain(lower = 0)

        Ω ∈ PDiagDomain(5)        # CL, Vc, Q, Vp, Ka
        σ_prop ∈ RealDomain(lower = 0)
        σ_add ∈ RealDomain(lower = 0)
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = tvcl * exp(η[1])
        Vc = tvvc * exp(η[2])
        Q = tvq * exp(η[3])
        Vp = tvvp * exp(η[4])
        Ka = tvka * exp(η[5])
    end

    @vars begin
        # Plasma concentration in ng/mL. Single named alias used by @derived.
        cp = 1000 * Central / Vc
    end

    @dynamics Depots1Central1Periph1   # 2-cmt + first-order absorption

    @derived begin
        # Pumas 2.8 canonical combined error: σ_add (ng/mL) + proportional σ_prop.
        DV ~ @. CombinedNormal(cp, σ_add, σ_prop)
    end
end

d200_init = (
    tvcl = 35.0,
    tvvc = 200.0,
    tvq = 25.0,
    tvvp = 600.0,
    tvka = 0.6,
    Ω = Diagonal([0.1, 0.1, 0.1, 0.1, 0.2]),
    σ_prop = 0.20,
    σ_add = 0.5,
)
