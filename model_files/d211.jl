# programs/50-base-pk/pumas/d211.jl
#
# d211 — Phase B model: d210 (2-cmt FO + dose-saturable F) plus food
# covariates on F and Ka. The publication's food-effect study at 200 mg
# reports a small fed/fasted Cmax/AUC bump (14–18 % per Table 3) and a
# tmax shift of ~1 h, attributable to (a) a small increase in
# bioavailability when fed and (b) slowed absorption (smaller Ka).
#
# d211 retains all d210 structural parameters and adds:
#   tvfood_f  — multiplicative effect on F when fed (FED == 1)
#   tvfood_ka — multiplicative effect on Ka when fed
# Both default to 1.0 (no food effect); positive food effect → tvfood_f > 1
# and tvfood_ka < 1 per the publication direction.

d211 = @model begin
    @metadata begin
        desc = "2-cmt FO absorption + saturable F + food covariates on F and Ka (Phase B)"
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

        # Food covariates — multiplicative; 1.0 = no effect.
        tvfood_f ∈ RealDomain(lower = 0)
        tvfood_ka ∈ RealDomain(lower = 0)

        Ω ∈ PDiagDomain(6)        # CL, Vc, Q, Vp, Ka, Fmax
        σ_prop ∈ RealDomain(lower = 0)
        σ_add ∈ RealDomain(lower = 0)
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates DOSE_MG FED

    @pre begin
        food_coef_ka = (FED == 1) ? tvfood_ka : 1.0
        CL = tvcl * exp(η[1])
        Vc = tvvc * exp(η[2])
        Q = tvq * exp(η[3])
        Vp = tvvp * exp(η[4])
        Ka = tvka * food_coef_ka * exp(η[5])
    end

    @dosecontrol begin
        # Saturable F · food multiplier (when fed) · log-normal IIV.
        bioav = (
            Depot = tvfmax * DOSE_MG / (tvkf50 + DOSE_MG) *
                    ((FED == 1) ? tvfood_f : 1.0) *
                    exp(η[6]),
        )
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

# Initial estimates: structural parameters seeded from d210; food
# coefficients at 1.0 (null hypothesis "no food effect").
d211_init = (
    tvcl = 35.0,
    tvvc = 357.0,
    tvq = 19.0,
    tvvp = 446.0,
    tvka = 0.52,
    tvfmax = 0.70,
    tvkf50 = 146.0,
    tvfood_f = 1.0,
    tvfood_ka = 1.0,
    Ω = Diagonal([0.08, 0.07, 0.13, 0.05, 0.09, 0.11]),
    σ_prop = 0.10,
    σ_add = 0.42,
)
