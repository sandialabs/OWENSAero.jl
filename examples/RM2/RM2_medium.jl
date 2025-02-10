using Plots
using DelimitedFiles: readdlm
import OWENSAero: setupTurb, steadyTurb
import FLOWMath: akima

# Turbine
airfoils_path = abspath("./airfoils")
airfoil = "NACA_0021.dat"
radius = 0.538
height = 0.807
chordmid = 0.066667
chordtip = 0.04
eta = 0.5
B = 3
ω = NaN # will vary this parameter later
ntheta = 30
Nslices = 30
afname = joinpath(airfoils_path, airfoil)
chord = akima([0.0, 0.5, 1.0], [chordtip, chordmid, chordtip], LinRange(0, 1, Nslices))
shapeX = ones(Nslices) .* radius
shapeZ = LinRange(0, height, Nslices)
area = height * 2radius

# Fluid
Vinf = 1.2
rho = 1000.0
mu = 1.792E-3

# Model
ifw = false
# AeroModel = "DMS"  # AeroModel ∈ ["DMS", "AC"]
# Aero_AddedMass_Active = true
# Aero_Buoyancy_Active = false
Aero_RotAccel_Active = true
DynamicStallModel = "BV"

results = Dict{String,Dict{String,Array}}()

min_tipspeed, max_tipspeed = 1, 5
n_tipspeeds = 15
tsr = range(min_tipspeed, max_tipspeed, n_tipspeeds)

for AeroModel ∈ ["DMS", "AC"], Aero_AddedMass_Active ∈ [false, true], Aero_Buoyancy_Active ∈ [false, true]
    # Setup
    _ = setupTurb(shapeX, shapeZ, B, chord, ω, Vinf;
        rho, mu, eta, afname, DynamicStallModel, Nslices, Aero_AddedMass_Active, AeroModel,
        Aero_RotAccel_Active, Aero_Buoyancy_Active
    )

    # Solve problem for range of tip speed ratios
    CP = Vector{Float64}(undef, n_tipspeeds)
    RpSteady = Array{Float64}(undef, (B, Nslices, ntheta, n_tipspeeds))
    TpSteady = similar(RpSteady)
    alphaSteady = similar(RpSteady)

    label = AeroModel * (Aero_AddedMass_Active ? "-AM" : "") * (Aero_Buoyancy_Active ? "-BY" : "")
    println("$(label)")
    for (iTSR, TSR) in enumerate(tsr)
        println("  TSR: $(round(TSR; digits=2))")
        omega = Vinf / radius * TSR
        i_results = steadyTurb(; omega, Vinf)

        RpSteady[:, :, :, iTSR] = i_results[2]
        TpSteady[:, :, :, iTSR] = i_results[3]
        alphaSteady[:, :, :, iTSR] = i_results[5]
        CP[iTSR] = i_results[18] / (0.5rho * Vinf^3 * area)
    end
    println("")
    results[label] = Dict(
        "cₚ" => CP,
        "rₚ" => RpSteady,
        "tₚ" => TpSteady,
        "α" => alphaSteady,
    )
end

# load experimental data
expdata_path = abspath("./exp_data")
exp_0p9 = readdlm(joinpath(expdata_path, "RM2_0.538D_RE_D_0.9E6.csv"), ',', Float64)
exp_1p3 = readdlm(joinpath(expdata_path, "RM2_0.538D_RE_D_1.3E6.csv"), ',', Float64)

# Plot
fig_path = abspath("./figures")
mkpath(fig_path)

p = plot(xlabel="tip speed ratio", ylabel="Cₚ")
for AeroModel ∈ ["DMS", "AC"], Aero_AddedMass_Active ∈ [false, true], Aero_Buoyancy_Active ∈ [false, true]
    if !(!Aero_AddedMass_Active && Aero_Buoyancy_Active)
        label = AeroModel * (Aero_AddedMass_Active ? "-AM" : "") * (Aero_Buoyancy_Active ? "-BY" : "")
        ls = Aero_Buoyancy_Active ? :dot : (Aero_AddedMass_Active ? :dash : :solid)
        lc = AeroModel == "DMS" ? :blue : :red
        plot!(tsr, results[label]["cₚ"]; ls, lc, label)
    end
end
# scatter!(exp_0p9[:, 1], exp_0p9[:, 2]; mc=:black, markershape=:circle, label="Exp. Re=0.9")
scatter!(exp_1p3[:, 1], exp_1p3[:, 2]; mc=:black, markershape=:diamond, label="Exp. Re=1.3")
savefig(p, joinpath(fig_path, "Cp.pdf"))

function plot_slice(islice, itsr, iblade=1, results=results)
    p_r = plot(xlabel="θ [°]", ylabel="rₚ")
    p_t = plot(xlabel="θ [°]", ylabel="tₚ")
    p_α = plot(xlabel="θ [°]", ylabel="α")

    θ = (1:ntheta) ./ ntheta .* 360

    for AeroModel ∈ ["DMS", "AC"], Aero_AddedMass_Active ∈ [false, true], Aero_Buoyancy_Active ∈ [false, true]
        if !(!Aero_AddedMass_Active && Aero_Buoyancy_Active)
            label = AeroModel * (Aero_AddedMass_Active ? "-AM" : "") * (Aero_Buoyancy_Active ? "-BY" : "")
            rₚ = results[label]["rₚ"][iblade, islice, :, itsr]
            tₚ = results[label]["tₚ"][iblade, islice, :, itsr]
            α = results[label]["α"][iblade, islice, :, itsr]
            ls = Aero_Buoyancy_Active ? :dot : (Aero_AddedMass_Active ? :dash : :solid)
            lc = AeroModel == "DMS" ? :blue : :red
            plot!(p_r, θ, rₚ; ls, lc, label)
            plot!(p_t, θ, tₚ; ls, lc, label)
            plot!(p_α, θ, α; ls, lc, label)
        end
    end
    savefig(p_r, joinpath(fig_path, "rp_$(islice)_$(itsr)_$(iblade).pdf"))
    savefig(p_t, joinpath(fig_path, "tp_$(islice)_$(itsr)_$(iblade).pdf"))
    savefig(p_α, joinpath(fig_path, "alpha_$(islice)_$(itsr)_$(iblade).pdf"))
end

islice = 15
itsr = 7

plot_slice(islice, itsr)
