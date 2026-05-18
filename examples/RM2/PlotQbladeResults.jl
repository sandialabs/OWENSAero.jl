using Plots
Plots.closeall()
using DelimitedFiles: readdlm

path,_ = splitdir(@__FILE__)

expdata_path = "$path/exp_data"
exp_0p9 = readdlm(joinpath(expdata_path, "RM2_0.538D_RE_D_0.9E6.csv"), ',', Float64)
exp_1p3 = readdlm(joinpath(expdata_path, "RM2_0.538D_RE_D_1.3E6.csv"), ',', Float64)

Qbld_0p9_TSR = [2.5,3.0]
Qbld_0p9_CP = [0.33,0.38]

Qbld_1p2_TSR = [2.0,3.0,4.0]
Qbld_1p2_CP = [0.19,0.42,0.35]

timeplot = Plots.scatter(Qbld_1p2_TSR,Qbld_1p2_CP,
xlabel = "TSR",
ylabel = "CP",
label = "QBlade",
markershape = :circle,
legend = true,
# title = "(CF), AOA: $(AOA[iaoa])",
)

Plots.scatter!(timeplot,exp_1p3[:,1],exp_1p3[:,2],
xlabel = "TSR",
ylabel = "CP",
label = "Exp",
legend = true,
markershape = :diamond,
color = :black,
# title = "(CF), AOA: $(AOA[iaoa])",
)

fig_path = joinpath(path, "figures")
mkpath(fig_path)
outfile = joinpath(fig_path, "qblade_experiment_cp.pdf")
Plots.savefig(timeplot, outfile)
println("Wrote $(outfile)")

if get(ENV, "OWENSAERO_DISPLAY_PLOTS", "false") == "true"
    Plots.display(timeplot)
end
