include("models.jl")
using LsqFit, Plots, Printf, LaTeXStrings

default(
  fontfamily = "Times Roman",
  guidefontsize = 12,
  tickfontsize = 10,
  legendfontsize = 10,
  markerstrokewidth = 2,
  margin = 5Plots.mm,
  size = (600, 400),
  grid = false,
  foreground_color_legend = nothing,
  background_color_legend = colorant"rgba(255,255,255,0.6)",
  linewidth = 1
)

"""
Do calculation and plot the results.
## Parameters
`χ`: max bond dimension of HOTRG.
`stepnum`: step number of HOTRG.
`eigvalnum`: number of operators to calculate one-point function.
"""
function main(;χ, stepnum, eigvalnum)
  #
  # HOTRG calculation
  # It takes about 3 seconds with χ = 16 and 35 seconds with χ = 24 for each calculation on my laptop (without compilation time).
  #
  @time _, cftdata_ising_sro = hotrg(bulk(weight(Ising(), inv(Tc(Ising()))))...; maxdim = χ, stepnum, eigvalnum)
  @time _, cftdata_ising_refl = hotrg(bulk(weight(Ising(), inv(Tc(Ising()))))...; maxdim = χ, stepnum, eigvalnum, withsro = false)
  @time _, cftdata_potts3_sro = hotrg(bulk(weight(Potts(3), inv(Tc(Potts(3)))))...; maxdim = χ, stepnum, eigvalnum)
  @time _, cftdata_potts3_refl = hotrg(bulk(weight(Potts(3), inv(Tc(Potts(3)))))...; maxdim = χ, stepnum, eigvalnum, withsro = false)


  #
  # plot rainbow free energy term
  #
  fr_ising = log.(abs.(cftdata_ising_sro["<R|i>"][1, :]))
  fr_ising_reflect = log.(abs.(cftdata_ising_refl["<R|i>"][1, :]))
  fr_potts3 = log.(abs.(cftdata_potts3_sro["<R|i>"][1, :]))
  fr_potts3_reflect = log.(abs.(cftdata_potts3_refl["<R|i>"][1, :]))
  β = 2 .^ (1:stepnum)

  fitrange_ising = 1:6
  fitrange_potts3 = 1:3

  f_ising(x, p) = @. (centralcharge(Ising()) / 4) * log(x) + p[1]
  fit_ising = curve_fit(f_ising, β[fitrange_ising], fr_ising[fitrange_ising], [0.])
  f_ising(x) = f_ising(x, fit_ising.param)
  f_potts3(x, p) = @. (centralcharge(Potts(3)) / 4) * log(x) + p[1]
  fit_potts3 = curve_fit(f_potts3, β[fitrange_potts3], fr_potts3[fitrange_potts3], [0.])
  f_potts3(x) = f_potts3(x, fit_potts3.param)

  p_fr = plot()
  scatter!(p_fr, β, fr_ising; xscale = :log2, xlabel = L"\beta", ylabel = L"F_{R}", legend = :topleft,
    label = "Ising", markershape = :x, markercolor = :darkcyan, xticks = [4, 16, 64] |> x -> (x, string.(x)), format = :png)
  scatter!(p_fr, β, fr_ising_reflect; label = "Ising (reflect)", markershape = :vline, markercolor = :darkcyan)
  scatter!(p_fr, β, fr_potts3; label = "Potts3", markershape = :+, markercolor = :lightcoral)
  scatter!(p_fr, β, fr_potts3_reflect; label = "Potts3", markershape = :x, markercolor = :lightcoral)
  plot!(p_fr, β, f_ising(β); line = :dash, linecolor = :orange, label = L"F_{R} = \frac{0.5}{4}\ln\beta + b")
  plot!(p_fr, β, f_potts3(β); line = :dot, linecolor = :mediumorchid, label = L"F_{R} = \frac{0.8}{4}\ln\beta + b")
  p_fr[:dpi] = round(Int, 4 * p_fr[:dpi])
  savefig(p_fr, "fig/fr.png")


  #
  # plot crosscap free energy term
  #
  fc_ising = log.(abs.(cftdata_ising_sro["<C|i>"][1, :]))
  fc_potts3 = log.(abs.(cftdata_potts3_sro["<C|i>"][1, :]))

  p_fc = plot()
  scatter!(p_fc, β, fc_ising; xscale = :log2, xlabel = L"\beta", ylabel = L"F_{C}",
    label = "Ising", markershape = :x, markercolor = :darkcyan, xticks = [4, 16, 64] |> x -> (x, string.(x)), legend = :right)
  scatter!(p_fc, β, fc_potts3; label = "Potts3", markershape = :+, markercolor = :lightcoral, format = :png)
  hline!(p_fc, [log(quantumdimension(Ising())) / 2]; line = :dash, linecolor = :orange, label = L"\frac{1}{2}\ln g_{\mathrm{Ising}}")
  hline!(p_fc, [log(quantumdimension(Potts(3))) / 2]; line = :dot, linecolor = :mediumorchid, label = L"\frac{1}{2}\ln g_{\mathrm{Potts3}}")
  p_fc[:dpi] = round(Int, 4 * p_fc[:dpi])
  savefig(p_fc, "fig/fc.png")


  #
  # plot one-point function of the Ising CFT
  #
  gamma2_ising = cftdata_ising_sro["<C|i>"][1:3, :] .^ 2

  p_gamma = plot()
  plot!(p_gamma, β, [gamma2_ising[i, :] for i in 1:3]; xscale = :log2, xlabel = L"\beta", ylabel = L"\Gamma_k^2",
    label = nothing, markershape = :x, color = :darkcyan, xticks = [4, 16, 64] |> x -> (x, string.(x)), legend = :right)
  hline!(p_gamma, [1 + 1 / √2], label = L"\Gamma_\mathrm{I}^2 = 1 + 1 / \sqrt{2}")
  hline!(p_gamma, [1 - 1 / √2], label = L"\Gamma_{\epsilon}^2 = 1 - 1 / \sqrt{2}")
  hline!(p_gamma, [0], label = L"\Gamma_{\sigma}^2 = 0")
  p_gamma[:dpi] = round(Int, 4 * p_gamma[:dpi])
  savefig(p_gamma, "fig/gamma2.png")
end

# main(χ = 16, stepnum = 7, eigvalnum = 3)

function impltrg(;χ, stepnum, eigvalnum, relT)
  f = []
  βs = inv.(Tc(Ising()) * relT)
  for β in βs
    norms, eigval = gilttnr(bulk(weight(Ising(), β))...; maxdim = χ, stepnum, eigvalnum)
    lnz = logpartfunc(norms, sitenum_per_step = 2)
    push!(f, lnz[end] / -β)
  end

  scatter(βs, f)
  xs = [minimum(βs):1e-3:maximum(βs);]
  plot!(xs, xs .|> x -> freeenergy(Ising(), x))
  savefig("fig/f.png")
end

function err4ϵ(;χ, stepnum, ϵs)
  f = []
  β = inv(Tc(Ising()))
  for ϵ in ϵs
    norms, _ = gilttnr(bulk(weight(Ising(), β))...; maxdim = χ, stepnum, eigvalnum = 1, ϵ)
    lnz = logpartfunc(norms, sitenum_per_step = 2)
    push!(f, lnz / -β)
  end
  plot()
  f_exact = freeenergy(Ising(), β)
  for i in eachindex(ϵs)
    plot!(abs.(f[i] .- f_exact), yscale = :log10, label = "\epsilon = $(ϵs[i])", xlabel = "RG step", ylabel = "error of free energy"; marker = :auto)
  end

  norms, _ = trg(bulk(weight(Ising(), β))...; maxdim = χ, stepnum, eigvalnum = 1)
  lnz = logpartfunc(norms, sitenum_per_step = 2)
  f_trg = lnz / -β
  plot!(abs.(f_trg .- f_exact), label = "trg")
  open("data/ferr.txt", "w") do fp
    for i in eachindex(ϵs)
      Base.print_array(fp, (f[i] .- f_exact))
    end
    Base.print_array(fp, (f_trg .- f_exact))
  end
  savefig("fig/ferr.png")
end

# @time impltrg(χ = 24, stepnum = 15, eigvalnum = 3, relT = [0.95:1e-2:1.05;])
@time err4ϵ(χ = 24, stepnum = 25, ϵs = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10])