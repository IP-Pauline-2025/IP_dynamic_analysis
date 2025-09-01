using Printf, Dates, Plots, JSON, JLD2, UncertaintyQuantification

filename = "results_load_and_plot/20250827_170605_model_1_samples.jld2"

# --- Load JLD2 Data ---
@load filename samples

t = collect(0:samples.dt[1]:samples.T[1])

nmc = rand(1:length(samples.wl_base))  # Randomly select one sample

# Normalized wind speed and displacement for sample nmc
norm_wind_speed = samples.wl_base[nmc] ./ maximum(abs.(samples.wl_base[nmc]))
norm_disp = samples.disp[nmc] ./ maximum(abs.(samples.disp[nmc]))

pnorm = plot(
    samples.sim_time[nmc], norm_disp;
    label="Normalized Displacement",
    xlabel="Time (s)",
    ylabel="Normalized Value",
    legend=:topright,
    linewidth=2
)
plot!(t, norm_wind_speed; label="Normalized Wind Speed", linewidth=2)
display(pnorm)