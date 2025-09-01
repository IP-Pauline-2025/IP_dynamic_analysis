############################################################
# Plotting from saved `samples` DataFrame (Model_2)
############################################################
# obtain dataframe from JLD2 file
using JLD2, DataFrames, JSON3, Plots, Statistics, CSV

# ---------------- paths ----------------
results_dir = raw"C:\IP\results\SuS\Model_2\plots"
jld_path    = raw"C:\IP\results\SuS\Model_2\20250828_162150_model_2_samples.jld2"   # <-- update if different
meta_path   = raw"C:\IP\results\SuS\Model_2\20250828_162150_model_2_metadata.json"       # <-- update if different

# ---------------- load DataFrame ----------------
function load_samples_df(jld_path)::DataFrame
    jldopen(jld_path, "r") do f
        for k in ("samples", "df", "dataframe", "results")
            if haskey(f, k) && f[k] isa DataFrame
                return f[k]
            end
        end
        for k in keys(f)
            try
                v = f[k]
                if v isa DataFrame
                    return v
                end
            catch
            end
        end
        error("No DataFrame found inside: $jld_path (expected a variable named `samples`)")
    end
end

samples = load_samples_df(jld_path)
meta = isfile(meta_path) ? JSON3.read(read(meta_path, String)) : nothing

# Quick sanity check (Model_2 columns)
for name in [:sim_time, :wl_base, :wl_node1, :wl_node2, :max_abs_disp, :E, :disp_node1]
    @assert hasproperty(samples, name) "Missing column: $name. Available: $(names(samples))"
end

# ---------------- helpers ----------------
getF(x) = x isa Number ? float(x) :
          hasproperty(x, :F_wp) ? float(getproperty(x, :F_wp)) :
          throw(ArgumentError("Cannot extract numeric value; missing :F_wp on $(typeof(x))."))
to_float_vec(v) = map(getF, v)

function align_to_window(full_vec::AbstractVector, window_len::Int)
    off = length(full_vec) - window_len
    if off >= 0
        @views return full_vec[(off+1):end]
    elseif !isempty(full_vec)
        return vcat(fill(first(full_vec), -off), full_vec)
    else
        return fill(0.0, window_len)
    end
end

# ---------------- columns & setup ----------------
times       = samples.sim_time                   # Vector{Vector{T}}
disp_node1  = samples.disp_node1                 # top-node displacement (as per your script)
wl_base     = samples.wl_base                    # wind speed signals (m/s)
wl_node1    = samples.wl_node1                   # kN
wl_node2    = samples.wl_node2                   # kN
wl_node3    = samples.wl_node3                   # kN
wl_node4    = samples.wl_node4                   # kN
wl_node5    = samples.wl_node5                   # kN
wl_node6    = samples.wl_node6                   # kN   
wl_node7    = samples.wl_node7                   # kN
wl_node8    = samples.wl_node8                   # kN
wl_node9    = samples.wl_node9                   # kN
wl_node10   = samples.wl_node10                  # kN
wl_node11   = samples.wl_node11                  # kN
wl_node12   = samples.wl_node12                  # kN
wl_node13   = samples.wl_node13                  # kN
wl_node14   = samples.wl_node14                  # kN
wl_node15   = samples.wl_node15                  # kN
wl_node16   = samples.wl_node16                  # kN

E_vals      = samples.E
max_abs_disp= samples.max_abs_disp
N_MC        = length(times)


# choose which load node to use in the “all-samples loads” plot
LOAD_NODE_FOR_ALL = 16     # 1..16
load_col = getproperty(samples, Symbol("wl_node$(LOAD_NODE_FOR_ALL)"))

# capacity for failure check (m) — set to match your run
capacity = 0.1

############################################################
# Plotting & summaries
############################################################

# ---------------- choose sample & align ----------------
nmc = min(5, N_MC)
t_n  = collect(times[nmc])     # truncated to t >= 20 s in your extractors
u_n  = collect(disp_node1[nmc])

v_full   = collect(wl_base[nmc])
wl1_full = collect(wl_node1[nmc])   # kN
wl2_full = collect(wl_node2[nmc])   # kN
wl3_full = collect(wl_node3[nmc])   # kN
wl4_full = collect(wl_node4[nmc])   # kN
wl5_full = collect(wl_node5[nmc])   # kN
wl6_full = collect(wl_node6[nmc])   # kN
wl7_full = collect(wl_node7[nmc])   # kN
wl8_full = collect(wl_node8[nmc])   # kN
wl9_full = collect(wl_node9[nmc])   # kN
wl10_full = collect(wl_node10[nmc]) # kN
wl11_full = collect(wl_node11[nmc]) # kN
wl12_full = collect(wl_node12[nmc]) # kN
wl13_full = collect(wl_node13[nmc]) # kN
wl14_full = collect(wl_node14[nmc]) # kN
wl15_full = collect(wl_node15[nmc]) # kN
wl16_full = collect(wl_node16[nmc]) # kN


wl1_full = eltype(wl1_full) <: Number ? wl1_full : to_float_vec(wl1_full)
wl2_full = eltype(wl2_full) <: Number ? wl2_full : to_float_vec(wl2_full)
wl3_full = eltype(wl3_full) <: Number ? wl3_full : to_float_vec(wl3_full)
wl4_full = eltype(wl4_full) <: Number ? wl4_full : to_float_vec(wl4_full)
wl5_full = eltype(wl5_full) <: Number ? wl5_full : to_float_vec(wl5_full)
wl6_full = eltype(wl6_full) <: Number ? wl6_full : to_float_vec(wl6_full)
wl7_full = eltype(wl7_full) <: Number ? wl7_full : to_float_vec(wl7_full)
wl8_full = eltype(wl8_full) <: Number ? wl8_full : to_float_vec(wl8_full)
wl9_full = eltype(wl9_full) <: Number ? wl9_full : to_float_vec(wl9_full)
wl10_full = eltype(wl10_full) <: Number ? wl10_full : to_float_vec(wl10_full)
wl11_full = eltype(wl11_full) <: Number ? wl11_full : to_float_vec(wl11_full)
wl12_full = eltype(wl12_full) <: Number ? wl12_full : to_float_vec(wl12_full)
wl13_full = eltype(wl13_full) <: Number ? wl13_full : to_float_vec(wl13_full)
wl14_full = eltype(wl14_full) <: Number ? wl14_full : to_float_vec(wl14_full)
wl15_full = eltype(wl15_full) <: Number ? wl15_full : to_float_vec(wl15_full)
wl16_full = eltype(wl16_full) <: Number ? wl16_full : to_float_vec(wl16_full)       

v_n   = collect(align_to_window(v_full,  length(t_n)))
wl1_n = collect(align_to_window(wl1_full, length(t_n)))  # kN
wl2_n = collect(align_to_window(wl2_full, length(t_n)))  # kN
wl3_n = collect(align_to_window(wl3_full, length(t_n)))  # kN
wl4_n = collect(align_to_window(wl4_full, length(t_n)))  # kN
wl5_n = collect(align_to_window(wl5_full, length(t_n)))  # kN
wl6_n = collect(align_to_window(wl6_full, length(t_n)))  # kN
wl7_n = collect(align_to_window(wl7_full, length(t_n)))  # kN
wl8_n = collect(align_to_window(wl8_full, length(t_n)))  # kN
wl9_n = collect(align_to_window(wl9_full, length(t_n)))  # kN
wl10_n = collect(align_to_window(wl10_full, length(t_n)))# kN
wl11_n = collect(align_to_window(wl11_full, length(t_n)))#  kN
wl12_n = collect(align_to_window(wl12_full, length(t_n)))#  kN
wl13_n = collect(align_to_window(wl13_full, length(t_n)))#  kN
wl14_n = collect(align_to_window(wl14_full, length(t_n)))#  kN
wl15_n = collect(align_to_window(wl15_full, length(t_n)))#  kN
wl16_n = collect(align_to_window(wl16_full, length(t_n)))#  kN      


# ---------------- 9) summaries + probability of failure ----------------

# --- report using the *critical* sample (max displacement), but keep plotting nmc ---
crit_val, crit_idx = findmax(max_abs_disp)
println(">> Reporting maxima for critical sample = $crit_idx (max_abs_disp = $(crit_val) m)")

# align critical series to its own time window
t_c  = collect(times[crit_idx])
u_c  = collect(disp_node1[crit_idx])

v_full_c   = collect(wl_base[crit_idx])
wl2_full_c = collect(wl_node2[crit_idx])

# handle object-valued loads
wl2_full_c = eltype(wl2_full_c) <: Number ? wl2_full_c : to_float_vec(wl2_full_c)

v_c   = collect(align_to_window(v_full_c,   length(t_c)))
wl2_c = collect(align_to_window(wl2_full_c, length(t_c)))

println("Young's modulus (critical sample): $(E_vals[crit_idx]) Pa")
println("Maximum wind speed (critical sample): $(maximum(abs.(v_c))) m/s")
println("Maximum wind load at node 2 (critical sample): $(maximum(abs.(wl2_c))) kN")


# ---------------- probability of failure + fails table ----------------
fail = max_abs_disp .> capacity
pf   = mean(fail)
println("Probability of failure (capacity = $(capacity) m): $(round(pf, digits=4))  ($(sum(fail)) of $(N_MC) samples)")

wind_max = map(v -> maximum(v), wl_base)  # max wind speed per sample (m/s)
df_fail = DataFrame(
    sample       = findall(fail),
    max_abs_disp = max_abs_disp[fail],
    E            = E_vals[fail],
    wind_max_mps = wind_max[fail],
)
println(df_fail)
CSV.write(joinpath(results_dir, "fails.csv"), df_fail)


# Combined time series (wind, loads at nodes 1 & 2)
p_disp = plot(t_n, v_n; 
    label="Wind speed (m/s)", 
    color=:grey, linewidth=1,
    xlabel="Time (s)", ylabel="Wind speed / Wind load",
    legend=:outerright, legendfontsize=10,
    xguidefont=12, yguidefont=12,
    xtickfont=10, ytickfont=10,
    size=(1000,700),          # higher resolution
    dpi=300
)

# Wind loads 1–5 (red/pink/yellow)
plot!(p_disp, t_n, wl1_n;  label="Wind load node 1 (kN)",  linewidth=1, color=:red)
plot!(p_disp, t_n, wl2_n;  label="Wind load node 2 (kN)",  linewidth=1, color=:orange)
plot!(p_disp, t_n, wl3_n;  label="Wind load node 3 (kN)",  linewidth=1, color=:orange)
plot!(p_disp, t_n, wl4_n;  label="Wind load node 4 (kN)",  linewidth=1, color=:orange)
plot!(p_disp, t_n, wl5_n;  label="Wind load node 5 (kN)",  linewidth=1, color=:red)

# Wind loads 6–16 (blue/green tones)
plot!(p_disp, t_n, wl6_n;  label="Wind load node 6 (kN)",  linewidth=1, color=:red)
plot!(p_disp, t_n, wl7_n;  label="Wind load node 7 (kN)",  linewidth=1, color=:orange)
plot!(p_disp, t_n, wl8_n;  label="Wind load node 8 (kN)",  linewidth=1, color=:orange)
plot!(p_disp, t_n, wl9_n;  label="Wind load node 9 (kN)",  linewidth=1, color=:orange)
plot!(p_disp, t_n, wl10_n; label="Wind load node 10 (kN)", linewidth=1, color=:orange)
plot!(p_disp, t_n, wl11_n; label="Wind load node 11 (kN)", linewidth=1, color=:blue)
plot!(p_disp, t_n, wl12_n; label="Wind load node 12 (kN)", linewidth=1, color=:blue)
plot!(p_disp, t_n, wl13_n; label="Wind load node 13 (kN)", linewidth=1, color=:green)
plot!(p_disp, t_n, wl14_n; label="Wind load node 14 (kN)", linewidth=1, color=:green)
plot!(p_disp, t_n, wl15_n; label="Wind load node 15 (kN)", linewidth=1, color=:green)
plot!(p_disp, t_n, wl16_n; label="Wind load node 16 (kN)", linewidth=1, color=:green)

display(p_disp)
savefig(p_disp, joinpath(results_dir, "timeseries_sample$(nmc).png"))



# Plot: influence of the hight
p_disp = plot(t_n, v_n; 
    label="Wind speed (m/s)", 
    color=:grey, linewidth=1,
    xlabel="Time (s)", ylabel="Wind speed / Wind load",
    legend=:right, legendfontsize=10,
    xguidefont=12, yguidefont=12,
    xtickfont=10, ytickfont=10,
    dpi=300
)
plot!(p_disp, t_n, wl15_n; label="Wind load node 15 (kN)", linewidth=1, color=:blue)
plot!(p_disp, t_n, wl16_n; label="Wind load node 16 (kN)", linewidth=1, color=:red)

display(p_disp)
savefig(p_disp, joinpath(results_dir, "timeseries1_sample$(nmc).png"))

# Normalized (wind vs top displacement)
norm_wind = v_n ./ max(maximum(abs.(v_n)), eps())
norm_disp = u_n ./ max(maximum(abs.(u_n)), eps())
p_norm = plot(t_n, norm_disp; label="Normalized Displacement",
    xlabel="Time (s)", ylabel="Normalized value", legend=:topright, linewidth=2)
plot!(p_norm, t_n, norm_wind; label="Normalized Wind Speed", linewidth=2)
savefig(p_norm, joinpath(results_dir, "normalized_sample$(nmc).png"))

# All-samples wind speeds & loads (first up to 100 samples) for chosen load node
K = min(100, N_MC)
legendval = K > 10 ? false : :bottomright

wl_all_aligned  = Vector{Vector{Float64}}(undef, K)  # kN (node = LOAD_NODE_FOR_ALL)
wb_all_aligned  = Vector{Vector{Float64}}(undef, K)  # m/s

for i in 1:K
    ti  = samples.sim_time[i]
    vi   = collect(samples.wl_base[i])
    wli = collect(load_col[i])        # kN
    wli = eltype(wli) <: Number ? wli : to_float_vec(wli)

    wb_all_aligned[i] = collect(align_to_window(vi,   length(ti)))
    wl_all_aligned[i] = collect(align_to_window(wli, length(ti)))
end

same_len = length(unique(length.(samples.sim_time[1:K]))) == 1
if same_len
    x = collect(samples.sim_time[1])
    p_all_loads = plot(x, wl_all_aligned; label="Wind load node $(LOAD_NODE_FOR_ALL) (kN)",
        xlabel="Time (s)", ylabel="Wind load (kN)",
        title="Wind Loads (first $(K) samples)", legend=legendval)
    p_all_wind = plot(x, wb_all_aligned; label="Wind speed (m/s)",
        xlabel="Time (s)", ylabel="Wind speed (m/s)",
        title="Wind Speeds (first $(K) samples)", legend=legendval)
else
    p_all_loads = plot([eachindex(y) for y in wl_all_aligned], wl_all_aligned;
        label="Wind load node $(LOAD_NODE_FOR_ALL) (kN)", xlabel="Index",
        ylabel="Wind load (kN)", title="Wind Loads (first $(K) samples)", legend=legendval)
    p_all_wind = plot([eachindex(y) for y in wb_all_aligned], wb_all_aligned;
        label="Wind speed (m/s)", xlabel="Index",
        ylabel="Wind speed (m/s)", title="Wind Speeds (first $(K) samples)", legend=legendval)
end
savefig(p_all_wind,  joinpath(results_dir, "all_wind_speeds.png"))
savefig(p_all_loads, joinpath(results_dir, "all_wind_loads_node$(LOAD_NODE_FOR_ALL).png"))

# Histograms
p_hist_max = histogram(max_abs_disp;
    bins=30, xlabel="Maximum Displacement (m)",
    ylabel="Frequency", title="Histogram of Maximum Absolute Displacement",
    legend=false)

# Concatenate all realized loads at nodes 1 & 2 (kN)
wl1_all_kN = vcat([ (eltype(wl_node1[i]) <: Number ? wl_node1[i] : to_float_vec(wl_node1[i])) for i in 1:N_MC ]...)
wl2_all_kN = vcat([ (eltype(wl_node2[i]) <: Number ? wl_node2[i] : to_float_vec(wl_node2[i])) for i in 1:N_MC ]...)
p_hist_wl = histogram(wl1_all_kN; bins=30,
    xlabel="Wind Load (kN)", ylabel="Frequency",
    label="Node 1", title="Histogram of Wind Loads", legend=true)
histogram!(p_hist_wl, wl2_all_kN; bins=30, label="Node 2")

p_hist = plot(p_hist_max, p_hist_wl; layout=(2,1), size=(900,700), title="Histograms")
savefig(p_hist, joinpath(results_dir, "histograms.png"))

println("Fails table saved to: ", joinpath(results_dir, "fails.csv"))
println("Plots saved to: ", results_dir)
