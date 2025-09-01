############################################################
# Plotting from saved `samples` DataFrame (Model_1)
############################################################
using JLD2, DataFrames, JSON3, Plots, Statistics, CSV

# ---------------- paths ----------------
results_dir = raw"C:\IP\results\MCS\Model_1\plots"
jld_path    = raw"C:\IP\results\MCS\Model_1\20250830_130353_model_1_samples.jld2"
meta_path   = raw"C:\IP\results\MCS\Model_1\20250830_130353_model_1_metadata.json"

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

# check dataframe columns (Model_1)
for name in [:sim_time, :disp, :wl_base, :wl_node2, :wl_node1, :max_abs_disp, :E]
    @assert hasproperty(samples, name) "Missing column: $name. Available: $(names(samples))"
end

# Safe numeric extraction for wind-load items
getF(x) = x isa Number ? float(x) :
          hasproperty(x, :F_wp) ? float(getproperty(x, :F_wp)) :
          throw(ArgumentError("Cannot extract numeric value; missing :F_wp on $(typeof(x))."))
to_float_vec(v) = map(getF, v)

# Align a full-length signal to the modified time window (t >= t_start)
function align_to_window(full_vec::AbstractVector, window_len::Int)
    offset = length(full_vec) - window_len
    if offset >= 0
        @views return full_vec[(offset+1):end]        # drop leading samples
    elseif !isempty(full_vec)
        return vcat(fill(first(full_vec), -offset), full_vec)  # pad if shorter (rare)
    else
        return fill(0.0, window_len)                  # empty -> zeros
    end
end

# ---------------- columns & setup ----------------
times    = samples.sim_time
disp     = samples.disp
wl_base  = samples.wl_base
wl_node2 = samples.wl_node2
wl_node1 = samples.wl_node1
E_vals       = samples.E
max_abs_disp = samples.max_abs_disp
N_MC = length(times)

# capacity for failure check (m)
capacity = 0.008


# ---------------- choose sample & align ----------------
nmc = min(544, N_MC)
t_n  = collect(times[nmc])     # already t >= 20 s
u_n  = collect(disp[nmc])

v_full   = collect(wl_base[nmc])
wl2_full = collect(wl_node2[nmc])   # kN
wl1_full = collect(wl_node1[nmc])   # kN

wl2_full = eltype(wl2_full) <: Number ? wl2_full : to_float_vec(wl2_full)
wl1_full = eltype(wl1_full) <: Number ? wl1_full : to_float_vec(wl1_full)

v_n   = collect(align_to_window(v_full,  length(t_n)))
wl2_n = collect(align_to_window(wl2_full, length(t_n)))  # kN
wl1_n = collect(align_to_window(wl1_full, length(t_n)))  # kN

# ---------------- 9) summaries + probability of failure ----------------
println("Maximum tip displacement (overall): $(maximum(max_abs_disp)) m")
println("Young's modulus (sample $nmc): $(E_vals[nmc]) Pa")
println("Maximum wind speed (sample $nmc): $(maximum(abs.(v_n))) m/s")
println("Maximum wind load at node 2 (sample $nmc): $(maximum(abs.(wl2_n))) kN")

# Probability of failure based on capacity
fail = samples.max_abs_disp .> capacity
pf   = mean(fail)  # fraction of failing samples
println("Probability of failure (capacity = $(capacity) m): $(round(pf, digits=4))  ($(sum(fail)) of $(N_MC) samples)")

# Failure table + CSV
wind_max = map(v -> maximum(v), samples.wl_base)
df_fail = DataFrame(
    sample       = findall(fail),
    max_abs_disp = samples.max_abs_disp[fail],
    E            = samples.E[fail],
    wind_max_mps = wind_max[fail],
)
println(df_fail)
CSV.write(joinpath(results_dir, "fails.csv"), df_fail)

# ---------------- plots ----------------
# Combined time series

# ---------------- plots ----------------
# Combined time series (without displacement)
p_wind = plot(t_n, v_n; 
    label="Wind speed (m/s)",
    xlabel="Time (s)", 
    ylabel="Wind speed / Wind load",
    legend=:right,            
    legendfontsize=12,        # bigger legend
    xguidefont=12, yguidefont=12,   # bigger axis labels
    xtickfont=12, ytickfont=12,     # bigger tick numbers        
    dpi=300             # better quality for saved figure
)
plot!(p_wind, t_n, wl2_n; 
    label="Wind load node 2 (kN)", 
    linewidth=2)
plot!(p_wind, t_n, wl1_n; 
    label="Wind load node 1 (kN)", 
    linewidth=2)

display(p_wind)
savefig(p_wind, joinpath(results_dir, "timeseries_sample$(nmc).png"))


# Normalized (wind vs displacement)
norm_wind = v_n ./ max(maximum(abs.(v_n)), eps())
norm_disp = u_n ./ max(maximum(abs.(u_n)), eps())
p_norm = plot(t_n, norm_disp; label="Normalized Displacement",
    xlabel="Time (s)", ylabel="Normalized value", legend=:topright, linewidth=2)
plot!(p_norm, t_n, norm_wind; label="Normalized Wind Speed", linewidth=2)
display(p_norm)
savefig(p_norm, joinpath(results_dir, "normalized_sample$(nmc).png"))

# All-samples wind speeds & node-2 loads (first up to 100 samples)
K = min(100, N_MC)
legendval = K > 10 ? false : :bottomright

wl2_all_aligned = Vector{Vector{Float64}}(undef, K)  # kN
wl_base_aligned = Vector{Vector{Float64}}(undef, K)

for i in 1:K
    ti  = samples.sim_time[i]
    vi   = collect(samples.wl_base[i])
    wl2 = collect(samples.wl_node2[i])               # kN
    wl2 = eltype(wl2) <: Number ? wl2 : to_float_vec(wl2)

    wl_base_aligned[i] = collect(align_to_window(vi,   length(ti)))
    wl2_all_aligned[i] = collect(align_to_window(wl2, length(ti)))
end

same_len = length(unique(length.(samples.sim_time[1:K]))) == 1
if same_len
    x = collect(samples.sim_time[1])
    p_wl2 = plot(x, wl2_all_aligned; label="Wind load node 2 (kN)",
        xlabel="Time (s)", ylabel="Wind load (kN)",
        title="Wind Loads at Node 2 (first $(K) samples)", legend=legendval)
    p_wlb = plot(x, wl_base_aligned; label="Wind speed (m/s)",
        xlabel="Time (s)", ylabel="Wind speed (m/s)",
        title="Wind Speeds (first $(K) samples)", legend=legendval)
else
    p_wl2 = plot([eachindex(y) for y in wl2_all_aligned], wl2_all_aligned;
        label="Wind load node 2 (kN)", xlabel="Index",
        ylabel="Wind load (kN)", title="Wind Loads at Node 2 (first $(K) samples)", legend=legendval)
    p_wlb = plot([eachindex(y) for y in wl_base_aligned], wl_base_aligned;
        label="Wind speed (m/s)", xlabel="Index",
        ylabel="Wind speed (m/s)", title="Wind Speeds (first $(K) samples)", legend=legendval)
end

display(p_wlb); display(p_wl2)
savefig(p_wlb, joinpath(results_dir, "all_wind_speeds.png"))
savefig(p_wl2, joinpath(results_dir, "all_wind_loads_node2.png"))


# Histograms
p_hist_max = histogram(max_abs_disp;
    bins=30, xlabel="Maximum Displacement (m)",
    ylabel="Frequency", title="Histogram of Maximum Absolute Displacement",
    legend=false)

wl2_all_kN = vcat([ (eltype(samples.wl_node2[i]) <: Number ?
                        samples.wl_node2[i] : to_float_vec(samples.wl_node2[i]))
                    for i in 1:N_MC ]...)
wl1_all_kN = vcat([ (eltype(samples.wl_node1[i]) <: Number ?
                        samples.wl_node1[i] : to_float_vec(samples.wl_node1[i]))
                    for i in 1:N_MC ]...)

p_hist_wl = histogram(wl2_all_kN; bins=30,
    xlabel="Wind Load (kN)", ylabel="Frequency",
    label="Node 2", title="Histogram of Wind Loads", legend=true)
histogram!(p_hist_wl, wl1_all_kN; bins=30, label="Node 1")

p_hist = plot(p_hist_max, p_hist_wl; layout=(2,1), size=(900,700), title="Histograms")
display(p_hist)
savefig(p_hist, joinpath(results_dir, "histograms.png"))

println("Fails table saved to: ", joinpath(results_dir, "fails.csv"))
println("Plots saved to: ", results_dir)
