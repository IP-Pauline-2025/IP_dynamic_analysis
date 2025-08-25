#load two column time history data from a file in a DataFrame
using DataFrames, CSV, Plots
filename = "displacement.out"
dfres_loaded = CSV.read(filename, DataFrame, delim=' ', header=false)

# --- Create the Data Matrix (samples x channels) ---
# Each column is a channel
my_signals = dfres_loaded.Column2
tt_loaded = dfres_loaded.Column1
# Convert to a matrix with one channel
my_signals = hcat(my_signals)
plot(tt_loaded, my_signals, label="Displacement Signal", xlabel="Time (s)", ylabel="Displacement", title="Loaded Displacement Signal")
duration = tt_loaded[end] - tt_loaded[1]
my_fs = round(1/(tt_loaded[2]-tt_loaded[1]))
#floaded = dfres_loaded.ω[1] # system eigenfrequencies

println("Generated test signal matrix with dimensions: ", size(my_signals))
println("Sampling frequency: ", my_fs, " Hz")
println("Signal duration: ", duration, " s")
#println("True Natural Frequencies (Hz): $floaded")

# Main execution
println("--- OMA FDD Example ---")

num_channels = 1
num_samples = size(my_signals)
println("Generated data: $num_samples samples, $num_channels channels.")

nperseg_fdd = 1024
noverlap_ratio_fdd = 0.75

include("fdd-oma-func.jl")
freq_vector, S_values_matrix, U_vectors_tensor = perform_fdd(
    my_signals, my_fs,
    nperseg=nperseg_fdd,
    noverlap_ratio=noverlap_ratio_fdd
)

if !isempty(freq_vector)
    identified_frequencies, peak_indices = find_fdd_peaks(
        freq_vector, S_values_matrix,
        prominence_threshold_factor=0.05,
        max_peaks=num_channels * 2
    )

    println("\n--- FDD Results ---")
    #println("True Natural Frequencies (Hz): $(dfres_loaded.ω[1]) Hz")
    println("Identified Natural Frequencies (Hz): $identified_frequencies")

    estimated_mode_shapes_picked = []
    for (i, p_idx) in enumerate(peak_indices)
        freq = identified_frequencies[i]
        mode_shape_complex = U_vectors_tensor[:, 1, p_idx] 
        max_abs_val, max_idx = findmax(abs.(mode_shape_complex))
        if max_abs_val > 1e-9 # Avoid division by zero if mode shape is null
            mode_shape_normalized = mode_shape_complex ./ mode_shape_complex[max_idx]
        else
            mode_shape_normalized = mode_shape_complex # Keep as is if null
        end
        push!(estimated_mode_shapes_picked, mode_shape_normalized)
        println("Mode $(i+1) at $(round(freq, digits=3)) Hz:")
        for (ch, val) in enumerate(mode_shape_normalized)
            mag = abs(val)
            phase_deg = rad2deg(angle(val))
            println("  Channel $ch: Magnitude = $(round(mag, digits=3)), Phase = $(round(phase_deg, digits=3)) deg")
        end
    end

    # --- 5. Plotting ---
    plot_fdd_1 = plot(freq_vector, S_values_matrix[1, :], label="1st Singular Value", yaxis=:log,
                      xlabel="Frequency (Hz)", ylabel="Singular Value (log scale)", legend=:topright)
    if size(S_values_matrix,1) >= 2
        plot!(plot_fdd_1, freq_vector, S_values_matrix[2, :], label="2nd Singular Value", linestyle=:dash)
    end
    if size(S_values_matrix,1) >= 3 && size(S_values_matrix,2) == length(freq_vector)
        plot!(plot_fdd_1, freq_vector, S_values_matrix[3, :], label="3rd Singular Value", linestyle=:dot)
    end
    scatter!(plot_fdd_1, identified_frequencies, S_values_matrix[1, peak_indices], color=:red, label="Identified Peaks", markersize=5)
    title!("FDD Singular Values")
    display(plot_fdd_1)

    for i = 1:length(identified_frequencies)
        freq_hz = round(identified_frequencies[i], digits=2)
        mode_shape = estimated_mode_shapes_picked[i]
        p_mag = bar(1:num_channels, abs.(mode_shape), title="Mode @ $(freq_hz) Hz (Mag)", legend=false, xlabel="Channel", ylabel="Magnitude", xticks=1:num_channels)
        p_phase = bar(1:num_channels, rad2deg.(angle.(mode_shape)), title="Mode @ $(freq_hz) Hz (Phase)", legend=false, xlabel="Channel", ylabel="Phase (deg)", xticks=1:num_channels)
        plot_mode = plot(p_mag, p_phase, layout=(1,2))
        display(plot_mode)
    end
else
    println("FDD did not produce a valid frequency vector. Skipping peak picking and plotting.")
end

println("\nNote: Mode shape phases are relative. One channel's phase can be set to zero for reference.")