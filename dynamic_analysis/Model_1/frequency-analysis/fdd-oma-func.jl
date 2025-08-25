using DSP, LinearAlgebra
using Peaks # For findmaxima (Pkg.add("Peaks") if you don't have it)
# using FFTW # DSP.jl usually re-exports fftfreq

# --- FDD Core Function ---
function perform_fdd(signals::Matrix{Float64}, fs::Real;
                     nperseg::Int=1024, noverlap_ratio::Float64=0.5,
                     window_func=hann)
    num_samples, num_channels = size(signals)
    noverlap = floor(Int, nperseg * noverlap_ratio)

    stfts_coefficient_matrices = [] # Store the actual coefficient matrices
    local freqs::Vector{Float64} # Ensure type stability
    
    if num_channels == 0
        error("No signals provided.")
    end

    println("Calculating STFTs...")
    for i = 1:num_channels
        # Perform STFT
        # S_output is what stft returns. We will check its type.
        S_output = stft(signals[:, i], nperseg, noverlap; fs=fs, window=window_func(nperseg))
        
        # Store the STFT coefficients (which should be a matrix)
        # If S_output is an STFT object, S_output[:,:] extracts the matrix.
        # If S_output is ALREADY a matrix, this just makes a copy.
        current_stft_coeffs = S_output[:,:] # Make sure this is a Matrix{ComplexF64}
        push!(stfts_coefficient_matrices, current_stft_coeffs)

        if i == 1
            println("Debug: Type of S_output from stft for channel 1: ", typeof(S_output))
            println("Debug: Size of current_stft_coeffs for channel 1: ", size(current_stft_coeffs))

            # --- MODIFICATION: Manually create frequency vector ---
            # The number of frequency bins in the output of stft (one-sided) is nfft/2 + 1.
            # DSP.stft uses nfft = nperseg by default if not specified otherwise.
            nfft_actual = nperseg # Assuming stft uses nperseg as nfft
            num_freq_bins_calculated = nfft_actual ÷ 2 + 1
            
            # Verify this matches the number of rows from the STFT output
            if size(current_stft_coeffs, 1) != num_freq_bins_calculated
                println("Warning: Mismatch in expected number of frequency bins.")
                # Fallback or adjust nfft_actual if stft pads differently
                num_freq_bins_calculated = size(current_stft_coeffs, 1)
            end

            # Use DSP.rfftfreq for robust frequency vector generation
            # fftfreq needs nfft (the length of the FFT window)
            freqs = DSP.rfftfreq(nfft_actual, fs) # This gives the correct one-sided frequencies
            # Ensure freqs matches the number of rows from STFT output
            if length(freqs) != size(current_stft_coeffs,1)
                println("Warning: Length of generated freqs does not match STFT output rows. Adjusting.")
                # This might happen if stft uses a different nfft internally, e.g. nextpow2(nperseg)
                # For now, let's assume DSP.rfftfreq with nperseg is correct or stft output rows define it
                # If DSP.stft output S_output was an STFT object, DSP.freq(S_output) would be best.
                # As a robust fallback if it's a matrix:
                freqs = LinRange(0.0, fs/2.0, size(current_stft_coeffs, 1))

            end
            println("Debug: Generated frequency vector length: ", length(freqs))
            if !isempty(freqs)
                 println("Debug: Freq vector starts $(freqs[1]), ends $(freqs[end])")
            end
        end
    end
    
    if isempty(stfts_coefficient_matrices)
        error("STFT processing failed to produce coefficient matrices.")
    end

    num_freq_bins = length(freqs)
    num_segments = size(stfts_coefficient_matrices[1], 2)

    max_sv_to_store = num_channels
    singular_values_over_freq = zeros(Float64, max_sv_to_store, num_freq_bins)
    left_singular_vectors_tensor = zeros(ComplexF64, num_channels, max_sv_to_store, num_freq_bins)

    println("Performing SVD for each frequency bin...")
    for f_idx = 1:num_freq_bins
        Y_fk = zeros(ComplexF64, num_channels, num_segments)
        for chan_idx = 1:num_channels
            Y_fk[chan_idx, :] = stfts_coefficient_matrices[chan_idx][f_idx, :]
        end

        G_yy_fk = (1.0 / num_segments) * (Y_fk * Y_fk')
        F_svd = svd(G_yy_fk)
        
        num_available_sv = length(F_svd.S)
        num_to_copy = min(max_sv_to_store, num_available_sv)

        if num_to_copy > 0 # Ensure there are singular values to copy
            singular_values_over_freq[1:num_to_copy, f_idx] = F_svd.S[1:num_to_copy]
            left_singular_vectors_tensor[:, 1:num_to_copy, f_idx] = F_svd.U[:, 1:num_to_copy]
        elseif num_available_sv == 0 && f_idx == 1 # Warning only for first freq bin if SVD is empty
            println("Warning: SVD resulted in no singular values at freq_idx $f_idx.")
        end
    end
    println("SVD computation complete.")
    return freqs, singular_values_over_freq, left_singular_vectors_tensor
end

# --- Simple Peak Picking Function (find_fdd_peaks) ---
# ... (find_fdd_peaks function remains the same as before) ...
function find_fdd_peaks(freqs::Vector{Float64}, singular_values_matrix::Matrix{Float64};
                        prominence_threshold_factor::Float64 = 0.1,
                        max_peaks::Int = 5)
    if isempty(freqs) || isempty(singular_values_matrix)
        println("Warning: Empty frequencies or singular values for peak picking.")
        return Float64[], Int[]
    end
    first_singular_values = singular_values_matrix[1, :]
    if isempty(first_singular_values)
        println("Warning: First singular value vector is empty.")
        return Float64[], Int[]
    end
    
    # Use an absolute threshold as well to avoid picking noise if max SV is tiny
    min_abs_threshold = 1e-8 * maximum(first_singular_values) # Heuristic
    min_prominence_val = prominence_threshold_factor * maximum(first_singular_values)

    peak_indices_all, props = findmaxima(first_singular_values) 
    
    significant_peaks = Tuple{Int, Float64}[] # Store (index, value)
    for idx in peak_indices_all
        # A more robust peak selection might involve checking actual prominence values
        # from a more advanced findmaxima if available, or comparing to neighbors.
        # Simple check: value > threshold_factor * global_max
        if first_singular_values[idx] > min_prominence_val && first_singular_values[idx] > min_abs_threshold
            # Further check if it's a clear "peak" relative to its local neighborhood (simplified)
            is_prominent_locally = true
            local_neighborhood = max(1, idx-5):min(length(first_singular_values), idx+5)
            if idx > 1 && idx < length(first_singular_values)
                 # Example: Check if it's significantly larger than its immediate neighbors
                 if first_singular_values[idx] < 1.1 * first_singular_values[idx-1] && first_singular_values[idx] < 1.1 * first_singular_values[idx+1]
                     # Potentially part of a flat top or minor ripple, could skip
                 end
            end
            if is_prominent_locally
                push!(significant_peaks, (idx, first_singular_values[idx]))
            end
        end
    end

    if isempty(significant_peaks)
        println("No significant peaks found with current criteria.")
        return Float64[], Int[]
    end
    
    sort!(significant_peaks, by = x -> x[2], rev=true)
    
    selected_peak_indices = [p[1] for p in significant_peaks]
    if length(selected_peak_indices) > max_peaks
        selected_peak_indices = selected_peak_indices[1:max_peaks]
    end
    
    sort!(selected_peak_indices) # Sort by frequency again

    peak_frequencies = freqs[selected_peak_indices]
    
    return peak_frequencies, selected_peak_indices
end

function estimated_mode_shapes(peak_indices, identified_frequencies, U_vectors_tensor)
    estimated_mode_shapes = []
    for (i, p_idx) in enumerate(peak_indices)
        freq = identified_frequencies[i]
        mode_shape_complex = U_vectors_tensor[:, 1, p_idx] 
        max_abs_val, max_idx = findmax(abs.(mode_shape_complex))
        if max_abs_val > 1e-9 # Avoid division by zero if mode shape is null
            mode_shape_normalized = mode_shape_complex ./ mode_shape_complex[max_idx]
        else
            mode_shape_normalized = mode_shape_complex # Keep as is if null
        end
        push!(estimated_mode_shapes, mode_shape_normalized)
        println("Mode $(i+1) at $(round(freq, digits=3)) Hz:")
        for (ch, val) in enumerate(mode_shape_normalized)
            mag = abs(val)
            phase_deg = rad2deg(angle(val))
            println("  Channel $ch: Magnitude = $(round(mag, digits=3)), Phase = $(round(phase_deg, digits=3)) deg")
        end
    end
    return estimated_mode_shapes
end

function plot_fdd(freq_vector, S_values_matrix, identified_frequencies, peak_indices,
                  estimated_mode_shapes, num_channels)
    # --- 5. Plotting ---
    plot_fdd = plot(freq_vector, S_values_matrix[1, :], label="1st Singular Value", yaxis=:log,
                        xlabel="Frequency (Hz)", ylabel="Singular Value (log scale)", legend=:topright)
    if size(S_values_matrix,1) >= 2
        plot!(plot_fdd, freq_vector, S_values_matrix[2, :], label="2nd Singular Value", linestyle=:dash)
    end
    if size(S_values_matrix,1) >= 3 && size(S_values_matrix,2) == length(freq_vector)
        plot!(plot_fdd, freq_vector, S_values_matrix[3, :], label="3rd Singular Value", linestyle=:dot)
    end
    scatter!(plot_fdd, identified_frequencies, S_values_matrix[1, peak_indices], color=:red, label="Identified Peaks", markersize=5)
    title!("FDD Singular Values")
    display(plot_fdd)

    for i = 1:length(identified_frequencies)
        freq_hz = round(identified_frequencies[i], digits=2)
        mode_shape = estimated_mode_shapes[i]
        p_mag = bar(1:num_channels, abs.(mode_shape), title="Mode @ $(freq_hz) Hz (Mag)", legend=false, xlabel="Channel", ylabel="Magnitude", xticks=1:num_channels)
        p_phase = bar(1:num_channels, rad2deg.(angle.(mode_shape)), title="Mode @ $(freq_hz) Hz (Phase)", legend=false, xlabel="Channel", ylabel="Phase (deg)", xticks=1:num_channels)
        plot_mode = plot(p_mag, p_phase, layout=(1,2))
        display(plot_mode)
    end
end