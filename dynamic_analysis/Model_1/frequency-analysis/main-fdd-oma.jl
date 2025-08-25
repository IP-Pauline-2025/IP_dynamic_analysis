using Distributed, Printf, Dates, Plots, CSV, Plots
time_stamp = Dates.format(now(),  "yyyymmdd_HHMMSS")
@printf "[%s]: " time_stamp
println("Reliability Analysis of a Dynamic Tower-like Structure under Wind loads")
tic_total = time()

addprocs(6)

@everywhere begin
    using UncertaintyQuantification, DelimitedFiles
    hostname = gethostname()
    if hostname == "IRZ-PC1337"
        include("C:/OMNISSIAH/Julia/2025-08-pauline-git/dynamic_analysis/calculate_wind_force.jl")
    else
        include("C:/IP/dynamic_analysis/calculate_wind_force.jl")
    end

    # parameter definition, they can be changed but are the same for each sample
    Δt = Parameter(0.02, :dt)                       # time step size
    T = Parameter(120.0, :T)                       # total simulation time
    E = RandomVariable(Uniform(2.0e8, 2.2e8), :E)   # Young's modulus
    Iz= RandomVariable(Uniform(2e-2, 1e-1), :Iz)    # second moment of area
    t = collect(0:Δt.value:T.value)                      # time vector
    timeSteps = Parameter(length(t), :timeSteps)    # number of time steps

    # contstant PSD function
    fdisc = collect(0.01:0.01:100)
    σ_u = 1 # contant value
    constantpsdfunc = (f) -> σ_u
    psdvalues = constantpsdfunc.(fdisc)

    #coinstruct the empirical PSD representing the von Karman spectrum
    constantPSD = EmpiricalPSD(fdisc, psdvalues)

    wl = SpectralRepresentation(constantPSD, t, :wl)
    wl_model = StochasticProcessModel(wl)

    # directory where the OpenSees input files are stored

    if hostname == "IRZ-PC1337"
        sourcedir = joinpath(pwd(), "C:/OMNISSIAH/Julia/2025-08-pauline-git/dynamic_analysis/Model_1/frequency-analysis/")
    else
        sourcedir = joinpath(pwd(), "C:/IP/dynamic_analysis/Model_1")
    end
    sourcefile = ["FEM_1_eigen_uq.tcl", "wind-speed.dat"]

    numberformats = Dict(:dt => ".8e", :T => ".8e", :wl => ".8e", :E => ".8e", :Iz => ".8e")

    # case where the results are stored
    workdir = joinpath(pwd(), "workdir-1")

    max_abs_disp = Extractor(base -> begin
        file = joinpath(base, "displacement3.out")
        data = DelimitedFiles.readdlm(file, ' ')

        return maximum(abs.(data[:, 2]))
    end, :max_abs_disp)

    # this is the displacement time history, only needed for plotting
    disp2 = Extractor(base -> begin
        file = joinpath(base, "displacement2.out")
        data = DelimitedFiles.readdlm(file, ' ')

        return data[:, 2]
    end, :disp2)

    # this is the displacement time history, only needed for plotting
    disp3 = Extractor(base -> begin
        file = joinpath(base, "displacement3.out")
        data = DelimitedFiles.readdlm(file, ' ')

        return data[:, 2]
    end, :disp3)

    # this is the simulation time, also only needed for plotting
    sim_time = Extractor(base -> begin
        file = joinpath(base, "displacement2.out")
        data = DelimitedFiles.readdlm(file, ' ')

        return data[:, 1]
    end, :sim_time)


    opensees = Solver(
        "OpenSees", # path to OpenSees binary
        "FEM_1_eigen_uq.tcl";
        args="", # (optional) extra arguments passed to the solver
    )

        
    # define the external model that runs OpenSees, this is the model that will be run for each sample
    # it uses the sourcedir and sourcefile to find the OpenSees input files, and the workdir to store the results
    ext = ExternalModel(
        sourcedir, sourcefile, [max_abs_disp, disp2, disp3, sim_time], opensees; workdir=workdir, formats=numberformats
    )

    # this needs to include the models that are needed and that you want to use
    models = [wl_model, ext]
end

# this runs the reliability analysis, in this case a Monte Carlo simulation with N_MC samples
N_MC = 100 # Number of Monte Carlo Samples
println("Running Monte Carlo simulation with $N_MC samples...")
# this part: df -> 200 .- df.max_abs_disp
# actually defines the performance function also known as limit state function which is evaluated for each of the samples, if you use the same record this of course does not make any sense
pf, mc_std, samples = probability_of_failure(models, df -> 1 .- df.max_abs_disp, [Δt, T, timeSteps, wl, E, Iz], MonteCarlo(N_MC))
println
println("Probability of failure: $pf")

rmprocs()

#pdisp = plot(t, samples.disp[:]; label="disp in m", xlabel="time in s", ylabel="displacement")

#OpenSees Eigenvalue Analysis
#Mode 1: Eigenvalue = 37.57503864804687054857, Frequency = 0.9755956374475411 Hz
#Mode 2: Eigenvalue = 960.07288597010165176471, Frequency = 4.931422745755171 Hz

# --- Create the Data Matrix (samples x channels) ---
# Each column is a channel
my_signals = hcat(samples.disp2[:]...)
tt_loaded = t

duration = tt_loaded[end] - tt_loaded[1]
my_fs = round(1/(tt_loaded[2]-tt_loaded[1]))
#floaded = dfres_loaded.ω[1] # system eigenfrequencies

println("Generated test signal matrix with dimensions: ", size(my_signals))
println("Sampling frequency: ", my_fs, " Hz")
println("Signal duration: ", duration, " s")
#println("True Natural Frequencies (Hz): $floaded")

# Main execution
println("--- OMA FDD Example ---")

num_channels = size(my_signals,1)
num_samples = size(my_signals,2)
println("Generated data: $num_samples samples, $num_channels channels.")

nperseg_fdd = 2024
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

####################################################################################################################################
toc_total = time() - tic_total
println("Total elapsed time: $(round(toc_total, digits=3))s")
@printf "[%s]: " Dates.format(now(), "dd/mm/yy, HH:MM:SS")
println("Script finished.")