using Distributed, Printf, Dates, Plots, JSON, JLD2
time_stamp = Dates.format(now(),  "yyyymmdd_HHMMSS")
@printf "[%s]: " time_stamp
println("Reliability Analysis of a Dynamic Tower-like Structure under Wind loads")
tic_total = time()

# Automatically detect number of workers from Slurm
jobid = get(ENV, "SLURM_JOB_ID", "local")
if jobid == "local"
    nworkerz = 6
    addprocs(nworkerz)
    println("Number of workers: $(nworkerz)")
else
    nprocs_requested = nprocs_requested = parse(Int, get(ENV, "SLURM_NTASKS", "1"))
    if nprocs_requested > nworkers()
        addprocs(nprocs_requested - (nworkers()-1))
    end
    println("Number of workers: $(nworkers())")
end

DOWEWANTplots = true
DOWEWANTseeplots = false
DOWEWANTsaveplots = true
DOWEWANTsaveresults = true

@everywhere begin
    using UncertaintyQuantification, DelimitedFiles
    hostname = gethostname()
    if hostname == "IRZ-PC1337"
        include("C:/OMNISSIAH/Julia/2025-08-pauline-git/dynamic_analysis/calculate_wind_force.jl")
    elseif Sys.islinux()
        include("/home/bittner/Documents/Julia/2025-08-22-Wind-Tower-Truss-Structure/dynamic_analysis/calculate_wind_force.jl")
    else
        include("C:/IP/dynamic_analysis/calculate_wind_force.jl")
    end

    # parameter definition, they can be changed but are the same for each sample
    Δt = Parameter(0.02, :dt)                       # time step size
    E = RandomVariable(Uniform(2.0e8, 2.2e8), :E)   # Young's modulus
    T = Parameter(100, :T)                          # total time of simulation
    t = collect(0:Δt.value:T.value)                 # time vector
    timeSteps = Parameter(length(t), :timeSteps)    # number of time steps

    # von Karman spectrum
    fdisc = collect(0.01:0.01:50)
    # roughly derived parameters from https://www.sciencedirect.com/science/article/pii/S0167610520302725
    σ_u = 1 # standard deviation of wind speed in m/s
    L_u = 50.8 # length scale in m
    v_mean = 30.0 # mean wind speed in m/s
    karmanfunc = (f) -> σ_u^2 * (4*L_u / v_mean) / ((1 + (2*π*f * L_u / v_mean)^2)^(5/6))
    psdvalues = karmanfunc.(fdisc)

    #coinstruct the empirical PSD representing the von Karman spectrum
    epWindPSD = EmpiricalPSD(fdisc, psdvalues)

    wl = SpectralRepresentation(epWindPSD, t, :wl)
    wl_model = StochasticProcessModel(wl)

    # directory where the OpenSees input files are stored

    if hostname == "IRZ-PC1337"
        sourcedir = joinpath(pwd(), "C:/OMNISSIAH/Julia/2025-08-pauline-git/dynamic_analysis/Model_2")
    elseif Sys.islinux()
        sourcedir = "/home/bittner/Documents/Julia/2025-08-22-Wind-Tower-Truss-Structure/dynamic_analysis/Model_2"
    else
        sourcedir = joinpath(pwd(), "C:/IP/dynamic_analysis/Model_2")
    end
    sourcefile = [
    "FEM_2.tcl",
    "wind-speed.dat",
    "wind-load-abs.dat",
    "wind-load-node1.dat",
    "wind-load-node2.dat",
    "wind-load-node3.dat",
    "wind-load-node4.dat",
    "wind-load-node5.dat",
    "wind-load-node6.dat",
    "wind-load-node7.dat",
    "wind-load-node8.dat",
    "wind-load-node9.dat",
    "wind-load-node10.dat",
    "wind-load-node11.dat",
    "wind-load-node12.dat",
    "wind-load-node13.dat",
    "wind-load-node14.dat",
    "wind-load-node15.dat",
    "wind-load-node16.dat"]
    
    numberformats = Dict(:dt => ".8e", :wl => ".8e", :E => ".8e", :T => ".8e")

    # case where the results are stored
    workdir = joinpath(pwd(), "workdir-2")

    t_start = 20.0  # nur Daten ab dieser Simulationszeit

    # define the models that extract the results from the OpenSees output files, we are interested in the maximum absolute displacement, the displacement at the top node and the simulation time
    max_abs_disp = Extractor(base -> begin
        file = joinpath(base, "displacement_node1.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return maximum(abs.(data[sel, 2]))
    end, :max_abs_disp)

    # this is the displacement time history, only needed for plotting

    disp_node1 = Extractor(base -> begin
        file = joinpath(base, "displacement_node1.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node1)

    disp_node2 = Extractor(base -> begin
        file = joinpath(base, "displacement_node2.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node2)

    disp_node3 = Extractor(base -> begin
        file = joinpath(base, "displacement_node3.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node3)

    disp_node4 = Extractor(base -> begin
        file = joinpath(base, "displacement_node4.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node4)

    disp_node5 = Extractor(base -> begin
        file = joinpath(base, "displacement_node5.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node5)

    disp_node6 = Extractor(base -> begin
        file = joinpath(base, "displacement_node6.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node6)

    disp_node7 = Extractor(base -> begin
        file = joinpath(base, "displacement_node7.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node7)

    disp_node8 = Extractor(base -> begin
        file = joinpath(base, "displacement_node8.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node8)

    disp_node9 = Extractor(base -> begin
        file = joinpath(base, "displacement_node9.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node9)

    disp_node10 = Extractor(base -> begin
        file = joinpath(base, "displacement_node10.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node10)

    disp_node11 = Extractor(base -> begin
        file = joinpath(base, "displacement_node11.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node11)

    disp_node12 = Extractor(base -> begin
        file = joinpath(base, "displacement_node12.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node12)

    disp_node13 = Extractor(base -> begin
        file = joinpath(base, "displacement_node13.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node13)

    disp_node14 = Extractor(base -> begin
        file = joinpath(base, "displacement_node14.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node14)

    disp_node15 = Extractor(base -> begin
        file = joinpath(base, "displacement_node15.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node15)

    disp_node16 = Extractor(base -> begin
        file = joinpath(base, "displacement_node16.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node16)

    disp_node17 = Extractor(base -> begin
        file = joinpath(base, "displacement_node17.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 2]
    end, :disp_node17)

    # this is the simulation time, also only needed for plotting
    sim_time = Extractor(base -> begin
        file = joinpath(base, "displacement_node1.out")
        data = DelimitedFiles.readdlm(file, ' ')
        sel  = data[:,1] .>= t_start
        return data[sel, 1]
    end, :sim_time)

    opensees = Solver(
        "OpenSees", # path to OpenSees binary
        "FEM_2.tcl";
        args="", # (optional) extra arguments passed to the solver
    )

    # alters the model of the wind speed to habe a baseline wind speed of v_mean
    base_wl_model = Model(df -> [v.+v_mean for v in df.wl], :wl_base)

    # each node has a different wind load, so we need to define a model for each node
    # for each model a wind-load.dat needs to be created

    wl_node1_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 50.0, 4, 0.78, 2.6425, 0.49, 0.29, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node1)

    wl_node2_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 47.5, 4, 1.25, 5.285, 0.68, 0.57, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node2)

    wl_node3_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 45.0, 4, 1.25, 8.005, 0.68, 0.57, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node3)

    wl_node4_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 42.5, 4, 1.25, 5.285, 0.68, 0.57, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node4)

    wl_node5_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 40.42, 4, 0.49, 2.6, 0.49, 0.0, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node5)

    wl_node6_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 39.58, 4, 0.77, 2.595933, 0.77, 0.0, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node6)

    wl_node7_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 37.52, 4, 1.1, 5.323556, 1.1, 0.0, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node7)

    wl_node8_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 35.10, 4, 1.13, 3.371424, 1.13, 0.0, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node8)

    wl_node9_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 32.68, 4, 1.16, 6.262081, 1.16, 0.0, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node9)

    wl_node10_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 30.21, 4, 1.37, 7.002028, 1.37, 0.0, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node10)

    wl_node11_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 27.08, 4, 2.3, 11.339223, 2.3, 0.0, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node11)

    wl_node12_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 22.69, 4, 2.92, 17.0680335, 2.92, 0.0, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node12)

    wl_node13_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 17.65, 4, 3.06, 19.102884, 3.06, 0.0, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node13)

    wl_node14_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 12.61, 4, 3.2, 21.1377345, 3.2, 0.0, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node14)

    wl_node15_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 7.56, 4, 3.82, 23.172585, 3.82, 0.0, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node15)

    wl_node16_model = Model(df -> [
        [calculate_wind_force(0, v, 53.0, 2.52, 4, 3.95, 25.2074355, 3.95, 0.0, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base
        ], :wl_node16)
        
    # define the external model that runs OpenSees, this is the model that will be run for each sample
    # it uses the sourcedir and sourcefile to find the OpenSees input files, and the workdir to store the results
    ext = ExternalModel(
    sourcedir, sourcefile, [max_abs_disp,
        disp_node1, disp_node2, disp_node3, disp_node4, disp_node5, disp_node6, disp_node7, disp_node8,
        disp_node9, disp_node10, disp_node11, disp_node12, disp_node13, disp_node14, disp_node15, disp_node16, disp_node17,
        sim_time], opensees; workdir=workdir, formats=numberformats)
        
    # this needs to include the models that are needed and that you want to use
    models = [
    wl_model, base_wl_model,
    wl_node1_model, wl_node2_model, wl_node3_model, wl_node4_model,
    wl_node5_model, wl_node6_model, wl_node7_model, wl_node8_model,
    wl_node9_model, wl_node10_model, wl_node11_model, wl_node12_model,
    wl_node13_model, wl_node14_model, wl_node15_model, wl_node16_model,
    ext]
end

# the capacity is part of the limit state function
capacity = 0.1 #maximum allowed displacements in m

DOWEWANT_MC = false
DOWEWANT_SUS = true

if DOWEWANT_MC
    # this runs the reliability analysis, in this case a Monte Carlo simulation with N samples
    N = 10 # Number of Monte Carlo Samples
    println("Running Monte Carlo simulation with $N samples...")
    pf, std, samples = probability_of_failure(models, df -> capacity .- df.max_abs_disp, [Δt, timeSteps, wl, E, T], MonteCarlo(N))
    println("Probability of failure according to MC: $pf")
end

if DOWEWANT_SUS
    N = 1000 # Number of inisital subset samples
    # Compute probability of failure using Subset Sampling
    subset = UncertaintyQuantification.SubSetSimulation(N, 0.1, 10, Uniform(-0.5, 0.5))
    println("Running Subset simulation with $N initial samples...")
    pf, std, samples = probability_of_failure(models, df -> capacity .- df.max_abs_disp, [Δt, timeSteps, wl, E, T], subset)
    println("Probability of failure according to SuS: $pf")
end

#remove all processes from the Distributed stuff
rmprocs()

######################################### plotting section ################################################################ 
if DOWEWANTplots
    using Plots

    # plotting of PSD function
    plotpsdq = false
    if plotpsdq && length(procs()) == 1
        plot(fdisc, psdvalues, 
            xscale=:log10, 
            #xticks=collect(0:5:50),
            yscale=:log10, 
            label="von Karman Spectrum", 
            xlabel="Frequency (Hz)", 
            ylabel="Power Spectral Density (m²/s²/Hz)", 
            title="von Karman PSD for Longitudinal Wind Component"
        )
    end

    #without normalization
    #which sample should be plotted?
    nmc = 5
    #overall maximum tip displacement that is reached
    max_tip_disp = maximum(samples.max_abs_disp)
    println("Maximum tip displacement: $max_tip_disp m")
    #With following material properties:
    println("Young's modulus: $(samples.E[nmc]) Pa")
    println("Second moment of area: according to section m^4")
    println("Cross sectional area according to section m^2")
    # Maximum occuring wind speed
    max_wind_speed = maximum(abs.(samples.wl_base[nmc]))
    println("Maximum wind speed: $max_wind_speed m/s")
    #Maximum occuring wind load at node 2
    max_wl_node2 = maximum(abs.(samples.wl_node2[nmc]))
    println("Maximum wind load at node 2: $max_wl_node2 N")

    pdisp = plot(t, samples.wl_base[nmc]; label="wind speed in m/s", xlabel="time in s", ylabel="wind speed and displacement", legend=:topright)
    plot!(samples.sim_time[nmc], samples.disp_node1[nmc]; label="Displacement at top node in m", linewidth=2)
    plot!(t, samples.wl_node2[nmc]; label="Wind load at node 2 in kN", linewidth=2)
    plot!(t, samples.wl_node1[nmc]; label="Wind load at node 1 (top) in KN", linewidth=2)

        # Normalized wind speed and displacement for sample nmc
        norm_wind_speed = samples.wl_base[nmc] ./ maximum(abs.(samples.wl_base[nmc]))
        norm_disp = samples.disp_node1[nmc] ./ maximum(abs.(samples.disp_node1[nmc]))

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

    if N > 10
        legendval = false
    else
        legendval = :bottomright
    end


    if N <= 1000
        # plot all wind loads at node 2 for all samples
        pwl_node2 = plot(t, samples.wl_node2[:];
            label="Wind Load at Node 2 (kN)", 
            xlabel="Time (s)", 
            ylabel="Wind Load (kN)", 
            title="Wind Loads at Node 2 for All Samples",
            legend=legendval
        )

        # plot all wind speed signals for all samples
        pwl_base = plot(t, samples.wl_base[:];
            label="Wind Speed (m/s)", 
            xlabel="Time (s)", 
            ylabel="Wind Speed (m/s)", 
            title="Wind Speeds for All Samples",
            legend=legendval
        )
        display(pwl_base)
        display(pwl_node2)
    end

    # Histogram of the maximum absolute displacement
    max_abs_disp_values = samples.max_abs_disp
    phist_maxdisp = histogram(max_abs_disp_values, 
        bins=30, 
        xlabel="Maximum Displacement (m)", 
        ylabel="Frequency", 
        title="Histogram of Maximum Absolute Displacement",
        legend=false
    )

    # Histograms for all realized wind loads at node 2 and node 1
    wl_node2_values = vcat(samples.wl_node2[:]...)
    wl_node1_values = vcat(samples.wl_node1[:]...)

    # Create a single plot with both histograms on it
    phist_wl = histogram(wl_node2_values, 
        bins=30, 
        xlabel="Wind Load (kN)", 
        ylabel="Frequency",
        label ="Node 2",
        title="Histogram of Wind Loads",
        legend=true # Or simply omit this line, as 'true' is often the default
    )

    # Add the second histogram to the same plot
    histogram!(phist_wl, wl_node1_values, 
        bins=30, 
        label="Node 1"
    )

    # aggregate the histograms for disp and wind loads at node 2 and node 1
    # aggregate the histograms for disp and wind loads at node 2 and node 1
    phist = plot(phist_maxdisp, phist_wl; layout=(2, 1), size=(800, 600), title="Histograms")

    if DOWEWANTseeplots
        display(phist)
        display(pdisp)
    end

    if DOWEWANTsaveplots
        png(pdisp, joinpath(workdir, time_stamp*"displacement_windload_timeseries.png"))
        png(phist, joinpath(workdir, time_stamp*"histograms.png"))
    end
end

if DOWEWANTsaveresults
    #### Save results
    # Save the hyperparameters
    hyperparameters = Dict(
        "Nsamples" => N,
        "pf" => pf,
        "jobid" => jobid,
        "MC" => DOWEWANT_MC,
        "SuS" => DOWEWANT_SUS,
        "capacity" => capacity,
    )
    open(joinpath(workdir, time_stamp*"_model_2_metadata.json"), "w") do f
        JSON.print(f, hyperparameters, 4)
    end
    # Define the filename for your JLD2 file
    filename = time_stamp*"_model_2_samples.jld2"
    # Save the DataFrame to the JLD2 file
    @save filename samples
end

####################################################################################################################################
toc_total = time() - tic_total
println("Total elapsed time: $(round(toc_total, digits=3))s")
@printf "[%s]: " Dates.format(now(), "dd/mm/yy, HH:MM:SS")
println("Script finished.")