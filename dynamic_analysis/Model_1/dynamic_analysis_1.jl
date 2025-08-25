using Distributed, Printf, Dates, Plots
time_stamp = Dates.format(now(),  "yyyymmdd_HHMMSS")
@printf "[%s]: " time_stamp
println("Reliability Analysis of a Dynamic Tower-like Structure under Wind loads")
tic_total = time()

addprocs(6)

DOWEWANTplots = true

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
    E = RandomVariable(Uniform(2.0e8, 2.2e8), :E)   # Young's modulus
    t = collect(0:Δt.value:40)                      # time vector
    timeSteps = Parameter(length(t), :timeSteps)    # number of time steps

    # von Karman spectrum
    fdisc = collect(0.01:0.01:50)
    # roughly derived parameters from https://www.sciencedirect.com/science/article/pii/S0167610520302725
    σ_u = 1 # standard deviation of wind speed in m/s
    L_u = 50.8 # length scale in m
    v_mean = 10.0 # mean wind speed in m/s
    karmanfunc = (f) -> σ_u^2 * (4*L_u / v_mean) / ((1 + (2*π*f * L_u / v_mean)^2)^(5/6))
    psdvalues = karmanfunc.(fdisc)

    #coinstruct the empirical PSD representing the von Karman spectrum
    epWindPSD = EmpiricalPSD(fdisc, psdvalues)

    wl = SpectralRepresentation(epWindPSD, t, :wl)
    wl_model = StochasticProcessModel(wl)

    # directory where the OpenSees input files are stored

    if hostname == "IRZ-PC1337"
        sourcedir = joinpath(pwd(), "C:/OMNISSIAH/Julia/2025-08-pauline-git/dynamic_analysis/Model_1")
    else
        sourcedir = joinpath(pwd(), "C:/IP/dynamic_analysis/Model_1")
    end
    sourcefile = ["FEM_1.tcl", "wind-speed.dat", "wind-load-abs.dat", "wind-load-node2.dat", "wind-load-node1.dat"]

    numberformats = Dict(:dt => ".8e", :wl => ".8e", :E => ".8e", :Iz => ".8e")

    # case where the results are stored
    workdir = joinpath(pwd(), "workdir-1")

    # define the models that extract the results from the OpenSees output files, we are interested in the maximum absolute displacement, the displacement at the top node and the simulation time
    max_abs_disp = Extractor(base -> begin
        file = joinpath(base, "displacement.out")
        data = DelimitedFiles.readdlm(file, ' ')

        return maximum(abs.(data[:, 2]))
    end, :max_abs_disp)

    # this is the displacement time history, only needed for plotting
    disp = Extractor(base -> begin
        file = joinpath(base, "displacement.out")
        data = DelimitedFiles.readdlm(file, ' ')

        return data[:, 2]
    end, :disp)

    # this is the simulation time, also only needed for plotting
    sim_time = Extractor(base -> begin
        file = joinpath(base, "displacement.out")
        data = DelimitedFiles.readdlm(file, ' ')

        return data[:, 1]
    end, :sim_time)


    opensees = Solver(
        "OpenSees", # path to OpenSees binary
        "FEM_1.tcl";
        args="", # (optional) extra arguments passed to the solver
    )


    # alters the model of the wind speed to habe a baseline wind speed of v_mean
    base_wl_model = Model(df -> [v.+v_mean for v in df.wl], :wl_base)

    # each node has a different wind load, so we need to define a model for each node
    # for each model a wind-load.dat needs to be created

    wl_node2_model = Model(df -> [
        [calculate_wind_force(
            0, v, 53.0, 40.0, 4, 24.78, 144.0, 24.78, 0.0, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base], 
        :wl_node2)  

    wl_node1_model = Model(df -> [
        [calculate_wind_force(                      
            0, v, 53.0, 50.0, 4, 5.02, 23.82, 3.02, 2.0, 0.0)
        for v in wl_vec]
        for wl_vec in df.wl_base],
        :wl_node1)

        
    # define the external model that runs OpenSees, this is the model that will be run for each sample
    # it uses the sourcedir and sourcefile to find the OpenSees input files, and the workdir to store the results
    ext = ExternalModel(
        sourcedir, sourcefile, [max_abs_disp, disp, sim_time], opensees; workdir=workdir, formats=numberformats
    )

    # this needs to include the models that are needed and that you want to use
    models = [wl_model, base_wl_model, wl_node2_model, wl_node1_model, ext]
end

# this runs the reliability analysis, in this case a Monte Carlo simulation with N_MC samples
N_MC = 100 # Number of Monte Carlo Samples
println("Running Monte Carlo simulation with $N_MC samples...")
# this part: df -> 200 .- df.max_abs_disp
# actually defines the performance function also known as limit state function which is evaluated for each of the samples, if you use the same record this of course does not make any sense
pf, mc_std, samples = probability_of_failure(models, df -> 1 .- df.max_abs_disp, [Δt, timeSteps, wl, E, Iz], MonteCarlo(N_MC))
println
println("Probability of failure: $pf")

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
    plot!(samples.sim_time[nmc], samples.disp[nmc]; label="Displacement at top node in m", linewidth=2)
    plot!(samples.sim_time[nmc], samples.wl_node2[nmc]; label="Wind load at node 2 in kN", linewidth=2)
    plot!(samples.sim_time[nmc], samples.wl_node1[nmc]; label="Wind load at node 1 (top) in KN", linewidth=2)
    # This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

    if N_MC > 10
        legendval = false
    else
        legendval = :bottomright
    end

    if N_MC <= 1000
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
    phist = plot(phist_maxdisp, phist_wl; layout=(2, 1), size=(800, 600), title="Histograms")

    display(phist)

    display(pdisp)

end

rmprocs()

####################################################################################################################################
toc_total = time() - tic_total
println("Total elapsed time: $(round(toc_total, digits=3))s")
@printf "[%s]: " Dates.format(now(), "dd/mm/yy, HH:MM:SS")
println("Script finished.")