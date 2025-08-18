using UncertaintyQuantification
using DelimitedFiles
using Plots

include("C:/IP/dynamic_analysis/calculate_wind_force.jl")

# parameter definition, they can be changed but are the same for each sample
Δt = Parameter(0.02, :dt)                       # time step size
E = RandomVariable(Uniform(2.0e11, 2.2e11), :E)   # Young's modulus
Iz= RandomVariable(Uniform(1e-6, 1e-3), :Iz)    # second moment of area
t = collect(0:Δt.value:40)                      # time vector
timeSteps = Parameter(length(t), :timeSteps)    # number of time steps


ω = collect(0:0.6:150)                          # frequency vector

cp = CloughPenzien(ω, 0.1, 0.8π, 0.6, 8π, 0.6) 

# this is the stochastic signal generation model, which is disabled in this script but we want to use it later
wl = SpectralRepresentation(cp, t, :wl)
wl_model = StochasticProcessModel(wl)

#C:\Users\pbett\OneDrive\Dokumente\IP\dymanic_analysis\tower_uq\wind_loads\wind-load-node16.dat

# directory where the OpenSees input files are stored
sourcedir = joinpath(pwd(), "C:/IP/dynamic_analysis/Model_2")
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
    "wind-load-node16.dat"
]

numberformats = Dict(:dt => ".8e", :wl => ".8e", :E => ".8e", :Iz => ".8e")

# case where the results are stored
workdir = joinpath(pwd(), "workdir-2")

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
    "FEM_2.tcl";
    args="", # (optional) extra arguments passed to the solver
)


# alters the model of the wind load to contain only absolut values
abs_wl_model = Model(df -> [abs.(v) for v in df.wl], :wl_abs)

# each node has a different wind load, so we need to define a model for each node
# for each model a wind-load.dat needs to be created

# factor to scale wl to a realistic size
s = 1.0

# Set a threshold for "extremely large" values
max_reasonable = 1e3*s  # adjust as needed

# structure of the model for each node  
 #create new wind load model for node 1
 #clean up the wind load, if it is not finite set to 0.0, if it is larger than the threshold set it to the threshold
 #for each value in the wind load vector
 #apply the wind force calculation function to each value in the wind load vector
 #for each wind load vector in the absolute wind load model
 #save in datafram as column wl_node1

# Model for node 1
wl_node1_model = Model(df -> [
    [isfinite(wl_node1) ? (wl_node1 > max_reasonable ? max_reasonable : wl_node1) : 0.0
     for v in wl_vec
     for wl_node1 = [calculate_wind_force(
         0, v*s, 53.0, 50.0, 4, 0.78, 2.6425, 0.49, 0.29, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node1)

# Model for node 2
wl_node2_model = Model(df -> [
    [isfinite(wl_node2) ? (wl_node2 > max_reasonable ? max_reasonable : wl_node2) : 0.0
    for v in wl_vec
    for wl_node2 = [calculate_wind_force(
    0, v*s, 53.0, 47.5, 4, 1.25, 5.285, 0.68, 0.57, 0.0).F_wp ]]
    for wl_vec in df.wl_abs], 
    :wl_node2)

# Model für node 3
wl_node3_model = Model(df -> [
    [isfinite(wl_node3) ? (wl_node3 > max_reasonable ? max_reasonable : wl_node3) : 0.0
     for v in wl_vec
     for wl_node3 = [calculate_wind_force(
         0, v*s, 53.0, 45.0, 4, 1.25, 8.005, 0.68, 0.57, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node3)

# Model für node 4
wl_node4_model = Model(df -> [
    [isfinite(wl_node4) ? (wl_node4 > max_reasonable ? max_reasonable : wl_node4) : 0.0
     for v in wl_vec
     for wl_node4 = [calculate_wind_force(
         0, v*s, 53.0, 42.5, 4, 1.25, 5.285, 0.68, 0.57, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node4)

# Model für node 5
wl_node5_model = Model(df -> [
    [isfinite(wl_node5) ? (wl_node5 > max_reasonable ? max_reasonable : wl_node5) : 0.0
     for v in wl_vec
     for wl_node5 = [calculate_wind_force(
         0, v*s, 53.0, 40.42, 4, 0.49, 2.6, 0.49, 0.0, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node5)

# Model für node 6
wl_node6_model = Model(df -> [
    [isfinite(wl_node6) ? (wl_node6 > max_reasonable ? max_reasonable : wl_node6) : 0.0
     for v in wl_vec
     for wl_node6 = [calculate_wind_force(
         0, v*s, 53.0, 39.58, 4, 0.77, 2.595933, 0.77, 0.0, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node6)

# Model für node 7
wl_node7_model = Model(df -> [
    [isfinite(wl_node7) ? (wl_node7 > max_reasonable ? max_reasonable : wl_node7) : 0.0
     for v in wl_vec
     for wl_node7 = [calculate_wind_force(
         0, v*s, 53.0, 37.52, 4, 1.1, 5.323556, 1.1, 0.0, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node7)

# Model für node 8
wl_node8_model = Model(df -> [
    [isfinite(wl_node8) ? (wl_node8 > max_reasonable ? max_reasonable : wl_node8) : 0.0
     for v in wl_vec
     for wl_node8 = [calculate_wind_force(
         0, v*s, 53.0, 35.10, 4, 1.13, 3.371424, 1.13, 0.0, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node8)

# Model für node 9
wl_node9_model = Model(df -> [
    [isfinite(wl_node9) ? (wl_node9 > max_reasonable ? max_reasonable : wl_node9) : 0.0
     for v in wl_vec
     for wl_node9 = [calculate_wind_force(
         0, v*s, 53.0, 32.68, 4, 1.16, 6.262081, 1.16, 0.0, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node9)

# Model für node 10
wl_node10_model = Model(df -> [
    [isfinite(wl_node10) ? (wl_node10 > max_reasonable ? max_reasonable : wl_node10) : 0.0
     for v in wl_vec
     for wl_node10 = [calculate_wind_force(
         0, v*s, 53.0, 30.21, 4, 1.37, 7.002028, 1.37, 0.0, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node10)

# Model für node 11
wl_node11_model = Model(df -> [
    [isfinite(wl_node11) ? (wl_node11 > max_reasonable ? max_reasonable : wl_node11) : 0.0
     for v in wl_vec
     for wl_node11 = [calculate_wind_force(
         0, v*s, 53.0, 27.08, 4, 2.3, 11.339223, 2.3, 0.0, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node11)

# Model für node 12
wl_node12_model = Model(df -> [
    [isfinite(wl_node12) ? (wl_node12 > max_reasonable ? max_reasonable : wl_node12) : 0.0
     for v in wl_vec
     for wl_node12 = [calculate_wind_force(
         0, v*s, 53.0, 22.69, 4, 2.92, 17.0680335, 2.92, 0.0, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node12)

# Model für node 13
wl_node13_model = Model(df -> [
    [isfinite(wl_node13) ? (wl_node13 > max_reasonable ? max_reasonable : wl_node13) : 0.0
     for v in wl_vec
     for wl_node13 = [calculate_wind_force(
         0, v*s, 53.0, 17.65, 4, 3.06, 19.102884, 3.06, 0.0, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node13)

# Model für node 14
wl_node14_model = Model(df -> [
    [isfinite(wl_node14) ? (wl_node14 > max_reasonable ? max_reasonable : wl_node14) : 0.0
     for v in wl_vec
     for wl_node14 = [calculate_wind_force(
         0, v*s, 53.0, 12.61, 4, 3.2, 21.1377345, 3.2, 0.0, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node14)

# Model für node 15
wl_node15_model = Model(df -> [
    [isfinite(wl_node15) ? (wl_node15 > max_reasonable ? max_reasonable : wl_node15) : 0.0
     for v in wl_vec
     for wl_node15 = [calculate_wind_force(
         0, v*s, 53.0, 7.56, 4, 3.82, 23.172585, 3.82, 0.0, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node15)

# Model für node 16
wl_node16_model = Model(df -> [
    [isfinite(wl_node16) ? (wl_node16 > max_reasonable ? max_reasonable : wl_node16) : 0.0
     for v in wl_vec
     for wl_node16 = [calculate_wind_force(
         0, v*s, 53.0, 2.52, 4, 3.95, 25.2074355, 3.95, 0.0, 0.0).F_wp]]
    for wl_vec in df.wl_abs
], :wl_node16)


# define the external model that runs OpenSees, this is the model that will be run for each sample
# it uses the sourcedir and sourcefile to find the OpenSees input files, and the workdir to store the results
ext = ExternalModel(
    sourcedir, sourcefile, [max_abs_disp, disp, sim_time], opensees; workdir=workdir, formats=numberformats
)

# this needs to include the models that are needed and that you want to use
models = [
    wl_model,
    abs_wl_model,
    wl_node1_model,
    wl_node2_model,
    wl_node3_model,
    wl_node4_model,
    wl_node5_model,
    wl_node6_model,
    wl_node7_model,
    wl_node8_model,
    wl_node9_model,
    wl_node10_model,
    wl_node11_model,
    wl_node12_model,
    wl_node13_model,
    wl_node14_model,
    wl_node15_model,
    wl_node16_model,
    ext
]

# this runs the reliability analysis, in this case a Monte Carlo simulation with 1 sample
# this part: df -> 200 .- df.max_abs_disp
# actually defines the performance function also known as limit state function which is evaluated for each of the samples, if you use the same record this of course does not make any sense
pf, mc_std, samples = probability_of_failure(models, df -> 1.4 .- df.max_abs_disp, [Δt, timeSteps, wl, E, Iz], MonteCarlo(10))
println("Probability of failure: $pf")

#with normalization
plot(t, samples.wl_abs[1]./(maximum(abs.(samples.wl_abs[1]))); label="wind speed", xlabel="time in s", ylabel="windspeed and displacement")
#plot!(samples.sim_time[1], samples.disp[1]./(maximum(abs.(samples.disp[1]))); label="Displacement at top node", linewidth=2)

#without normalization
#plot(t, samples.wl_sum[1]; label="wind load", xlabel="time in s", ylabel="windload and displacement")
plot!(samples.sim_time[1], samples.disp[1]; label="Displacement at top node", linewidth=2, alpha=0.8)
# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
