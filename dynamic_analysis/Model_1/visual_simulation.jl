using DelimitedFiles
using Plots

# Y-coordinates of the nodes
# Base at 0 m, Node 2 at 40 m, Node 3 at 50 m
node_y = [50.0, 40.0, 0.0]   
n_nodes = 3

# sample index to visualize
nmc = 1

# Get the time series and displacements from the DataFrame
times = samples.sim_time[nmc]

displacements = [
    i == 1 ? samples[!, :disp][nmc] : samples[!, Symbol("disp_node$(i)")][nmc]
    for i in 1:n_nodes]

anim = @animate for t_idx in 1:length(times)
    node_x = [displacements[i][t_idx] for i in 1:n_nodes]
    plot(node_x, node_y, marker=:circle, xlabel="X (m)", ylabel="Y (m)", title="Tower Deformation",
         xlims=(-0.01,0.01), ylims=(0, 55), legend=false, size=(600,800))
end

gif(anim, "C:/IP/dynamic_analysis/Model_1/tower_simulation.gif", fps=20)
