using Plots

node_y = [50.00, 47.50, 45.00, 42.50, 40.42, 39.58, 37.52, 35.10, 32.68, 30.21, 27.08, 22.69, 17.65, 12.61, 7.56, 2.52]
n_nodes = 16

# Beispiel: nmc ist der Index des gewünschten Samples
nmc = 1

# Hole die Zeitreihe und die Verschiebungen aus dem DataFrame
times = samples.sim_time[nmc]
displacements = [samples[!, Symbol("disp_node$(i)")][nmc] for i in 1:n_nodes]

anim = @animate for t_idx in 1:length(times)
    node_x = [displacements[i][t_idx] for i in 1:n_nodes]
    plot(node_x, node_y, marker=:circle, xlabel="X (m)", ylabel="Y (m)", title="Tower Deformation",
         xlims=(-0.1,0.1), ylims=(0, 55), legend=false, size=(600,800))
end

gif(anim, "tower_simulation.gif", fps=20)