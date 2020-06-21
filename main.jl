#!/usr/bin/env julia
using JuMP, Mosek, MosekTools, MathOptInterface
using PowerModels, DataStructures, DataFrames, CSV, ArgParse
using CDDLib, Polyhedra, LinearAlgebra
using Distributions

# load scripts
include("scr/build_opt.jl")
include("scr/exp_settings.jl")
include("scr/load_data.jl")
include("scr/out_of_sample.jl")
include("scr/tools.jl")
include("scr/parsing.jl")

# extract inputs
inputs = parse_commandline()
# extract network data
caseID = "testbeds/" * inputs["case"] *".m"; data = load_data(caseID) # test function
# extract experiment setting
set = experiment_setting(inputs)
# def output frames
results_exp = DataFrame(it = Int[], opt_loss_a = Any[], opt_loss_s = Any[], infeas_de = Any[], infeas_cc_a = Any[], infeas_cc_s = Any[])
# run exp_num number of experiments
for iter in 1:inputs["exp_num"]
    # load random network data
    data = load_data(caseID)
    # get node release set
    node_release_set=get_release_set(data,set)
    # get identity matrix for the indentity query
    I_m = identity_matrix(node_release_set)
    # get node clusters for the linear query
    cluster_sets, cluster_set = cluster_node_allocation(set,node_release_set)
    # get noise parameters
    Σ, Σ½ = noise_data(set,I_m,cluster_sets)
    # solve deterministic model
    model_de, status_de, solution_de = build_and_solve_opt_model(data, Σ, Σ½, I_m, node_release_set, cluster_sets, cluster_set, set, "deterministic")
    # # solve chance-constrained model with analytic reformulation
    set[:ref_type] = "analytic"
    model_cc_a, status_cc_a, solution_cc_a = build_and_solve_opt_model(data, Σ, Σ½, I_m, node_release_set, cluster_sets, cluster_set, set, "chance_constrained")
    # run out-of-sample analysis
    inf_de, inf_cc_a = out_of_sample(data, model_de, node_release_set, solution_de, solution_cc_a, set, cluster_sets)
    # # solve chance-constrained model with sample reformulation
    set[:ref_type] = "sample"
    model_cc_s, status_cc_s, solution_cc_s = build_and_solve_opt_model(data, Σ, Σ½, I_m, node_release_set, cluster_sets, cluster_set, set, "chance_constrained")
    # run out-of-sample analysis
    inf_de, inf_cc_s = out_of_sample(data, model_de, node_release_set, solution_de, solution_cc_s, set, cluster_sets)
    # save the results
    push!(results_exp,[iter,(abs(solution_de["cost"]-solution_cc_a["cost"])/solution_de["cost"]*100),(abs(solution_de["cost"]-solution_cc_s["cost"])/solution_de["cost"]*100),inf_de,inf_cc_a,inf_cc_s])
    @info("done experiment $(iter)")
end
outdir = "output/" * inputs["case"]; mkpath(outdir);
CSV.write("$(outdir)/results_avr_query_$(set[:query_type])_cluster_num_$(inputs["cluster_num"])_cost_$(set[:cost_function])_alpha_$(set[:α])_epsilon_$(set[:ε])_violation_$(set[:η]).csv",results_exp)
