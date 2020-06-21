function experiment_setting(inputs)
    # def share of nodes to be released
    γ = inputs["share_node"]
    # def query type
    query_type = inputs["query"]
    # def privacy parameters
    α = inputs["adjacency"]; ε = inputs["epsilon"]; Δ = 1;
    # def method for chance constraints reformulation
    ref_type = inputs["reformulation"]
    # def constraint violation tolerance
    η = inputs["violation"];
    # def cost function type
    cost_function = inputs["cost_fun"]
    # def number of clusters for linear query analysis
    cluster_num = inputs["cluster_num"]
    # def variance settings
    obj_var_control = inputs["obj_var_control"]
    sol_var_control = inputs["sol_var_control"]
    var_penalty = inputs["var_penalty"]
    set = Dict(:γ => γ, :query_type => query_type, :α => α, :ε => ε,
                    :Δ => Δ, :ref_type => ref_type, :η => η,
                    :cost_function => cost_function, :β => 0.01,
                    :cluster_num => cluster_num,
                    :obj_var_control => obj_var_control,
                    :sol_var_control => sol_var_control,
                    :φ => var_penalty)
    return set
end
