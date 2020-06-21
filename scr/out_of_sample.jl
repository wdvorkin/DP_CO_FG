function out_of_sample(data, model_de, node_release_set, solution_de, solution_cc, set, cluster_sets)
    S = 1000
    infeas_de = zeros(S)
    infeas_cc = zeros(S)
    # get noise samples
    if set[:query_type] == "identity"
        ξ = zeros(data[:Nb],S)
        for s in 1:S, i in 1:data[:Nb]
            i ∈ node_release_set ? ξ[i,s] = rand(Laplace(0,set[:α]*set[:Δ]/set[:ε]),1)[1] : NaN
        end
    elseif set[:query_type] == "linear"
        ξ = rand(Laplace(0,set[:α]*set[:Δ]/set[:ε]),length(cluster_sets),S)
    end
    # chance-constrained solution analysis
    for k in 1:S
        model_de_ofs = copy(model_de)
        set_optimizer(model_de_ofs, Mosek.Optimizer)
        set_optimizer_attributes(model_de_ofs, "LOG" => 0)
        p = model_de_ofs[:p]
        set[:query_type] == "identity" ? @constraint(model_de_ofs, var_fix[i=node_release_set], p[i] == solution_cc["p"][i] + solution_cc["𝐏"][i,:]'*ξ[node_release_set,k]) : NaN
        set[:query_type] == "linear"   ? @constraint(model_de_ofs, var_fix[c=1:length(cluster_sets)], sum(p[i] for i in cluster_sets[c]) == sum(solution_cc["p"][i] + solution_cc["𝐏"][i,c]*ξ[c,k] for i in cluster_sets[c])) : NaN
        optimize!(model_de_ofs)
        "$(termination_status(model_de_ofs))" != "OPTIMAL" ? infeas_cc[k] = 1 : infeas_cc[k] = 0
    end
    # deterministic solution analysis
    for k in 1:S
        model_de_ofs = copy(model_de)
        set_optimizer(model_de_ofs, Mosek.Optimizer)
        set_optimizer_attributes(model_de_ofs, "LOG" => 0)
        p = model_de_ofs[:p]
        set[:query_type] == "identity" ? @constraint(model_de_ofs, var_fix[i=node_release_set], p[i] == solution_de["p"][i] + ξ[i,k]) : NaN
        set[:query_type] == "linear"   ? @constraint(model_de_ofs, var_fix[c=1:length(cluster_sets)], sum(p[i] for i in cluster_sets[c]) == sum(solution_de["p"][i] for i in cluster_sets[c]) + ξ[c,k]) : NaN
        optimize!(model_de_ofs)
        "$(termination_status(model_de_ofs))" != "OPTIMAL" ? infeas_de[k] = 1 : infeas_de[k] = 0
    end
    return sum(infeas_de)/S, sum(infeas_cc)/S
end
