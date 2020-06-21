function build_and_solve_opt_model(data, Σ, Σ½, I_m, node_release_set, cluster_sets, cluster_set, set, model_type)
    # noise dimention
    n = size(Σ,1)
    # initialize model
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0))
    # variable declaration
    @variable(model, p[1:data[:Nb]])
    @variable(model, θ[1:data[:Nb]])
    model_type == "chance_constrained" ? @variable(model, 𝐏[1:data[:Nb],1:n]) : NaN
    model_type == "chance_constrained" ? @variable(model, 𝚯[1:data[:Nb],1:n]) : NaN
    model_type == "chance_constrained" ? @variable(model, obj_std) : NaN
    model_type == "chance_constrained" ? @variable(model, sol_std) : NaN
    # deterministic constraints
    @constraint(model, node_bal_det,            data[:B]*θ .== p - data[:d])
    @constraint(model, ref_node_angle_det,      θ[data[:refnode]] == 0)
    @constraint(model, ref_node_supply_det,     p[data[:refnode]] == 0)
    if model_type == "deterministic"
        @constraint(model, p_max_det,                       p .<= data[:p̅])
        @constraint(model, p_min_det,                       p .>= data[:p̲])
        @constraint(model, flow_min_det[l=1:data[:Nl]],     -data[:f̅][l] <= data[:β][l]*(θ[s(l)]-θ[r(l)]))
        @constraint(model, flow_max_det[l=1:data[:Nl]],     data[:β][l]*(θ[s(l)]-θ[r(l)]) <= data[:f̅][l])
    end
    if model_type == "chance_constrained"
        # chance constraints
        if set[:ref_type] == "analytic"
            @constraint(model, p_max_cc[i=1:data[:Nb]],
                [data[:p̅][i]-p[i];sqrt(2/(9*set[:η]))*Σ½*𝐏[i,:]] in SecondOrderCone())
            @constraint(model, p_min_cc[i=1:data[:Nb]],
                [p[i]-data[:p̲][i];sqrt(2/(9*set[:η]))*Σ½*𝐏[i,:]] in SecondOrderCone())
            @constraint(model, flow_min_cc[l=1:data[:Nl]],
                [data[:f̅][l] - data[:β][l]*(θ[s(l)]-θ[r(l)]);sqrt(2/(9*set[:η]))*data[:β][l]*Σ½*(𝚯[s(l),:]-𝚯[r(l),:])] in SecondOrderCone())
            @constraint(model, flow_max_cc[l=1:data[:Nl]],
                [data[:f̅][l] - data[:β][l]*(θ[r(l)]-θ[s(l)]);sqrt(2/(9*set[:η]))*data[:β][l]*Σ½*(𝚯[r(l),:]-𝚯[s(l),:])] in SecondOrderCone())
        end
        if set[:ref_type] == "sample"
            data[:Nb] ∉ [57, 118] && set[:query_type] == "identity" ? ξ̄ = get_sample_vertices(n,set,I_m,cluster_sets) : NaN                        # use (Margellos et al., 2014)  for small and medium systems
            data[:Nb] ∈ [57, 118] && set[:query_type] == "identity" ? ξ̄ = rand(Laplace(0,set[:α]*set[:Δ]/set[:ε]),1000,n) : NaN                    # use (Campi and Garatti, 2008) for large systems with 1000 samples
            set[:query_type] == "linear" ? ξ̄ = get_sample_vertices(n,set,I_m,cluster_sets) : NaN
            @constraint(model, p_con_sample[i=1:data[:Nb], v=1:size(ξ̄)[1]], data[:p̲][i] <= p[i] + 𝐏[i,:]'*ξ̄[v,:] <= data[:p̅][i])
            @constraint(model, f_con_sample[l=1:data[:Nl], v=1:size(ξ̄)[1]], -data[:f̅][l] <= data[:β][l]*(θ[s(l)] + (𝚯[s(l),:]'*ξ̄[v,:])' - θ[r(l)] - (𝚯[r(l),:]'*ξ̄[v,:])') <= data[:f̅][l])
        end
        # equality stochastic constraints
        @constraint(model, as_con,               data[:B]*𝚯 .== 𝐏)
        @constraint(model, ref_node_angle_cc,    𝚯[data[:refnode],:] .== 0)
        @constraint(model, ref_node_supply_cc,   𝐏[data[:refnode],:] .== 0)
        # query constraints
        if set[:query_type] == "identity"
            @constraint(model, q_con,    I_m*diag(𝐏[node_release_set,:]) .== diag(I_m))
            @constraint(model, q_con_,   I_m*(𝐏[node_release_set,:] .- diagm(diag(𝐏[node_release_set,:]))) .== 0)
        end
        if set[:query_type] == "linear"
            # in-cluster recourse constraint
            @constraint(model, q_con[p=1:n], sum(𝐏[i,p] for i in cluster_set if i ∈ cluster_sets[p]) == 1)
            # not-in-cluster recourse constraint
            for k in 1:n, i in cluster_set
                i ∉ cluster_sets[k] ? @constraint(model, 𝐏[i,k] .== 0) : NaN
            end
        end
        # objective standard deviation constraint
        @constraint(model, [obj_std;Σ½*𝐏'*data[:c1]] in SecondOrderCone())
        # solution standard deviation constraint
        @constraint(model, [sol_std;(𝐏+𝚯)*Σ½*ones(n)] in SecondOrderCone())
    end

    # objective function
    if model_type == "deterministic"
        set[:cost_function] == "linear" ? @objective(model, Min, data[:c1]'p) : @objective(model, Min, data[:c1]'p + p'diagm(data[:c2])p)
    else
        if set[:obj_var_control] == false && set[:sol_var_control] == false
            set[:cost_function] == "linear" ? @objective(model, Min, data[:c1]'p) : @objective(model, Min, data[:c1]'p + p'diagm(data[:c2])p + tr(𝐏'*diagm(data[:c2])*𝐏*Σ))
        end
        if set[:obj_var_control] == true
            @objective(model, Min, (1-set[:φ]) * data[:c1]'p + set[:φ]*obj_std)
        end
        if set[:sol_var_control] == true
            @objective(model, Min, (1-set[:φ]) * data[:c1]'p + set[:φ]*sol_std)
        end
    end
    # solve model
    optimize!(model)
    status = termination_status(model)
    @info("$(model_type) model terminates with status $(status) in $(round(MOI.get(model, MOI.SolveTime()),digits=5)) seconds")
    # prepare the results
    model_type == "deterministic" && set[:cost_function] == "linear"            ?   cost = data[:c1]'JuMP.value.(p) : NaN
    model_type == "deterministic" && set[:cost_function] == "quadratic"         ?   cost = data[:c1]'JuMP.value.(p) + JuMP.value.(p)'diagm(data[:c2])JuMP.value.(p) : NaN
    model_type == "chance_constrained" && set[:cost_function] == "linear"       ?   cost = data[:c1]'JuMP.value.(p) : NaN
    model_type == "chance_constrained" && set[:cost_function] == "quadratic"    ?   cost = data[:c1]'JuMP.value.(p) + JuMP.value.(p)'diagm(data[:c2])JuMP.value.(p) + tr(JuMP.value.(𝐏)'*diagm(data[:c2])*JuMP.value.(𝐏)*Σ) : NaN
    solution = Dict("p" => JuMP.value.(p),
                    "θ" => JuMP.value.(θ),
                    "𝐏" => model_type == "chance_constrained" ? JuMP.value.(𝐏) : NaN,
                    "𝚯" => model_type == "chance_constrained" ? JuMP.value.(𝚯) : NaN,
                    "cost" => cost,
                    "cost_var" => model_type == "chance_constrained" ? JuMP.value.(obj_std)^2 : NaN,
                    "sol_var" => model_type == "chance_constrained" ? JuMP.value.(sol_std)^2 : NaN
                )
    return model, status, solution
end
