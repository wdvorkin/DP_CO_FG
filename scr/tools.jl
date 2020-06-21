function get_sample_vertices(n,set,I_m,cluster_sets)
    if set[:query_type] == "identity"
        S = Int(round(1/(set[:η]) * ℯ/(ℯ-1) * (log(1/set[:β]) + 2*n - 1)))
        sample_model = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0))
        ξ = zeros(n,S); ξ̅ = zeros(n); ξ̲ = zeros(n);
        for s in 1:S, i in 1:n
            I_m[i,i] == 1 ? ξ[i,s] = rand(Laplace(0,set[:α] * set[:Δ] / set[:ε]),1)[1] : NaN
        end
        for i in 1:n
            ξ̅[i] = maximum(ξ[i,:])
            ξ̲[i] = minimum(ξ[i,:])
        end
        @variable(sample_model, ξ[1:n])
        @constraint(sample_model, con[i=1:n], ξ̲[i] <= ξ[i] <= ξ̅[i])
    end
    if set[:query_type] == "linear"
        S = Int(round(1/(set[:η]) * ℯ/(ℯ-1) * (log(1/set[:β]) + 2*length(cluster_sets) - 1)))
        sample_model = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0))
        ξ̅ = zeros(length(cluster_sets)); ξ̲ = zeros(length(cluster_sets));
        ξ = rand(Laplace(0,set[:α] * set[:Δ] / set[:ε]),length(cluster_sets),S)
        for i in 1:length(cluster_sets)
            ξ̅[i] = maximum(ξ[i,:])
            ξ̲[i] = minimum(ξ[i,:])
        end
        @variable(sample_model, ξ[1:length(cluster_sets)])
        @constraint(sample_model, con[i=1:length(cluster_sets)], ξ̲[i] <= ξ[i] <= ξ̅[i])
    end
    # create a polyhedron from the sample set
    poly = polyhedron(sample_model, CDDLib.Library(:exact))
    # obtain V-representation of the polyhedron
    vr = vrep(poly)
    vr = MixedMatVRep(vr)
    p = vr.V
    return float(p)
end
function get_release_set(data,set)
    sett=[1:data[:Nb];][rand(1:data[:Nb],max(Int(floor(data[:Nb]*set[:γ]*1.1)),1))]
    set_unique = unique!(sett)
    set_unique_no_ref = setdiff!(set_unique,data[:refnode])
    isempty(set_unique_no_ref) == true ? push!(set_unique_no_ref,1) : NaN
    set_γ = set_unique[1:min(max(Int(floor(data[:Nb]*set[:γ])),1),length(set_unique_no_ref))]
    return set_γ
end
function identity_matrix(node_release_set)
    return diagm(ones(length(node_release_set)))
end
function noise_data(set,I_m,cluster_sets)
    σ = sqrt(2*(set[:α]*set[:Δ]/set[:ε])^2)
    set[:query_type] == "identity" ? Σ  = zeros(size(I_m)) .+ I_m * σ^2 : NaN
    set[:query_type] == "identity" ? Σ½ = zeros(size(I_m)) .+ I_m * σ : NaN
    set[:query_type] == "linear"   ? Σ  = diagm(ones(length(cluster_sets)) .* σ^2) : NaN
    set[:query_type] == "linear"   ? Σ½ = diagm(ones(length(cluster_sets)) .* σ) : NaN
    return Σ, Σ½
end
function cluster_node_allocation(set,node_release_set)
    if set[:query_type] == "linear"
        nodes_in_cluster = Int(floor(length(node_release_set)/set[:cluster_num]))
        cluster_sets=collect(Iterators.partition(node_release_set, nodes_in_cluster))
        length(cluster_sets) > set[:cluster_num] ? cluster_sets = cluster_sets[1:end-1] : NaN

        cluster_set = []
        for k in 1:length(cluster_sets), i in 1:length(cluster_sets[k])
            push!(cluster_set,cluster_sets[k][i])
        end

        return cluster_sets, cluster_set
    else
        return NaN, NaN
    end
end
s(l) = Int(data[:s_end][l])
r(l) = Int(data[:r_end][l])
