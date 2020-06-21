function load_data(caseID)
    PowerModels.silence()
    net_data = PowerModels.parse_file(caseID)
    # network dimension
    Nb = length(net_data["bus"])
    Nl = length(net_data["branch"])
    # initialize network parameters
    p̅ = zeros(Nb); p̲ = zeros(Nb); d = zeros(Nb); c₁ = zeros(Nb); c₂ = zeros(Nb);
    s_end = zeros(Nl); r_end = zeros(Nl); β = zeros(Nl); f̅ = zeros(Nl); B = zeros(Nb,Nb);
    # extract node data
    c₁ = rand(Uniform(1,3),Nb); c₂ = rand(Uniform(0.1, 0.3),Nb);
    d = rand(Uniform(0.5,1),Nb); p̅ = p̅ .+ 3
    # extract edge data
    for l in 1:Nl
        s_end[l] = net_data["branch"][string(l)]["f_bus"]
        r_end[l] = net_data["branch"][string(l)]["t_bus"]
        β[l] = -imag(1/(net_data["branch"][string(l)]["br_r"] + net_data["branch"][string(l)]["br_x"]im))
        f̅[l] = net_data["branch"][string(l)]["rate_a"]
    end
    # extract weighted Laplace matrix
    B = zeros(Nb,Nb); Z_eq = zeros(Complex{Float64},Nb,Nb)
    for l in 1:Nl
        Z_eq[net_data["branch"][string(l)]["f_bus"],net_data["branch"][string(l)]["t_bus"]] == 0 ? Z_eq[net_data["branch"][string(l)]["f_bus"],net_data["branch"][string(l)]["t_bus"]] = net_data["branch"][string(l)]["br_r"] + net_data["branch"][string(l)]["br_x"]im : Z_eq[net_data["branch"][string(l)]["f_bus"],net_data["branch"][string(l)]["t_bus"]] = (Z_eq[net_data["branch"][string(l)]["f_bus"],net_data["branch"][string(l)]["t_bus"]]*(net_data["branch"][string(l)]["br_r"] + net_data["branch"][string(l)]["br_x"]im))/(Z_eq[net_data["branch"][string(l)]["f_bus"],net_data["branch"][string(l)]["t_bus"]]+(net_data["branch"][string(l)]["br_r"] + net_data["branch"][string(l)]["br_x"]im))
        Z_eq[net_data["branch"][string(l)]["t_bus"],net_data["branch"][string(l)]["f_bus"]] == 0 ? Z_eq[net_data["branch"][string(l)]["t_bus"],net_data["branch"][string(l)]["f_bus"]] = net_data["branch"][string(l)]["br_r"] + net_data["branch"][string(l)]["br_x"]im : Z_eq[net_data["branch"][string(l)]["t_bus"],net_data["branch"][string(l)]["f_bus"]] = (Z_eq[net_data["branch"][string(l)]["t_bus"],net_data["branch"][string(l)]["f_bus"]]*(net_data["branch"][string(l)]["br_r"] + net_data["branch"][string(l)]["br_x"]im))/(Z_eq[net_data["branch"][string(l)]["t_bus"],net_data["branch"][string(l)]["f_bus"]]+(net_data["branch"][string(l)]["br_r"] + net_data["branch"][string(l)]["br_x"]im))
    end
    for i in 1:Nb, j in 1:Nb
        Z_eq[i,j] != 0 ? B[i,j] = imag(1/Z_eq[i,j]) : 0
    end
    for i in 1:Nb
        val = -sum(B[i,:])
        B[i,i] = val
    end
    # choose reference node
    refnode = Nb
    # save network data to dictionary
    data = Dict(:p̅ => p̅, :p̲ => p̲, :d => d, :c1 => c₁, :c2 => c₂,
                :s_end => s_end, :r_end => r_end, :β => β, :f̅ => f̅,
                :B => B, :refnode => refnode, :Nb => Nb, :Nl => Nl)
    return data
end
