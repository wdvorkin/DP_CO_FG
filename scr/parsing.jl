function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--case", "-c"
            help = "chose case: case3_lmbd case5_pjm case14_ieee case39_epri case57_ieee case118_ieee"
            arg_type = String
            default = "case3_lmbd"
        "--query", "-q"
            help = "chose query: identity or linear"
            arg_type = String
            default = "identity"
        "--reformulation", "-r"
            help = "chose reformulation: sample or analytic"
            arg_type = String
            default = "analytic"
        "--cost_fun", "-f"
            help = "chose cost function: linear or quadratic"
            arg_type = String
            default = "quadratic"
        "--adjacency", "-a"
            help = "chose adjacency α ∈ (0,1]"
            arg_type = Float64
            default = 0.1
        "--epsilon", "-e"
            help = "chose privacy loss ε ∈ (0,1]"
            arg_type = Float64
            default = 1.0
        "--violation", "-v"
            help = "chose constraint violation tolerance η ∈ (0,1)"
            arg_type = Float64
            default = 0.025
        "--share_node", "-s"
            help = "chose release node share γ ∈ [0,1]"
            arg_type = Float64
            default = 0.3
        "--exp_num", "-n"
            help = "chose number of random network samples (experiments)"
            arg_type = Int64
            default = 3
        "--cluster_num", "-k"
            help = "chose number of clusters for linear queries"
            arg_type = Int64
            default = 9
        "--obj_var_control", "-o"
            help = "obj variance control: true or false"
            arg_type = Any
            default = false
        "--sol_var_control", "-y"
            help = "solution variance control: true or false"
            arg_type = Any
            default = false
        "--var_penalty", "-p"
            help = "variance penalty factor φ"
            arg_type = Float64
            default = 0.99
    end
    return parse_args(s)
end
