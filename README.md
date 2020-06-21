# DP_CO_FG

This repository contains the code for the paper __Differentially Private Convex Optimization with Feasibility Guarantees__ submitted to the 34th Conference on Neural Information Processing Systems

The optimization models were implemented in [Julia](https://juliacomputing.com/products/juliapro) (v.1.4) using [JuMP](https://github.com/JuliaOpt/JuMP.jl) modeling language for mathematical optimization embedded in Julia. The models require [Mosek](https://www.mosek.com) commercial optimization solver, which needs to be installed and licensed. 

To activate the packages in ```Project.toml```, ```cd``` to the project directory and call
```
$ julia 
julia> ]
(@v1.4) pkg> activate .
(code_compainion) pkg> instantiate
```

where ```julia``` is an alias to Julia installation. To run the code, ```cd``` to the project directory and call
```
$ julia main.jl
```

By default, the program returns the solution of the ```PIQ-a``` and ```PIQ-s``` algorithms for the ```case3_lmbd``` test case  for 100 simulation runs (Table 1). The results will be stored in ```output``` folder. To alter the default settings, specify the arguments in the command line. For example, 
```
$ julia main.jl -q "linear"
```
returns the same results for the ```PSQ-a``` and ```PSQ-s``` algorithms (Table 2). To access all available options, call 
```
$ julia main.jl --help
```

The simulations were carried out using the standard PC with Intel Core i5 3.4 GHz processor and 8 GB memory. Solving optimization problems with the analytic reformulation requires less than a few seconds on average, whereas the sample approximation of the chance constraints requires by at most 78 seconds for ```PIQ-s``` algorithm on average. 

