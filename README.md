# DP_CO_FG

This repository contains supplementary materials for the paper __Differentially Private Convex Optimization with Feasibility Guarantees__ by V. Dvorkin, F. Fioretto, P. Van Hentenryck, J. Kazempour, and P. Pinson.

The optimization models were implemented in [Julia](https://juliacomputing.com/products/juliapro) (v.1.4) using [JuMP](https://github.com/JuliaOpt/JuMP.jl) modeling language for mathematical optimization embedded in Julia. The models require [Mosek](https://www.mosek.com) comercial optimization solver, which needs to be installed and licensed. 

To activate the packages in ```Project.toml```, clone the project using e.g. ```git clone```, ```cd``` to the project directory and call
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

By default, the program returns the solution of the ```PIQ-a``` and ```PIQ-s``` algorothms for the ```case118_ieee``` and 100 siulation runs. The results will be stored in ```output``` folder. To alter the default settings, specify the arguments in the command line. For example, 
```
$ julia main.jl -q "linear"
```
returns the same results for the ```PSQ-a``` and ```PSQ-s``` algorithms. 

