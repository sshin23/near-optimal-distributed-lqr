# Near-Optimal Distributed Linear-Quadratic Regulator for Networked Systems

Code for "Near-Optimal Distributed Linear-Quadratic Regulator for Networked Systems" by  Sungho Shin, Yiheng Lin, Guannan Qu, Adam Wierman, and Mihai Anitescu.

# How to run
1. Install [Julia v1.7.2](https://julialang.org/downloads/)
2. Clone the repository and run the script
```
git clone https://github.com/sshin23/near-optimal-distributed-lqr
cd near-optimal-distributed-lqr
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. example.jl
```
