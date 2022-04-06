using LinearAlgebra, Printf
using MatrixEquations 
using LightGraphs, SimpleWeightedGraphs
using Plots, PGFPlotsX, LaTeXStrings
using PowerModels
pgfplotsx()
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{amsmath}")
PowerModels.silence()

function truncate(K,g,Is,Js,kappa=0)
    for i in 1:nv(g)
        d = gdistances(g,i)
        for j in 1:nv(g)
            if d[j] > kappa
                for ii in Is[i]
                    for jj in Js[j]
                        K[jj,ii] = 0.
                    end
                end
            end
        end
    end
    
    return K
end

function simple_2d_graph(n)
    
    g = SimpleWeightedGraph(n)
    c = .05
    
    for i=1:n
        for j=1:n
            if i!= n
                add_edge!(g,n*(i-1)+j,n*i+j,c)
            end
            if j!= n
                add_edge!(g,n*(i-1)+j,n*(i-1)+j+1,c)
            end
        end
    end
    
    return g
end

function dc_network()

    data = parse_file("pglib_opf_case118_ieee.m")
    c = (132e3)^2/1e5 * (1e-9) # (voltage^2/inertia) * (1 s/1e3 ms)^3

    n = length(data["bus"])
    g = SimpleWeightedGraph(n)
    
    for (~,branch) in data["branch"]
        i = branch["f_bus"]
        j = branch["t_bus"]
        ~, b = calc_branch_y(branch)

        add_edge!(g,i,j,abs(b)*c)
    end

    return g
end

function get_block_norm(K,g,uloc)
    D = maximum(maximum(gdistances(g,i)) for i in vertices(g))
    block_norm = zeros(D+1)

    for ii in eachindex(uloc)
        i = uloc[ii]
        d = gdistances(g,i)
        for j in 1:size(K,2)
            block_norm[d[j]+1] = max(block_norm[d[j]], K[ii,j])
        end
    end
    
    return block_norm
end
