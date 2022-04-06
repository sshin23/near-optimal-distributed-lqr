include("include.jl")

gs = [simple_2d_graph(10), dc_network()]
names = ["2d-mesh","dc-118ieee"]
eta_press= [-1:-1:-4, -1:-1:-4]
dts = [1 .005] # 1 sec / .005 ms

for (g,name,eta_pres,dt) in zip(gs,names,eta_press,dts)

    println(name)
    n = nv(g)
    L = laplacian_matrix(g);

    A =  if name == "2d-mesh"
        # heat
        Matrix( [I dt*I; 0*I I-dt*L])
    elseif name == "dc-118ieee"
        # power
        Matrix( [I dt*I; -dt*L I])
    end
    
    R = I

    Is = [[i,n+i] for i=1:n]
    Js = [[i] for i=1:n]
    Ks = [[i] for i=1:n]

    kappas = 0:5

    
    regrets = []
    sradiuss = []
    errors = []
    
    etas = 2 .^Float64.(eta_pres)
    lbls = permutedims([L"\eta = 2^{%$eta_pre}" for eta_pre in eta_pres])
    Ks = Matrix[]

    for eta in etas
        println("eta = $eta")

        B = [
            0*I;
            dt*diagm(
                0 => eta .* ones(n)
            )
        ]

        C = [
            eta*I diagm(
                0 => 0*ones(n)
            ) 
        ]

        Q = C'*C


        P,~ = ared(A,B,R,Q,rtol=1e-12)
        K = inv(R+B'*P*B)*B'*P*A

        regret = Float64[]
        sradius= Float64[]
        error= Float64[]
        
        for kappa = kappas
            _K = truncate(copy(K),g,Is,Js,kappa)
            _P = lyapd((A-B*_K)',1.,Q + _K'*R*_K)
            s = maximum(abs.(eigen(A-B*_K).values))
            push!(sradius,s)
            if s < 1
                push!(regret,(opnorm(_P-P))/opnorm(P))
            else
                push!(regret,Inf)
            end
            push!(error,opnorm(K-_K)/opnorm(K))
        end
        
        push!(errors,error)
        push!(sradiuss,sradius)
        push!(regrets,regret)
        push!(Ks,K)
    end

    plt1 = plot(
        kappas, hcat(regrets...);
        yscale = :log10,
        markershape = :auto,
        framestyle = :box,
        ylabel = L"\|\boldsymbol{P}^\kappa-\boldsymbol{P}^\star\|/\|\boldsymbol{P}^\star\|",
        xlabel = L"\kappa",
        labels = lbls,
        legend = :bottomleft,
        size = (450,300),
    )

    plt2 = plot(
        kappas, hcat(errors...);
        yscale = :log10,
        markershape = :auto,
        framestyle = :box,
        ylabel = L"\|\boldsymbol{K}^\kappa-\boldsymbol{K}^\star\|/\|\boldsymbol{K}^\star\|",
        xlabel = L"\kappa",
        labels = lbls,
        legend = :bottomleft,
        size = (450,300), 
    )
    
    savefig(plt1,"$name-P.pdf")
    savefig(plt2,"$name-K.pdf")
end
