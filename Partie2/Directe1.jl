
    #import Pkg
    #Pkg.activate(".")
    #Pkg.update();
    #Pkg.add("NLPModelsIpopt")
    using NLPModelsIpopt
    using OptimalControl
    using Plots
    using JuMP
    using Ipopt
    using DifferentialEquations
    using QuadGK  # For numerical integration
    include("Optimum.jl") #For the optimum solutions 

    t0 = 0     
    tf = 40.0

    α1 = 1
    α2 = 1.5

    a1 = 2.0
    a2 = 1.0

    # Définition des fonctions dynamiques
    O(x1, x2) = 1 / (1 + x1 + x2 )
    F11(x1, x2) = α1 * O(x1, x2)
    F21(x1, x2) = α2 * O(x1, x2)
    F12(x) = a1 * x
    F22(x) = a2 * x

    # Définir le problème de contrôle optimal
    Problem(x0) = begin 
        @def ocp begin
            t ∈ [t0, tf], time
            x ∈ R², state
            u ∈ R, control
            m1 = x₁
            m2 = x₂

            x(t0) == x0
            -m1(t) ≤ 0.0
            -m2(t) ≤ 0.0
            -1.0 ≤ u(t) ≤ 1.0
            
            ẋ(t) == [(1 + u(t)) / 2 * F11(m1(t), m2(t)) - (1 - u(t)) / 2 * F12(m1(t)),
                    (1 + u(t)) / 2 * F21(m1(t), m2(t)) - (1 - u(t)) / 2 * F22(m2(t))]
            
            ∫(u(t) * O(m1(t), m2(t))) → max
        end
    end
    #Optimal solution: m1 = 1.2492350871367335, m2 = 0.8328233914244892, u = 0.660476461880673
    x0 = [0.7,0.0]
    # Résoudre le problème de contrôle optimal
    prob = Problem(x0)
    sol_Directe = OptimalControl.solve(prob,grid_size=4000)
    #plotlyjs()
    #plot(sol_Directe, title="  α1=$α1, α2=$α2, a1 = $a1, a2 = $a2", size=(950, 600))  # Décommenter pour afficher le graphique de la solution

    ##########################################################################################
    ##########################################################################################
    
    #Code pour trouver les t1 et t2 optimaux à partir de "Optimum.jl"

    #Solution optimale : 
    u_bar = 0.6636975627412482

    res = optimum(x0, u_bar)
    # Get the results
    t1_val, t2_val, objective_val = res

    println("Optimal solution: t1 = $t1_val, t2 = $t2_val")
    println("Objective value: $objective_val")

    
    ###############################################################################################
    ###############################################################################################

    #Tracé de la solution sous-optimale et directe  et comparaison des résultats: 

    #u0 = [0.0, 0.0]  # Initial conditions for m1, m2
    tspan = (0.0, tf)
    p = [t1_val, t2_val]  # Parameters to pass t1 and t2

    prob = ODEProblem(my_dynamics!, x0, tspan, p)
    sol_ODE = solve(prob, AutoVern7(Rodas5()), saveat=0.01)  # Solve the ODE

    plt_m1 = plot()
    plt_m2 = plot()
    axe = 0:(40/4000):40
    axe = collect(axe)
    time = sol_Directe.time_grid
    tableau_m1 = [sol_Directe.state(t)[1] for t in time]
    tableau_m2 = [sol_Directe.state(t)[2] for t in time]

    plot!(plt_m1, axe,tableau_m1,label="m1 Directe")
    plot!(plt_m1,axe, sol_ODE[1,:], label = "m1 ODE")

    plot!(plt_m2, sol_ODE[2,:], label="m2 ODE")
    plot!(plt_m2, tableau_m2, label="m2 Directe")

    display(plt_m1)
    display(plt_m2)

    #Difference entre solution exacte et sous-optimale
    plt = plot()
    diff = abs.((tableau_m1 - sol_ODE[1,:]).^2 ./tableau_m1) .+ abs.((tableau_m2 - sol_ODE[2,:]).^2 ./tableau_m2)   
    plot(plt , axe, diff)

    #Difference entre m_bar et solution exacte : 
    #Optimal solution to the static problem: m1_bar = 0.6712433647609589, m2_bar = 2.0137300942828773, u_bar = 0.6636975627412119
    m1_bar = 0.6712433647609589
    m2_bar = 2.0137300942828773
    plotlyjs()
    tableau_m1_bar = [m1_bar for t in time]
    tableau_m2_bar = [m2_bar for t in time]
    diff_m = (tableau_m1 - tableau_m1_bar).^2 + (tableau_m2 - tableau_m2_bar).^2
    plot(axe, diff_m, label="diff_m_bar")

    ##########################################################################################
    ##########################################################################################
    # Evolution of the objective function in terms of t1 chosen for a certain m0
    
    times_t1 = 0.0:0.0001:5.0
    obj = [integrate_objective(x0,t1,t2_val) for t1 in times_t1]
    plotlyjs()
    plot(times_t1,obj, title = "Evolution de la fonction objectif en fonction de t1",xlabel = "t1",ylabel = "Cout",label = false)