#Faut ecrire un commentaire pour cette partie.
    """
    Methode de tir multiple avec p0,t1,t2 comme inconnus
    avec  des conditions de transversalité et condition sur la fonction de commutation
    """

    import Pkg;
    # Packages
    using DifferentialEquations
    using Plots
    using NLsolve
    using ForwardDiff
    using LinearAlgebra
    include("utils.jl"); # plot_traj!, plot_flow!, Flow

    #Hamiltonien du systeme
    tf = 40.0
    t0 = 0.0
    si1 = 1
    si2 = 1
    a1  = 1
    a2  = 1
    alpha1 = 1
    alpha2 = 1
    u_bar = 0.62

    function hv(m,p, u)
        m1,m2 = m[1],m[2]
        p1,p2 = p[1],p[2]
        psi = 1/(1+m1/si1 + m2/si2)
        dm1 = alpha1*psi*(1+u)/2 - a1*m1*(1-u)/2
        dm2 = alpha2*psi*(1+u)/2 - a2*m2*(1-u)/2
        dp1 = 1/si1*(u + p1*alpha1*(1+u)/2 + p2*alpha2*(1+u)/2)*psi^2 + (1-u)/2*p1*a1
        dp2 =  1/si2*(u + p1*alpha1*(1+u)/2 + p2*alpha2*(1+u)/2)*psi^2 + (1-u)/2*p2*a2
        return [dm1,dm2],[dp1,dp2]
    end
    # systèmes hamiltoniens
    hv_min(m,p) = hv(m,p, -1) ## Valeur u = -1
    hv_max(m,p) = hv(m,p, 1) ## Valeur u = 1
    hv_sin(m,p) = hv(m,p, u_bar) ## Valeur u = u_bar 

    # flots
    fmin = Flow(hv_min)
    fmax = Flow(hv_max)
    fsin = Flow(hv_sin)
    m0 = [0.0,0.0]

    # fonction de tir
    function shoot(params)

        k1,k2,t1,t2 = params[1:4]
        p0 = [k1,k2]
        # Integration
        m1, p1 = fmax(t0, m0, p0, t1)
        m2, p2 = fsin(t1, m1, p1, t2)
        m3, p3 = fmax(t2, m2, p2, tf)

        return [p3[1], p3[2],0.0,0.0]
    end
    # itéré initiale pour la méthode indirecte de tir
    y =  [-0.4085;-0.4085;3.168;35.64]
    println("Itéré initial:\n", y)

    # fonction de tir et sa jacobienne
    foo(y) = shoot(y...)
    jfoo(y) = ForwardDiff.jacobian(shoot, y)

    # Résolution de shoot(p0, t1, t2) = 0.
    nl_sol1 = NLsolve.nlsolve(foo, jfoo, y; xtol=1e-6, method=:trust_region, show_trace=true);

    # Retrieves solution
    if converged(nl_sol1)
        p01 = nl_sol1.zero[1]
        p02 = nl_sol1.zero[2]
        t1 = nl_sol1.zero[3]
        t2 = nl_sol1.zero[4]

        println("\nFinal solution:\n", nl_sol1.zero)
    else
        error("Not converged")
    end

    #Plots 
    # affichage de la solution
    p0 = [p01,p02]

    ode_sol = fmax((t0, t1), m0, p0)
    tt0 = ode_sol.t
    xx0 = [ ode_sol[1:2, j] for j in 1:size(tt0, 1) ]
    pp0 = [ ode_sol[3:4, j] for j in 1:size(tt0, 1) ]
    uu0 = ones(size(tt0, 1))


    ode_sol = fsin((t1, t2), xx0[end], pp0[end])
    tt1 = ode_sol.t
    xx1 = [ ode_sol[1:2, j] for j in 1:size(tt1, 1) ]
    pp1 = [ ode_sol[3:4, j] for j in 1:size(tt1, 1) ]
    uu1 = zeros(size(tt1, 1)) .+ u_bar

    ode_sol = fmax((t2, tf), xx1[end], pp1[end])
    tt2 = ode_sol.t
    xx2 = [ ode_sol[1:2, j] for j in 1:size(tt2, 1) ]
    pp2 = [ ode_sol[3:4, j] for j in 1:size(tt2, 1) ]
    uu2 = ones(size(tt2, 1))

    t_shoot = [ tt0 ; tt1 ; tt2 ]
    x_shoot = [ xx0 ; xx1 ; xx2 ]
    p_shoot = [ pp0 ; pp1 ; pp2 ]
    u_shoot = [ uu0 ; uu1 ; uu2 ]    


    x_1 = [e[1] for e in x_shoot[1:end] ]
    x_2 = [e[2] for e in x_shoot[1:end] ]
    p_1 = [e[1] for e in p_shoot[1:end] ]
    p_2 = [e[2] for e in p_shoot[1:end] ]

    #Tracer les solutions
    x_plot_1 = plot(t_shoot, x_1, xlabel = "t", ylabel = "m1", legend = false)
    x_plot_2 = plot(t_shoot, x_2, xlabel = "t", ylabel = "m2", legend = false)
    px_plot_1 = plot(t_shoot, p_1, xlabel = "t", ylabel = "px1", legend = false)
    px_plot_2 = plot(t_shoot, p_2, xlabel = "t", ylabel = "px2", legend = false)
    u_plot = plot(t_shoot, u_shoot, xlabel = "t", ylabel = "u", legend = false, size=(800,300), linetype=:steppre)

    #Affichage des solutions
    display(plot(x_plot_1,px_plot_1,x_plot_2,px_plot_2, layout = (2,2), size=(800, 300)))
    display(u_plot)







