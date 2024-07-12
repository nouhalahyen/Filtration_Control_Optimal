


#Fonction de commutation
function commut(m)
    return -2 / (b + a * m)
end

#Hamiltonien du systeme
function hv(m, p, u)
    g = 1 / (e + m)
    f2 = a * m / (e + m)
    f1 = b / (e + m)
    
    dm = (1 + u) / 2 * f1 - (1 - u) / 2 * f2
    dp = (p * (a * e + b) + u * (2 - p * (a * e - b))) / (2 * ((e + m)^2))
    return [dm, dp]
end
# systèmes hamiltoniens
hv_min(m, p) = hv(m, p, -1) ## Valeur u = -1
hv_max(m, p) = hv(m, p, 1) ## Valeur u = 1
hv_sin(m, p) = hv(m, p, u_bar) ## Valeur u = u_bar 

# flots
fmin = Flow(hv_min)
fmax = Flow(hv_max)
fsin = Flow(hv_sin)


#Fonction du systeme (dynamique)
function closed_loop!(dw, w, param, t)
    m, p = w[1], w[2]
    u_val = 1.0  # Valeur par défaut pour u
    dw[1], dw[2] = hv(m, p, u_val)
end

# Callback pour s'assurer que m reste positif et reste nulle une fois qu'elle s'annule
function condition(u, t, integrator)
    return u[1]  # Condition sur m
end

#Pour les call backs : arreteer l'integration
 function affect!(integrator)
    terminate!(integrator) # Forcer m à être au moins 0
 end

 function condition_commutation(u,t,integrator)
     return u[2] - commut(u[1])
 end 

 function affect_commutation!(integrator)
     global masse
     global temps
     global adjoint
     push!(masse, integrator.u[1])
     push!(temps, integrator.t)
     push!(adjoint, integrator.u[2])
     terminate!(integrator)
 end
 cb1 = ContinuousCallback(condition, affect!)
 cb2 = ContinuousCallback(condition_commutation, affect_commutation!)

 #Trcer les trajectoires issues de la methode directe
function plt_traj_directe!(plt_direct, times, Problem)

    for i =0.0 : 0.2:2.99
        direct_sol = Problem(i)
        xplot = [direct_sol.state(t) for t in times]
        plot!(times,xplot , size=(800, 800),label=false,lw = 2)
    end

    direct_sol = Problem(3.0)
    xplot = [direct_sol.state(t) for t in times]
    plot!(times,xplot , size=(800, 800),label=false,lw = 4,color =:green)

    for i =3.1 : 0.3:21.0
        direct_sol = Problem(i)
        xplot = [direct_sol.state(t) for t in times]
        plot!(times,xplot , size=(800, 800),label=false,lw = 2)
    end
    return nothing
end 
#Pour le tracé de la trajectoire dans la methode d'integration retrograde
function plt_traj_rev!(plt_direct, plt_adj, final_masses, pf=0.0 )
    for mf in final_masses
        wf = [mf, pf]
        prob = ODEProblem(closed_loop!, wf, tspan)
        sol = DifferentialEquations.solve(prob, Tsit5(), callback = CallbackSet(cb1,cb2 ),abstol=1e-10, reltol=1e-10)
        
        # Inverser les valeurs de temps et de l'état pour le tracé
        times_reversed = reverse(sol.t)
        m_values_reversed = reverse(sol[1, :])
        p_values_reversed = reverse(sol[2, :])
        
        # Visualiser les résultats avec l'axe des temps inversé
        plot!(plt_direct, times_reversed, m_values_reversed, label="", lw=2)
        plot!(plt_adj, times_reversed, p_values_reversed, label="", lw=2)
        scatter!(plt_direct, temps, masse, label="", color=:gray, markerstrokecolor=:gray)
        scatter!(plt_adj, temps, adjoint, label="", color=:gray, markerstrokecolor=:gray)

    end
end

function plt_traj_tir(x_plot, u_plot, px_plot,t0,t1,t2,m0,p0)

    ode_sol = fmin((t0, t1), m0, p0)

    tt0 = ode_sol.t
    xx0 = [ ode_sol[1, j] for j in 1:size(tt0, 1) ]
    pp0 = [ ode_sol[2, j] for j in 1:size(tt0, 1) ]
    uu0 = ones(size(tt0, 1))


    ode_sol = fsin((t1, t2), xx0[end], pp0[end])
    tt1 = ode_sol.t
    xx1 = [ ode_sol[1, j] for j in 1:size(tt1, 1) ]
    pp1 = [ ode_sol[2, j] for j in 1:size(tt1, 1) ]
    uu1 = zeros(size(tt1, 1)) .+ u_bar

    ode_sol = fmax((t2, tf), xx1[end], pp1[end])
    tt2 = ode_sol.t
    xx2 = [ ode_sol[1, j] for j in 1:size(tt2, 1) ]
    pp2 = [ ode_sol[2, j] for j in 1:size(tt2, 1) ]
    uu2 = ones(size(tt2, 1))

    t_shoot = [ tt0 ; tt1 ; tt2 ]
    x_shoot = [ xx0 ; xx1 ; xx2 ]
    p_shoot = [ pp0 ; pp1 ; pp2 ]
    u_shoot = [ uu0 ; uu1 ; uu2 ]    

    plot!(x_plot, t_shoot, x_shoot, xlabel = "t", ylabel = "m", legend = false)
    plot!(u_plot, t_shoot, u_shoot, xlabel = "t", ylabel = "u", legend = false, size=(800,300), linetype=:steppre)
    plot!(px_plot, t_shoot, p_shoot, xlabel = "t", ylabel = "px", legend = false)
    
end 