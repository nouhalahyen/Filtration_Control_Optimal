#########################################################################################################################################
# Troisiéme modéle d'étude visant à voir l'effet de a2 - a1 et alpha2 - alpha1 :

using OptimalControl
using Plots

t0 = 0     
tf = 40.0
#parametres du modele 1:
σ1 = 1
σ2 = 1

#parametres du modele 2 :

γ1 = 1/σ1
γ2 = 1/σ2



# Définition des fonctions dynamiques
O(x) = 1 / (1 + x)
F11(x) = b_bar * O(x)
F12(x,δ) = a_bar*x + eps_a*δ*c1 
F21(x) = eps_b * O(x)
F22(x,δ) = c2*δ - eps_a* x /(γ1 + γ2)

x0 = [0.0, 0.0]
# Définir le problème de contrôle optimal
Problem(x0,a1,a2,α1,α2) = begin 

    a_bar = (γ1*a1 + γ2*a2)/(γ1+γ2)
    b_bar = γ1*α1 + γ2*α2
    eps_a = a2 - a1
    eps_b = α2 - α1

    c1 = γ1*γ2/(γ1 + γ2)
    c2 = (γ1*a2 + γ2*a1)/(γ1 + γ2)

    @def ocp begin
        t ∈ [t0, tf], time
        x ∈ R², state
        u ∈ R, control
        m = x₁
        δ = x₂

        x(t0) == x0
        -m(t) ≤ 0.0
        
        -1.0 ≤ u(t) ≤ 1.0
        
        ẋ(t) == [(1 + u(t)) / 2 * F11(m(t)) - (1 - u(t)) / 2 * F12(m(t),δ(t)),
                 (1 + u(t)) / 2 * F21(m(t)) - (1 - u(t)) / 2 * F22(m(t),δ(t))]
        
        ∫(u(t) * O(m(t))) → max
    end
end
# Résoudre le problème de contrôle optimal
prob = Problem(x0,a1,a2,α1,α2)
sol = solve(prob,grid_size=1500)
plotlyjs()
plot(sol, title="σ1 = $σ1, σ2 = $σ2, α1=$α1, α2=$α2, a1 = $a1, a2 = $a2", size=(950, 600))  # Décommenter pour afficher le graphique de la solution



##########################################################
#Probleme quand a2=a1 et alpha2 = alpha1
F11(x) = b_bar * O(x)
F12(x,δ) = a_bar*x 
F21(x) = 0
F22(x,δ) = c2*δ 
x0 = [0.0, 0.0]
# Définir le problème de contrôle optimal
Problem(x0) = begin 
    @def ocp begin
        t ∈ [t0, tf], time
        x ∈ R², state
        u ∈ R, control
        m = x₁
        δ = x₂

        x(t0) == x0
        -m(t) ≤ 0.0
        
        -1.0 ≤ u(t) ≤ 1.0
        
        ẋ(t) == [(1 + u(t)) / 2 * F11(m(t)) - (1 - u(t)) / 2 * F12(m(t),δ(t)),
                 (1 + u(t)) / 2 * F21(m(t)) - (1 - u(t)) / 2 * F22(m(t),δ(t))]
        
        ∫(u(t) * O(m(t))) → max
    end
end

# Résoudre le problème de contrôle optimal
prob = Problem(x0)
sol = solve(prob,grid_size=1500)
plotlyjs()
plot(sol, title="σ1 = $σ1, σ2 = $σ2, α1=$α1, α2=$α2, a1 = $a1, a2 = $a2", size=(950, 600))  # Décommenter pour afficher le graphique de la solution
##############################################################
#Faire plusieurs courbes de solutions pour quand a2 s'éloigne de a1.