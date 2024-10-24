
using Plots

t0 = 0     
tf = 40.0

σ1 = 1
σ2 = 1

α1 = 1
α2 = 1

# Dynamics
a1 = 1
a2 = 1

F11(x1,x2) = begin
    return  (α1*O(x1,x2))
end
F21(x1,x2) = begin  
    return  (α2*O(x1,x2))
end
F12(x) = begin
    return  (a1*x) 
end;
F22(x) = begin
    return  (a2*x) 
end;

O(x1,x2) = begin 
    return 1/(1 + x1/σ1 +x2/σ2)
end;

function commut(m1,m2,p1,p2)
    return (1 + 1/2*(p1*α1 + p2*α2))*O(m1,m2) + 0.5*(p1*a1*m1 + p2*a2*m2)
end

# Function decrivant la dynamique du modele
function dynamique(m1,m2,p1,p2)

    if commut(m1,m2,p1,p2)>0
        u = 1.0
    elseif  commut(m1,m2,p1,p2)<0 
        u = -1.0
    end

    dm1 = (1+u)/2 * F11(m1,m2) - (1 - u)/2 * F12(m1)
    dm2 =  (1+u)/2 * F21(m1,m2) - (1 - u)/2 * F22(m2)
    dp1 = 1/σ1 * (u + (p1*α1 + p2*α2)*(1+u)/2)*O(m1,m2)^2 + (1-u)/2*p1*a1
    dp2 = 1/σ2 * (u + (p1*α1 + p2*α2)*(1+u)/2)*O(m1,m2)^2 + (1-u)/2*p2*a2

    return dm1,dm2,dp1, dp2, u
end

# Grid of points for m and p
m_values1 = (0.0:0.1:1.0)
m_values2 = (0.0:0.1:1.0)
p_values1 = (-3.0:0.5:0.0)
p_values2 = (-3.0:0.5:0.0)

#Initialize arrays for vector components and colors

#For the masses
DM11, DM21 = [],[]
DM12, DM22 = [],[]

#Pour u = 1
X11, Y11= [], []
#Pour u = -1
X12, Y12= [], []


# Calculate the vector field and separate based on the value of u
for m1 in m_values1
    for m2 in m_values2
        for p1 in p_values1
            for p2 in p_values2
                dm1,dm2, u = dynamique(m1,m2,p1,p2)    
                if u == 1.0
                    push!(DM11, dm1)#m1
                    push!(DM21, dm2)
                    push!(X11, m1) #m1
                    push!(Y11, m2)
                                    
                else
                    push!(DM12, dm1)
                    push!(DM22, dm2)
                    push!(X12, m1)                    
                    push!(Y12, m2)

                end
            end
        end
    end
end      

# Plot the vector fields

plotlyjs()
#plt2 = plot()
#quiver!(X21, Y21,Y22, quiver=(0.02.*DM21, 0.02.*DP21,0.02.*DP22), color=:blue, label="u = -1")
quiver!( X11, Y11,quiver=(DM11, DM21), color=:red, label="u = 1")
quiver!( X12, Y12,quiver=(DM12,DM22), color=:blue, label="u = -1")

xlabel!("m1")
ylabel!("m2")
#quiver!(plt2,X22, Y22, quiver=(0.2.*DM22, 0.2.*DP22), color=:blue, label="u = -1")
#quiver!(plt2,X12, Y12, quiver=(0.2.*DM12, 0.2.*DP12), color=:red, label="u = 1")

