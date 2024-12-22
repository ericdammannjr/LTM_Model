using OrdinaryDiffEq,Plots
pyplot()

include("fns/events.jl")
include("fns/lagoon_bed_shear_stress.jl")
include("fns/marsh_bed_shear_stress.jl")
include("fns/system_equations.jl")
include("fns/wave_height.jl")
include("fns/wave_number.jl")
include("fns/wave_period.jl")
include("fns/wave_power.jl")

# Computational Parameters

t0 = 0
tf = 1000
tspan = (t0,tf)
n = 100000
h = (tf-t0)/n

# Sea Level Rise Rate

zdot = 0.005

# Barrier Dynamics Parameters

αe = 0.02
bbmc = 1000
Dt = 10
He = 2
K = 2000
Qow_max = 100
We = 800
Vd_max = He*We

# Marsh-Lagoon Dynamics Parameters

β = 10^-3
Bpeak = 2.5
χref = 0.158
Co = 0.03
Dmin = 0 
g = 9.80171
ka = 2 
ke = 0.15
ko = 0.001
λ = 0.0001
νGp = 0.0138
P = 12.5/(24*365)
por = 1000/2650
r = 1.4
ρ = 1000
ρo = 1000
τcr = 0.1
U = 10
ws = 0.5*10^-3*(60*60*24*365)
x = 10
Dmax = 0.7167*r-0.0483

# Sea Level Initial Condition

Z0 = Dt

# Barrier Island Intital Conditions

xt0 = 0
xs0 = Dt/αe
H0 = He
xb0 = Dt/αe+We

# Marsh-Lagoon Initial Conditions

bbm0 = 1000
bL0 = 10000
bim0 = 2000
zm0 = Dmax/2
zL0 = 2
xbm0 = xb0+bbm0
xim0 = xb0+bbm0+bL0
xmm0 = xb0+bbm0+bL0+bim0-(zm0-r/2)/β

# Solution to System of ODEs

p = [zdot αe bbmc Dt He K Qow_max We Vd_max β Bpeak χref Co Dmin g ka ke ko λ νGp P por r ρ ρo τcr U ws x]
u0 = [Z0 xt0 xs0 H0 xb0 xbm0 xim0 xmm0 zm0 zL0]
prob = ODEProblem(system_equations,u0,tspan,p)
sol = solve(prob,Tsit5(),adaptive = false,dt = h,callback = events()) 

# Solution Vectors
    
Z = sol[1,:]
xt = sol[2,:]
xs = sol[3,:]
H = sol[4,:]
xb = sol[5,:]
xbm = sol[6,:]
xim = sol[7,:]
xmm = sol[8,:]
zm = sol[9,:]
zL = sol[10,:]

# System Evolution Animation

anim = @animate for k in 1:round(Int,n/100):length(Z)
    
    plot(
        
        title = "t = $(round(Int,sol.t[k])) years",
        label = "HH",
        xlims = (0,(xmm[length(Z)]/10^3+1)),
        xlabel = "Kilometers",
        ylims = (0,(Z[length(Z)]+H[length(Z)]+5)),
        ylabel = "Meters",
        framestyle = :box,
        legend = :bottomright,
        legend_columns = 1,
        legendfontsize = 8,
        foreground_color_legend = nothing,
        background_color_legend = nothing,
        grid = false,

    )

    # Time Scale

    annotate!(0.025*((xmm[length(Z)]+1000)/10^3), 0.9375*(Z[length(Z)]+H[length(Z)]+5), 
    
        text("1 second = $((tf-t0)/10) years", :left, 8)
    
    )

    # Time Step Size

    annotate!(0.975*((xmm[length(Z)]+1000)/10^3), 0.9375*(Z[length(Z)]+H[length(Z)]+5), 
    
        text("dt = $h years", :right, 8)

    )

    # Sea Level 

    plot!([0,xs[k]]/10^3,[Z[k],Z[k]],

        label = "Mean Sea Level",
        linestyle = :solid,
        linewidth = 1,
        linecolor = :lightblue
    
    )

    plot!([xb[k],xmm[k]+zm[k]/β]/10^3,[Z[k]+r/2,Z[k]+r/2],
    
        label = false,
        linestyle = :dash,
        linewidth = 1,
        linecolor = :lightblue
        
    )

    plot!([xbm[k],xim[k]]/10^3,[Z[k]-r/2,Z[k]-r/2],
    
    label = false,
    linestyle = :dash,
    linewidth = 1,
    linecolor = :lightblue
    
    )

    # Shoreface Toe Migration

    plot!(xt[1:k]./10^3,(Z[1:k].-Dt),
    
        label = "Trajectory of Shoreface Toe",
        linestyle = :dash,
        linewidth = 1,
        linecolor = :black
        
    )

    # Barrier Island

    plot!([xt[k],xs[k],xs[k],xb[k],xb[k]]/10^3,[Z[k]-Dt,Z[k],Z[k]+H[k],Z[k]+H[k],Z[k]+r/2-zm[k]],
    
        label = false,
        linestyle = :solid,
        linewidth = 1,
        linecolor = :tan
        
    )

    # Marshes

    plot!([xb[k],xbm[k]]/10^3,[Z[k]+r/2-zm[k],Z[k]+r/2-zm[k]],
    
        label = "Marshes",
        linestyle = :solid,
        linewidth = 1,
        linecolor = :green
        
    )

    plot!([xim[k],xmm[k]]/10^3,[Z[k]+r/2-zm[k],Z[k]+r/2-zm[k]],
    
        label = false,
        linestyle = :solid,
        linewidth = 1,
        linecolor = :green
        
    )

    # Lagoon

    plot!([xbm[k],xbm[k],xim[k],xim[k]]/10^3,[Z[k]+r/2-zm[k],Z[k]+r/2-zL[k],Z[k]+r/2-zL[k],Z[k]+r/2-zm[k]],
    
        label = false,
        linestyle = :solid,
        linewidth = 1,
        linecolor = :tan
        
    )

    # Mainland

    plot!([xmm[k],xmm[length(Z)]+1000]/10^3,[Z[k]+r/2-zm[k],(Z[length(Z)]+r/2-zm[length(Z)])+β*(1000)],
    
        label = false,
        linestyle = :solid,
        linewidth = 1,
        linecolor = :black

    )

end

fig1 = gif(anim, "figs/fig1.gif", fps = 10)

display(fig1)

# System Evolution Plots

plt1 = plot(sol.t,xbm-xb,

    title = "Back Barrier Marsh Width",
    xlabel = "Years",
    ylabel = "Meters",
    legend = false,
    grid = false,
    linestyle = :solid,
    linewidth = 1,
    linecolor = :black,

)

plt2 = plot(sol.t,(xim-xbm)/10^3,

    title = "Lagoon Width",
    xlabel = "Years",
    ylabel = "Kilometers",
    legend = false,
    grid = false,
    linestyle = :solid,
    linewidth = 1,
    linecolor = :black,

)

plt3 = plot(sol.t,xb-xs,

    title = "Barrier Width",
    xlabel = "Years",
    ylabel = "Meters",
    legend = false,
    grid = false,
    linestyle = :solid,
    linewidth = 1,
    linecolor = :black,

)

plt4 = plot(sol.t,zL,

    title = "Lagoon Depth",
    xlabel = "Years",
    ylabel = "Meters",
    legend = false,
    grid = false,
    linestyle = :solid,
    linewidth = 1,
    linecolor = :black,

)

fig2 = plot(plt1,plt2,plt3,plt4,layout=(2,2))

savefig(fig2,"figs/fig2.png")

display(fig2)