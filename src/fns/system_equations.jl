function system_equations(du,u,p,t)

    # Rate of Sea Level Rise

    zdot = p[1]

    # Barrier Island Parameters

    αe = p[2]
    bbmc = p[3]
    Dt = p[4]
    He = p[5]
    K = p[6]
    Qow_max = p[7]
    We = p[8]
    Vd_max = p[9]

    # Marsh-Lagoon Parameters

    β = p[10]
    Bpeak = p[11]
    χref = p[12]
    Co = p[13]
    Dmin = p[14] 
    g = p[15]
    ka = p[16]
    ke = p[17]
    ko = p[18]
    λ = p[19]
    νGp = p[20]
    P = p[21]
    por = p[22]
    r = p[23]
    ρ = p[24]
    ρo = p[25]
    τcr = p[26]
    U = p[27]
    ws = p[28] 
    x = p[29]

    # Sea Level 

    Z = u[1]

    # Barrier Island State Variable

    xt = u[2]
    xs = u[3]
    H = u[4]
    xb = u[5]

    # Marsh-Lagoon State Variables

    xbm = u[6]
    xim = u[7]
    xmm = u[8]
    zm = u[9]
    zL = u[10]

    # Additional Variable Calculations

    α = Dt/(xs-xt)
    W = xb-xs
    bbm = xbm-xb
    bL = xim-xbm
    bim = xmm-xim
    Dmax = 0.7167*r-0.0483

    # Barrier Island

        # Deficit Volume Calculations 

        ϕ = min(1,bbm/bbmc)
        Vd_B = max(0,(We-W)*(H+ϕ*(zm-r/2)+(1-ϕ)*(zL-r/2)))
        Vd_H = max(0,(He-H)*(W))
        Vd = Vd_B+Vd_H

        # Overwash Calculations

        if Vd<Vd_max

            Qow_H = Qow_max*Vd_H/Vd_max
            Qow_B = Qow_max*Vd_B/Vd_max

        else 

            Qow_H = Qow_max*Vd_H/Vd
            Qow_B = Qow_max*Vd_B/Vd

        end

        Qow = Qow_H+Qow_B
        Qow_Bl = (1-ϕ)*Qow_B
        Qow_Bm = ϕ*Qow_B

        # Sediment Flux at the Shoreface
        
        Qsf = K*(αe-α)

    # Marsh-Lagoon

        # Lagoon 

        ZL = (zL+(zL-min(r,zL)))/2
        τL = lagoon_bed_shear_stress(bL,ZL,g,ko,U)
        SL = max((τL-τcr)/τcr,0)
        Cr = ρ*λ*SL/(1+SL*λ)
        Fc = (Cr-Co)*min(r,zL)/P/ρ

        # Marshes

        Zm = (zm+(zm-min(r,zm)))/2
        B = Bpeak*(Dmax-zm)*(zm-Dmin)/(0.25*(Dmax-Dmin)^2)

        if B <= 1*ℯ^-3

            B = 0

        end

        Bfrac = B/Bpeak
        AMC = 180*νGp*B
        Rref = AMC*χref
        O = 1/por*(Rref/ρo)

        if Zm > 1*ℯ^-4

            τm = marsh_bed_shear_stress(bL,Zm,Bfrac,g,ko,U)

        else

            τm = 0 

        end

        Sm = max((τm-τcr)/τcr,0)
        Cm = ρ*λ*Sm/(1+Sm*λ)
        Fm = (Cr-Cm)*min(r,zm)/P/ρ

        # Marsh-Lagoon Edges

        zs = zm+(zL-zm)*(1-exp(-x*0.1/zL))
        Zs = (zs+(zs-min(r,zs)))/2
        WP = wave_power(bL,Zs,g,U)
        E = ke*WP/(zs-zm)-ka*Cr*ws/ρ

    # State Equations

        # Sea Level

        du[1] = zdot

        # Barrier Island

        du[2] = 4*Qsf*((H+Dt)/(Dt*(2*H+Dt)))+2*zdot/α
        du[3] = (2*Qow/(2*H+Dt)-4*Qsf*((H+Dt)/(2*H+Dt)^2))
        du[4] = Qow_H/W-zdot
        du[5] = Qow_Bm/(H+zm-r/2)

        # Marsh-Lagoon

        du[6] = -E+Qow_Bl/(zL-zm)
        du[7] = E
        du[8] = (Fm+O)/β
        du[9] = -Fm-O+zdot
        du[10] = -2*E*(zL-zm)/bL+Fm*(bbm+bim)/bL+Fc+zdot

    nothing

end