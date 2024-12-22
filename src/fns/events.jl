function events()

    function condition1(out, u, t, integrator)

        out[1] = u[5]-u[3]  # Barrier Marsh Drowning
        out[2] = u[6]-u[5]  # Back Barrier Marsh Erosion
        out[3] = u[8]-u[7]  # Inland Marsh Erosion
        out[4] = Dmax-u[9]  # Marsh Drowning
        out[5] = u[7]-u[6]  # Marsh Filling 
    

    end

    function affect1!(integrator, idx)

        if idx == 1
            
            terminate!(integrator)
            println("Barrier Drowned at t = ", integrator.t)

        elseif idx == 2
            
            println("Barrier Marsh Eroded at t = ", integrator.t)

        elseif idx == 3

            println("Inland Marsh Eroded at t = ", integrator.t)

        elseif idx == 4

            terminate!(integrator)
            println("Marsh Drowned at t = ", integrator.t)

        elseif idx == 5

            terminate!(integrator)
            println("Marsh Filled at t = ", integrator.t)
            
        end

    end

    cb1 = VectorContinuousCallback(condition1,affect1!,5)

    function condition2(u, t, integrator)

        u[6] <= u[5]

    end

    function affect2!(integrator)

        integrator.u[6] = integrator.u[5]

    end

    cb2 = DiscreteCallback(condition2, affect2!)

    function condition3(u, t, integrator)

        u[7] >= u[8]

    end

    function affect3!(integrator)

        integrator.u[7] = integrator.u[8]

    end

    cb3 = DiscreteCallback(condition3, affect3!)

    # Creating Callback Set

    cbs = CallbackSet(cb1,cb2,cb3)

    return cbs

end