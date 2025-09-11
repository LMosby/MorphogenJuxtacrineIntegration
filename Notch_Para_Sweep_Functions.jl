#########################################
# List of functions for SweepNotchXX.jl #
#########################################

# Define the discretised PDE for Notch and Delta only
function solverNotch!(du, u, p, t)
	sigma, epsilon = p

    # Species (rows): (1 - nCells) Notch, (nCells+1 - 2*nCells) Delta

    du[1] = gammaN * ((notchProd[1] / (1 + ((sigma / u[nCells2]) ^ m))) - u[1])
    @inbounds for i in 2:nCells_1
        du[i] = gammaN * ((notchProd[i] / (1 + ((sigma / (u[nCells_1 + i] + u[nCells1 + i])) ^ m))) - u[i])
    end
    du[nCells] = gammaN * ((notchProd[nCells] / (1 + ((sigma / u[n2Cells_1]) ^ m))) - u[nCells])

    @inbounds for i in nCells1:n2Cells
        du[i] = gammaD * ((dlProd[i - nCells] * (1 / (1 + ((u[i - nCells] / epsilon) ^ n)))) - u[i])
    end

end

# Define function that solves the discretised PDE for Notch, Delta and Notch signalling
function solveSys(pNew::Vector{Float64}, prob, cbSet)

    # Reset glia position
    global gliaXInit = -4.
    global currGliaCell1 = 1

    # Reset differentiation and Delta and Notch production states
    global cellDiff = [0 for i in 1:nCells]
    global dlProd = [0 for i in 1:nCells]
    global notchProd = [1 for i in 1:nCells]
    global currDlProdDelayTime = 1e6 * ones(Float64, nCells)

    # Zero flags and arrays storing save times
    global finNotch = [-ones(Float64, 6) for i in 1:6]

    # println("Notch-Delta parameters: $(pNew)")
    sol = solve(remake(prob, p = pNew),
            Euler(),
            dt = dtSim,
            maxiters = 1e8,
            callback = cbSet, # Cell fate decisions occur in callback functions (see below)
            abstol = absTol,
            reltol = relTol,
            save_everystep = false,
            verbose = false)

    # Check if simulation completed successfully
    if(sol.retcode != :MaxIters)
        return 0
    else
        return -1
    end

end

# Define the modified callback function that defines when the Notch-Delta system has reached steady-state
function TerminateSteadyStateCondition(u, t, integrator)

    # Extract time derivatives of each concentration element
    intFromCache = first(get_tmp_cache(integrator))
    DiffEqBase.get_du!(intFromCache, integrator)

    # Check if the first 4 cells have decided their fates
    if((t < tauDiff[4]) || (any(cellDiff[1:4] .== 0)))
        return false
    else
        # println("System has reached steady-state after glia has interacted with all cells")
        return true
    end

end
TerminateSteadyStateAffect!(integrator) = terminate!(integrator)

# Define the callback condition and function to update glia position and Delta production
gliaCondition(u, t, integrator) = true
function gliaAffect!(integrator)

    # Check if glia has reached next cell (assume monotonically increasing glia position and that glia cannot 'skip' cells by moving too fast)
    if(currGliaCell1 <= nCells)
        if((integrator.t * gliaV) + gliaXInit >= cellGliaIntXs[currGliaCell1])

            # Set Delta production activation
            global dlProd[currGliaCell1] = 1
            global currDlProdDelayTime[currGliaCell1] = integrator.t + dlProdDelayTime

            # Check if activation occurred at cell 1 (cell 1 and 6 are activated simultaneously)
            if(currGliaCell1 == 1)
                global dlProd[nCells] = 1
                global currDlProdDelayTime[nCells] = integrator.t + dlProdDelayTime
            end

            # Increment activation flag
            global currGliaCell1 += 1

        end
    end
        
    # Check whether Dl production should be switched off in any cells
    if(any(dlProd .== 1))

        # Loop over all cells
        @inbounds for i in 1:nCells

            # Check if Dl (and Notch) production should be switched off in this cell
            if((dlProd[i] == 1) && (integrator.t > currDlProdDelayTime[i]))
                global dlProd[i] = 0
                global currDlProdDelayTime[i] = -1.
                global notchProd[i] = 0
            end

        end

    end

    # DEBUG: Print important parameters
    #= if((mod(round(integrator.t / dtSim), round(dtSaveTime / dtSim)) == 0))
        println("\nt = $(integrator.t), x = $((integrator.t * gliaV) + gliaXInit)")
        println(cellGliaIntXs)
        println(cellDiff)
        println(dlProd)
    end =#

    # Loop over all cells
    @inbounds for i in 1:nCells

        # Check for cell differentiation
        if((cellDiff[i] == 0) && (integrator.t >= tauDiff[i]))

            # Store the current Notch profile
            finNotch[i] = integrator.u[1:nCells]
            global cellDiff[i] = 1 # Cell fate is not important, only need to change value from 0

        end
        
    end

end

# Define function to calculate the robustness of the two sets of Notch concentrations at the times of cell fate decisions being made
function robustCalculation(allNotch::Vector{Float64}, threshs::Vector{Float64})

    # Calculate distance between Hh concentration and nearest threshold for each cell
    rob = Vector{Float64}(undef, 3)
    @inbounds for i in 1:3
        rob[i] = minimum([abs(allNotch[i] - threshs[j]) / allNotch[i] for j in 1:nThreshs])
    end

    return rob

end
