########################################
# List of functions for Delay_Sweep.jl #
########################################

# Define the discretised PDE for only Hh
function solverHh!(du, u, p, t)
    DHh, kHh = p

    du[1] = (2 * DHh * u[2]) - (((2 * DHh) + kHh) * u[1]) + vHhVec[1]
    @inbounds for i in 2:nX_1
        du[i] = (DHh * (u[i - 1] + u[i + 1])) - (((2 * DHh) + kHh) * u[i]) + vHhVec[i]
    end
    du[nX] = (2 * DHh * u[nX_1]) - (((2 * DHh) + kHh) * u[nX]) + vHhVec[nX]

end

# Define the discretised PDE for Notch and Delta with Notch repressing Delta production
function solver!(du, u, p, t)
	betaN0, betaN, sigma, m, gammaN, betaD, epsilon, n, gammaD, DHh, kHh = p

    # Species (rows): (1-N) Hh, (N+1-N+6) Notch, (N+7-N+12) Delta

    # Notch and Delta evolution is on 1:nCells domain
    du[1] = betaN0 + (betaN * (1 / (1 + ((sigma / u[nCells2]) ^ m)))) - (gammaN * u[1])
    @inbounds for i in 2:nCells_1
        du[i] = betaN0 + (betaN * (1 / (1 + ((sigma / (u[nCells_1 + i] + u[nCells1 + i])) ^ m)))) - (gammaN * u[i])
    end
    du[nCells] = betaN0 + (betaN * (1 / (1 + ((sigma / u[n2Cells_1]) ^ m)))) - (gammaN * u[nCells])

    @inbounds for i in nCells1:n2Cells
        du[i] = (betaD * dlProd[i - nCells] * (1 / (1 + ((u[i - nCells] / epsilon) ^ n)))) - (gammaD * u[i])
    end

    # Hh evolution is on 1:N domain
    du[n2Cells1] = (2 * DHh * u[n2Cells2]) - (((2. * DHh) + kHh) * u[n2Cells1]) + vHhVec[1]
    @inbounds for i in n2Cells2:nTot_1
        du[i] = (DHh * (u[i - 1] + u[i + 1])) - (((2. * DHh) + kHh) * u[i]) + vHhVec[i - n2Cells]
    end
    du[nTot] = (2 * DHh * u[nTot_1]) - (((2. * DHh) + kHh) * u[nTot]) + vHhVec[nX]

end

# Define function that solves the discretised PDE for Notch, Delta and Notch signalling
function solveSys(pNew::Vector{Float64}, prob, cbSet)

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

# Define function to calculate times when random events occur
function ranTimes(rate::Float64)

    # Sample six random values from a uniform distribution
    ranVals::Vector{Float64} = rand(6)

    # Use inverse transform sampling to extract event times
    return .-(1 ./ rate) .* log.(1. .- ranVals)

end

# Define the modified callback function that defines when the Notch-Delta system has reached steady-state
function TerminateSteadyStateCondition(u, t, integrator)

    # Extract time derivatives of each concentration element
    intFromCache = first(get_tmp_cache(integrator))
    DiffEqBase.get_du!(intFromCache, integrator)

    # Check if all cells have differentiated
    if((t < maxT) || (any(cellDiff .== 0)))
        return false
    else
        return true
    end

end
TerminateSteadyStateAffect!(integrator) = terminate!(integrator)

# Define the callback condition and function to update the glia position and MAPK distribution
gliaMAPKCondition(u, t, integrator) = true
function gliaMAPKAffect!(integrator)

    # Extract current Notch and Delta concentrations
    @views currNotch = integrator.u[1:nCells]
    @views currDelta = integrator.u[nCells1:n2Cells]

    # Calculate the mean (and standard deviation) of the Hh signal observed by each cell
    @views currCellMeanHh::Vector{Float64} = [mean(integrator.u[cellMinIndsN2Cells[i]:cellMaxIndsN2Cells[i]]) for i in 1:nCells]
    # @views cellStdsHh::Vector{Vector{Float64}} = [cellStds(Hh[i], cellMinInds, cellMaxInds, nCells) for i in 1:nRegions]

    # Loop over all cells
    @inbounds for i in 1:nCells

        # Check if Dl production should be switched on in this cell
        if((dlProd[i] == 0) && (currDlProdDelayTime[i] > 0.) && (cellDiff[i] != -1) && (integrator.t > deltaT[i]))
            global dlProd[i] = 1
        end

        # Check if Dl (and Notch) production should be switched off in this cell
        if((dlProd[i] == 1) && (integrator.t > currDlProdDelayTime[i]))
            global dlProd[i] = 0
            global currDlProdDelayTime[i] = -1.
            global notchProd[i] = 0
        end

    end

    # DEBUG: Print important parameters
    #= if((mod(round(integrator.t / dtSim), round(waitTime / (20. * dtSim))) == 0)) # || (any(MAPK .== 1)))
        println("\nt = $(integrator.t), x = $((integrator.t * gliaV) + gliaXInit)")
        println(cellGliaIntXs)
        println(cellDiffHhProd)
        println(cellHhProd)
        println(dlProd)
        println(cellDiff)
        println(currCellMeanHh ./ (vHh / kHh))
        println(currNotch ./ (betaN / gammaN))
        println(currDelta ./ (betaD / gammaD))
    end =#
    #= @inbounds for i in 1:nCells
        if(round(integrator.t / dtSim) == round(fateT[i] / dtSim))
            println("\nt = $(integrator.t), cell $(i) can now differentiate")
        end
    end =#

    # Check for cell differentiation
    @inbounds for i in 1:nCells

        # Check for Hh production in cell (can happen for any cell at any time)
        if((cellDiffHhProd[i] == 0) && (integrator.t >= hhT[i]) && (currCellMeanHh[i] >= threshHh_Low) && (currNotch[i] <= threshNotchS_Low))
            # println("\nInitiating Hh production in cell $(i) at t = $(integrator.t) (Hh = $(currCellMeanHh[i]), N = $(currNotch[i]))\n")

            # Save cell fate
            global cellDiffHhProd[i] = 1

            # Update Hh production array
            global cellHhProd[i] = 1
            @inbounds for j in 1:nX

                # Photo-receptor production
                if(j < sIndex)
                    vHhVec[j] = vHh
                else

                    # Lamina cell production
                    @inbounds for k in 1:nCells
                        if((j >= cellMinInds[k]) && (j <= cellMaxInds[k]) && (cellHhProd[k] == 1))
                            vHhVec[j] = vHhCell
                        end
                    end

                end

            end

            # DEBUG: Check for incorrect differentiation ordering
            #= if(any(cellDiff[i] .!= 0))
                println("Turned on Hh production after differentiation in cell $(i)")
            end =#

        end

        # Check if cell fate decision has been made for cell already (turning on Hh production does not count)
        # Also check if specified "time of differentiation" has been reached
        if((cellDiff[i] == 0) && (integrator.t >= fateT[i]))

            # Check for L2 or L3 cell
            if(currCellMeanHh[i] >= threshHh_High)

                # Check for L2 cell
                if(currNotch[i] < threshNotchS_Low)
                    # println("\nL2 differentiation in cell $(i) at t = $(integrator.t) (Hh = $(currCellMeanHh[i]), N = $(currNotch[i]), D = $(currDelta[i]))\n")

                    # Save cell fate
                    global cellDiff[i] = 2
                    global cellDiffHhProd[i] = 1
                    global cellDiffHh[i] = currCellMeanHh[i]
                    global currDlProdDelayTime[i] = integrator.t + dlProdDelayTime
                    global cellDiffNotch[i] = currNotch[i]

                # Check for L3 cell
                elseif(currNotch[i] >= threshNotchS_Low)
                    # println("\nL3 differentiation in cell $(i) at t = $(integrator.t) (Hh = $(currCellMeanHh[i]), N = $(currNotch[i]), D = $(currDelta[i]))\n")

                    # Save cell fate
                    global cellDiff[i] = 3
                    global cellDiffHhProd[i] = 1
                    global cellDiffHh[i] = currCellMeanHh[i]
                    global currDlProdDelayTime[i] = integrator.t + dlProdDelayTime
                    global cellDiffNotch[i] = currNotch[i]

                end

            # Check for L1 or L4 cell
            elseif((currCellMeanHh[i] >= threshHh_Low) && (currCellMeanHh[i] < threshHh_High))

                # Check for L1 cell
                if(currNotch[i] < threshNotchS_High) # (currNotch[i] >= threshNotchS_Low)
                    # println("\nL1 differentiation in cell $(i) at t = $(integrator.t) (Hh = $(currCellMeanHh[i]), N = $(currNotch[i]), D = $(currDelta[i]))\n")

                    # Save cell fate
                    global cellDiff[i] = 1
                    global cellDiffHhProd[i] = 1
                    global cellDiffHh[i] = currCellMeanHh[i]
                    global currDlProdDelayTime[i] = integrator.t + dlProdDelayTime
                    global cellDiffNotch[i] = currNotch[i]

                # Check for L4 cell
                elseif(currNotch[i] >= threshNotchS_High)
                    # println("\nL4 differentiation in cell $(i) at t = $(integrator.t) (Hh = $(currCellMeanHh[i]), N = $(currNotch[i]), D = $(currDelta[i]))\n")

                    # Save cell fate
                    global cellDiff[i] = 4
                    global cellDiffHhProd[i] = 1
                    global cellDiffHh[i] = currCellMeanHh[i]
                    global currDlProdDelayTime[i] = integrator.t + dlProdDelayTime
                    global cellDiffNotch[i] = currNotch[i]

                end

            # Check for L5 cell
            elseif((currCellMeanHh[i] < threshHh_Low) && (currNotch[i] < threshNotchS_Low))
                # println("\nL5 differentiation in cell $(i) at t = $(integrator.t) (Hh = $(currCellMeanHh[i]), N = $(currNotch[i]), D = $(currDelta[i]))\n")

                # Save cell fate
                global cellDiff[i] = 5
                global cellDiffHhProd[i] = 1
                global cellDiffHh[i] = currCellMeanHh[i]
                global currDlProdDelayTime[i] = integrator.t + dlProdDelayTime
                global cellDiffNotch[i] = currNotch[i]

            end

        end

        # Check for apoptotic cell fate (does not depend on MAPK activation state)
        if((currCellMeanHh[i] < threshHh_Low) && (currNotch[i] >= threshNotchS_Low)) # threshNotchS_High

            # Check if cell has already differentiated
            if(cellDiff[i] == 0)

                # Check if cell has waited long enough before undergoing apoptosis
                if((tauDiffVaryFlag[i] == 0) && (tauDiffVary[i] < 0))

                    # Set new delay time to reach before undergoing apoptosis
                    tauDiffVary[i] = integrator.t + waitTime

                elseif((tauDiffVaryFlag[i] == 0) && (tauDiffVary[i] > 0))

                    # Check if new delay time has yet been reached
                    if(integrator.t > tauDiffVary[i])
                        tauDiffVaryFlag[i] = 1
                    end

                end
                
                # Cell undergoes apoptosis
                if(tauDiffVaryFlag[i] == 1)
                    # println("\nApoptosis in cell $(i) at t = $(integrator.t) (Hh = $(currCellMeanHh[i]), N = $(currNotch[i]), D = $(currDelta[i]))\n")

                    # Save cell fate
                    global cellDiff[i] = -1
                    global cellDiffHhProd[i] = 1
                    global cellDiffHh[i] = currCellMeanHh[i]
                    global cellDiffNotch[i] = currNotch[i]
                    
                    # Also turn off Delta and Notch production
                    global dlProd[i] = 0
                    global currDlProdDelayTime[i] = -1.
                    global notchProd[i] = 0

                end

            end

        elseif(tauDiffVary[i] > 0)

            # Reset delay times for apoptosis
            tauDiffVaryFlag[i] = 0
            tauDiffVary[i] = -1.

        end

    end

end
