############################################
# List of functions for Individual_Sims.jl #
############################################

# Define struct for storing simulation data
struct SimData

	paras::Vector{Float64}                 # Store input parameters for Notch, Delta and Hh
    threshs::Vector{Float64}               # Store the thresholds used to make cell fate decisions
    x::Vector{Float64}                     # Store grid for Hh simulations
    cellMinInds::Vector{Int}               # Store cell positions

    allHh::Vector{Vector{Float64}}         # Store steady-state Hh profiles for all source widths
    allCellMeanHh::Vector{Vector{Float64}} # Store steady-state mean Hh concentrations at each cell for all source widths

    cellDiff::Vector{Float64}              # Store the cell states after differentiation
    cellDiffHhProd::Vector{Float64}        # Store which cells produce Hh after differentiation
    notch::Vector{Float64}                 # Store final Notch profile
    delta::Vector{Float64}                 # Store final Delta profile
    hh::Vector{Float64}                    # Store final Hh profile
    cellMeanHh::Vector{Float64}            # Store final mean Hh concentrations at each cell

end

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
    du[1] = betaN0 + (betaN * notchProd[1] * (1 / (1 + ((sigma / u[nCells2]) ^ m)))) - (gammaN * u[1])
    @inbounds for i in 2:nCells_1
        du[i] = betaN0 + (betaN * notchProd[i] * (1 / (1 + ((sigma / (u[nCells_1 + i] + u[nCells1 + i])) ^ m)))) - (gammaN * u[i])
    end
    du[nCells] = betaN0 + (betaN * notchProd[nCells] * (1 / (1 + ((sigma / u[n2Cells_1]) ^ m)))) - (gammaN * u[nCells])

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

    # Reset/zero system parameters if necessary
    if(nSims > 1)

        # Reset Hh production
        if(!hhOver)
            global vHhVec = [i < sIndex ? vHh : 0 for i in 1:nX]
        else
            global vHhVec = [vHh for i in 1:nX]
        end

        # Solve for the initial steady-state Hh profile (can use adaptive timestepping here as do not need to update cell state etc. during simulation)
        initSol = solve(probSolverHh, Rodas4P(), maxiters = 1e8, callback = cbSetHh, abstol = absTol, reltol = relTol, save_everystep = false, verbose = false)

        # Zero initial solver states
        global u0Hh = zeros(nX)
        global u0 = zeros(nTot)
        u0[n2Cells1:nTot] = initSol.u[end]

        # Reset glia position
        global gliaXInit = -4.
        global currGliaCell1 = 1

        # Reset MAPK activation state
        if(MAPKOver)
            global MAPK = [1 for i in 1:nCells] # [1 for i in 1:nCells]
        else
            global MAPK = [0 for i in 1:nCells] # [1 for i in 1:nCells]
        end

        # Reset Delta and Notch production states
        global dlProd = deepcopy(MAPK)
        global notchProd = [1 for i in 1:nCells]
        global currDlProdDelayTime = 1e6 * ones(Float64, nCells)

        # Zero arrays storing cell fate decisions
        global cellDiff = [0 for i in 1:nCells]
        global cellDiffHhProd = [0 for i in 1:nCells]
        global cellHhProd = [0 for i in 1:nCells]

    end

    # Zero flags and arrays storing save times
    global currSaveTime = 0.
    empty!(allT)
    empty!(allSols)
    global allDiffHh = [-ones(Float64, 6) for i in 1:6]
    global allDiffDelta = [-ones(Float64, 6) for i in 1:6]
    global allDiffNotch = [-ones(Float64, 6) for i in 1:6]

    # Store simulation results
    notch = Vector{Float64}(undef, nX)
    delta = Vector{Float64}(undef, nX)
    hh = Vector{Float64}(undef, nX)

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

    # Store simulation results if system was terminated before solver reached MaxIters
    if(sol.retcode != :MaxIters)
        notch = [sol.u[end][i] for i in 1:nCells]
        delta = [sol.u[end][i] for i in nCells1:n2Cells]
        hh = [sol.u[end][i] for i in n2Cells1:nTot]
    end
    
    return SimData([pNew..., gliaV], 
                    threshs, 
                    x, 
                    cellMinInds, 
                    allHh, 
                    [[mean(allHh[i][cellMinInds[j]:cellMaxInds[j]]) for j in 1:nCells] for i in eachindex(allHh)], 
                    cellDiff, 
                    cellHhProd, 
                    notch, 
                    delta, 
                    hh, 
                    [mean(hh[cellMinInds[i]:cellMaxInds[i]]) for i in 1:nCells])

end

# Define the modified callback function that defines when the Notch-Delta system has reached steady-state
function TerminateSteadyStateCondition(u, t, integrator)

    # Extract time derivatives of each concentration element
    intFromCache = first(get_tmp_cache(integrator))
    DiffEqBase.get_du!(intFromCache, integrator)

    # Check if all cells have differentiated
    if(t < 6.5) # any(cellDiff .== 0))
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
    if(notchUnder)
        currNotch = zeros(nCells)
    end
    @views currDelta = integrator.u[nCells1:n2Cells]

    # Calculate the mean (and standard deviation) of the Hh signal observed by each cell
    @views currCellMeanHh::Vector{Float64} = [mean(integrator.u[cellMinIndsN2Cells[i]:cellMaxIndsN2Cells[i]]) for i in 1:nCells]
    # @views cellStdsHh::Vector{Vector{Float64}} = [cellStds(Hh[i], cellMinInds, cellMaxInds, nCells) for i in 1:nRegions]

    # Check if glia has reached next cell (assume monotonically increasing glia position and that glia cannot 'skip' cells by moving too fast)
    if((!MAPKOver) && (currGliaCell1 <= nCells))
        if((integrator.t * gliaV) + gliaXInit >= cellGliaIntXs[currGliaCell1])

            # Set MAPK activation
            # Do not allow MAPK activation in cell 5 at any time
            # MAPK activation occurs in cell 6 at the same time as cell 1, not when glia position >= cellGliaIntXs[6]
            # if((currGliaCell1 < 5) || ((currGliaCell1 == 5) && (hhOver)))
                global MAPK[currGliaCell1] = 1
                global dlProd[currGliaCell1] = 1
                global currDlProdDelayTime[currGliaCell1] = integrator.t + dlProdDelayTime
            # end

            # Check if activation occurred at cell 1 (cell 1 and 6 are activated simultaneously)
            if(currGliaCell1 == 1)
                global MAPK[nCells] = 1
                global dlProd[nCells] = 1
                global currDlProdDelayTime[nCells] = integrator.t + dlProdDelayTime
            end

            # Increment MAPK activation flag
            global currGliaCell1 += 1

        end
    end

    # Check whether Dl production should be switched off in any cells
    if(stopDlProd)
        
        # Check whether Dl production is also switched off in cell 3
        if((!noStopDlProd3) && (any(dlProd .== 1)))

            # Loop over all cells
            @inbounds for i in 1:nCells

                # Check if Dl (and Notch) production should be switched off in this cell
                if((dlProd[i] == 1) && (integrator.t > currDlProdDelayTime[i]))
                    global dlProd[i] = 0
                    global currDlProdDelayTime[i] = -1.
                    if(stopNotchProd)
                        global notchProd[i] = 0
                    end
                end

            end

        elseif((noStopDlProd3) && (any(dlProd[[1, 2, 4, 5, 6]] .== 1)))

            # Loop over all cells
            @inbounds for i in [1, 2, 4, 5, 6]

                # Check if Dl production should be switched off in this cell
                if((dlProd[i] == 1) && (integrator.t > currDlProdDelayTime[i]))
                    global dlProd[i] = 0
                    global currDlProdDelayTime[i] = 0.
                    if(stopNotchProd)
                        global notchProd[i] = 0
                    end
                end

            end

        end

    end

    # DEBUG: Print important parameters
    #= if((mod(round(integrator.t / dtSim), round(dtSaveTime / dtSim)) == 0)) # || (any(MAPK .== 1)))
        println("\nt = $(integrator.t), x = $((integrator.t * gliaV) + gliaXInit)")
        println(cellGliaIntXs)
        println(cellDiffHhProd)
        println(cellHhProd)
        println(MAPK)
        println(dlProd)
        println(cellDiff)
        println(currCellMeanHh ./ (vHh / kHh))
        println(currNotch ./ (betaN / gammaN))
        println(currDelta ./ (betaD / gammaD))
    end =#

    # Check for cell differentiation
    @inbounds for i in 1:nCells

        # Check if MAPK is activated inside cell
        if(MAPK[i] == 1)

            # Check for Hh production in cell (can happen for any cell at any time)
            if((!hhProdKO) && (!hhOver) && (currCellMeanHh[i] >= threshHh_Low) && (currNotch[i] <= threshNotchS_Low) && (cellDiffHhProd[i] == 0))
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

                # DEBUG: Set Hh profile to steady-state profile immediately
                if(instantHh)
                    integrator.u[n2Cells1:nTot] = allHh[findall(cellHhProd .== 1)[end] + 1]
                    currCellMeanHh = [mean(integrator.u[cellMinIndsN2Cells[i]:cellMaxIndsN2Cells[i]]) for i in 1:nCells]
                end

            end

            # Check if cell fate decision has been made for cell already (turning on Hh production does not count)
            # Also check if specified "time of differentiation" has been reached
            if((cellDiff[i] == 0) && (integrator.t >= tauDiff[i]))

                # Check for L2 or L3 cell
                if(currCellMeanHh[i] >= threshHh_High)

                    # Check for L2 cell
                    if(currNotch[i] < threshNotchS_Low)
                        # println("\nL2 differentiation in cell $(i) at t = $(integrator.t) (Hh = $(currCellMeanHh[i]), N = $(currNotch[i]), D = $(currDelta[i]))\n")

                        # Save cell fate
                        global cellDiff[i] = 2
                        global cellDiffHhProd[i] = 1
                        global allDiffHh[i] = currCellMeanHh
                        global allDiffDelta[i] = currDelta
                        global allDiffNotch[i] = currNotch

                    # Check for L3 cell
                    elseif(currNotch[i] >= threshNotchS_Low)
                        # println("\nL3 differentiation in cell $(i) at t = $(integrator.t) (Hh = $(currCellMeanHh[i]), N = $(currNotch[i]), D = $(currDelta[i]))\n")

                        # Save cell fate
                        global cellDiff[i] = 3
                        global cellDiffHhProd[i] = 1
                        global allDiffHh[i] = currCellMeanHh
                        global allDiffDelta[i] = currDelta
                        global allDiffNotch[i] = currNotch

                    end

                # Check for L1 or L4 cell
                elseif((currCellMeanHh[i] >= threshHh_Low) && (currCellMeanHh[i] < threshHh_High))

                    # Check for L1 cell
                    if(currNotch[i] < threshNotchS_High) # (currNotch[i] >= threshNotchS_Low)
                        # println("\nL1 differentiation in cell $(i) at t = $(integrator.t) (Hh = $(currCellMeanHh[i]), N = $(currNotch[i]), D = $(currDelta[i]))\n")

                        # Save cell fate
                        global cellDiff[i] = 1
                        global cellDiffHhProd[i] = 1
                        global allDiffHh[i] = currCellMeanHh
                        global allDiffDelta[i] = currDelta
                        global allDiffNotch[i] = currNotch

                    # Check for L4 cell
                    elseif(currNotch[i] >= threshNotchS_High)
                        # println("\nL4 differentiation in cell $(i) at t = $(integrator.t) (Hh = $(currCellMeanHh[i]), N = $(currNotch[i]), D = $(currDelta[i]))\n")

                        # Save cell fate
                        global cellDiff[i] = 4
                        global cellDiffHhProd[i] = 1
                        global allDiffHh[i] = currCellMeanHh
                        global allDiffDelta[i] = currDelta
                        global allDiffNotch[i] = currNotch

                    end

                # Check for L5 cell
                elseif(currCellMeanHh[i] < threshHh_Low)
                    # println("\nL5 differentiation in cell $(i) at t = $(integrator.t) (Hh = $(currCellMeanHh[i]), N = $(currNotch[i]), D = $(currDelta[i]))\n")

                    # Save cell fate
                    global cellDiff[i] = 5
                    global cellDiffHhProd[i] = 1
                    global allDiffHh[i] = currCellMeanHh
                    global allDiffDelta[i] = currDelta
                    global allDiffNotch[i] = currNotch

                end

            end

        # Check for apoptotic cell fate (does not depend on time since MAPK activation)
        elseif((currCellMeanHh[i] < threshHh_Low) && (currNotch[i] >= threshNotchS_Low)) # threshNotchS_High

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
                    global allDiffHh[i] = currCellMeanHh
                    global allDiffDelta[i] = currDelta
                    global allDiffNotch[i] = currNotch

                end

            end

        elseif(tauDiffVary[i] > 0)

            # Reset delay times for apoptosis
            tauDiffVaryFlag[i] = 0
            tauDiffVary[i] = -1.

        end

    end

end

# Define callback condition function that checks for when the current simulation time passes regular thresholds
function saveCondition(u, t, integrator)

    # Check if time has passed current threshold
    if((t > tSaveMax) || (t < currSaveTime - dtEps))
		return false
	else
		return true
	end

end

# Define affect of successful callback condition that saves current simulation state if simulation time passes regular threshold
function saveAffect!(integrator)

	# Store current total morphogen and expander concentrations
	# println("$(currSaveTime), $(integrator.t)")
	push!(allT, integrator.t)
	# @views push!(allSols, [integrator.u[1:nCells], integrator.u[nCells1:n2Cells], integrator.u[n2Cells1:nTot]])
    @views push!(allSols, [integrator.u[1:nCells], integrator.u[nCells1:n2Cells], [mean(integrator.u[cellMinIndsN2Cells[i]:cellMaxIndsN2Cells[i]]) for i in 1:nCells], [std(integrator.u[cellMinIndsN2Cells[i]:cellMaxIndsN2Cells[i]]) for i in 1:nCells]])
	
	# Update time threshold
	global currSaveTime += dtSaveTime
	while currSaveTime < integrator.t
		global currSaveTime += dtSaveTime
	end

end
