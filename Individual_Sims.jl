####################################
# Initialise simulation parameters #
####################################

# Check whether input argument(s) have been passed to code
inputs = ( if isempty(ARGS) ; ["1"] ; else ARGS ; end )

# Read in parameter(s) from input arguments
const ID = parse(UInt, inputs[1])

# Set number of simulation loops
const nSims::Int = 1

# Uncomment options below (and final line recalculating nSims) to sweep over different parameter values
# paraVaryArr::Vector{Float64} = [0., 0.25, 0.5, 0.75, 1., 1.25] # waitTime
# paraVaryArr::Vector{Float64} = [0.7, 0.8, 0.9, 0.95, 1., 1. / 0.95, 1. / 0.9, 1. / 0.8, 1. / 0.7] # waitTime
# paraVaryArr::Vector{Float64} = [0., 1., 2., 3., 4.] # dlProdDelayTime
# paraVaryArr::Vector{Float64} = [0.5, 1., 2.] # sigma
# paraVaryArr::Vector{Float64} = [0.5, 1., 2.] # epsilon
# paraVaryArr::Vector{Float64} = [10 ^ (-3), 10 ^ (-2.5), 10 ^ (-2), 10 ^ (-1.5), 10 ^ (-1), 10 ^ (-0.5), 1.] # gammaN, gammaD
# paraVaryArr::Vector{Float64} = [1., 10 ^ (1), 10 ^ (2), 10 ^ (3)] # gliaV
# paraVaryArr::Vector{Float64} = [0.8, 0.9, 0.95, 1., 1.05, 1.1, 1.2, 1.3] # gliaV
# paraVaryArr::Vector{Float64} = [0.8, 1.] # gliaV
# const nSims::Int = length(paraVaryArr)

# Initialise the file name for saving data
const folder = @__DIR__
fileName_Start = folder * "/Data/$(ID)_"

# Initialise spatial variables
const nX::Int = 281
const nX_1::Int = nX - 1
const L::Float64 = 24. # Microns
const LTot::Float64 = 28. # Microns
const dx::Float64 = LTot / nX_1
const dx2::Float64 = dx ^ 2
x::Vector{Float64} = [(i - 1) * dx for i in 1:nX] .- (LTot - L)

# Initialise widths of source region and cells (7) total regions
const nRegions::Int = 7
const nCells::Int = nRegions - 1
const nCells1::Int = nCells + 1
const nCells2::Int = nCells + 2
const nCells_1::Int = nCells - 1
const n2Cells::Int = nCells * 2
const n2Cells1::Int = n2Cells + 1
const n2Cells2::Int = n2Cells + 2
const n2Cells_1::Int = n2Cells - 1
const nTot::Int = n2Cells + nX
const nTot_1::Int = nTot - 1
const cellWidth::Float64 = LTot / nRegions
cellMinInds::Vector{Int} = [findall((x .> (i - 2) * cellWidth) .&& (x .<= (i - 1) * cellWidth))[1] for i in 2:nRegions]
cellMinIndsN2Cells::Vector{Int} = cellMinInds .+ n2Cells
cellMaxInds::Vector{Int} = [findall((x .> (i - 2) * cellWidth) .&& (x .<= (i - 1) * cellWidth))[end] for i in 2:nRegions]
cellMaxIndsN2Cells::Vector{Int} = cellMaxInds .+ n2Cells
cellGliaIntXs::Vector{Float64} = [(i - 1.5) * cellWidth for i in 2:nRegions]

# Initialise indexes for the end of the intrinsic and L2 cell source regions
const sIndex::Int = cellMinInds[1]
const sIndex_1::Int = sIndex - 1

# Inititalise tolerances for steady-state detection
const absTol::Float64 = 1e-6
const relTol::Float64 = 1e-4

# Set input flags
const instantHh::Bool = false
const fastHh::Bool = false
const hhProdKO::Bool = false
const hhOver::Bool = false
const MAPKOver::Bool = false
const notchUnder::Bool = false
const notchOverMild::Bool = false
const notchOverStrong::Bool = false
const notchOverMildDlUnder::Bool = false
const notchOverStrongDlUnder::Bool = false
const stopDlProd::Bool = true
const stopNotchProd::Bool = true
const noStopDlProd3::Bool = false



####################################
# Dimensional parameters           #
####################################

# Define Hh parameters
DHh::Float64 = 0.2844                             # Diffusivity
DHh /= dx2                                        # Non-dimensionalise diffusivity
kHh::Float64 = 0.008                              # Degradation rate
vHh::Float64 = 0.1                                # Production rate in the source region
vHhCell::Float64 = 0.178                          # Production rate of the L2 cell

# Define Notch-Delta parameters
betaN::Float64 = 0.1                              # Notch production rate
gammaN::Float64 = 0.008                           # Notch degradation rate
sigma::Float64 = 0.0509                           # Dissociation constant in Notch production
# sigma::Float64 = 0.1019                           # Generates on/off Notch concentrations
# sigma::Float64 = 28.3075                          # Generates correct phenotypes for slow degradation rate simulations
m::Int = 2                                        # Hill coefficient in Notch production

if((notchOverMildDlUnder) || (notchOverStrongDlUnder))
    betaD::Float64 = 0.                           # Delta production rate
else
    betaD::Float64 = 0.1
end
gammaD::Float64 = 0.008                           # Delta degradation rate
epsilon::Float64 = 1.6487                         # Dissociation constant in Delta production
# epsilon::Float64 = 0.0138                         # Generates on/off Delta concentrations
# epsilon::Float64 = 228.3762                       # Generates correct phenotypes for slow degradation rate simulations
n::Int = 3                                        # Hill coefficient in Delta production

# Initialise glial growth speed and position
gliaV0::Float64 =  7.0833 * (10 ^ (-4))
gliaV::Float64 = deepcopy(gliaV0)

# Initialise Hh and Notch signalling thresholds
const threshHh_Low::Float64 = 1.2710
const threshHh_High::Float64 = 4.2559
threshNotchS_Low::Float64 =  1.3441
# threshNotchS_Low::Float64 =  125.3039 # For slow degradation rate simulations
threshNotchS_High::Float64 = 4.7129
# threshNotchS_High::Float64 = 427.4269 # For slow degradation rate simulations

# Set parameters for slow simulations
# gammaN = 8e-5; sigma = 28.3075; gammaD = 8e-5; epsilon = 228.3762; threshNotchS_Low = 125.3039; threshNotchS_High = 427.4269;

# Combine threshold variables into vector
threshs::Vector{Float64} = [threshHh_Low, threshHh_High, threshNotchS_Low, threshNotchS_High]

# Define base-line Notch production (if any)
if((notchOverMild) || (notchOverMildDlUnder))
    betaN0::Float64 = ((threshNotchS_High + threshNotchS_Low) / 2.) * gammaN
    println("Base-line Notch production rate = $(betaN0)");
elseif((notchOverStrong) || (notchOverStrongDlUnder))
    betaN0::Float64 = (threshNotchS_High + (0.1 * (betaN / gammaN))) * gammaN
    println("Base-line Notch production rate = $(betaN0)");
else
    betaN0::Float64 = 0.
end

# Non-dimensionalise time by a characteristic glia time-scale
gliaT::Float64 = cellWidth / gliaV
DHh *= gliaT
kHh *= gliaT
vHh *= gliaT
vHhCell *= gliaT
betaN0 *= gliaT
betaN *= gliaT
gammaN *= gliaT
betaD *= gliaT
gammaD *= gliaT
gliaV *= gliaT

# Combine parameters into input vectors for solver
pHh::Vector{Float64} = [DHh, kHh]
pHND::Vector{Float64} = [betaN0, betaN, sigma, m, gammaN, betaD, epsilon, n, gammaD, DHh, kHh]

# Define maximum timestep for stability
maxT::Float64 = 1. / ((2. * DHh) + kHh)



####################################
# Other parameters                 #
####################################

# Initialise vector storing Hh production at each position
vHhVec::Vector{Float64} = [i < sIndex ? vHh : 0 for i in 1:nX]

# Initialise temporal variables
const dtTarget::Float64 = 2.5e-6
dtSim::Float64 = dtTarget
while dtSim > maxT
    println("Time-step too large for stability: $(dtSim) > $(1. / ((2. * DHh) + kHh)), $(dtSim) -> $(dtSim / 2.)")
    global dtSim /= 2.
end

# Initialise glial position
gliaXInit::Float64 = -4.
currGliaCell1::Int = 1

# Calculate times of cell differentiation from glia velocities
waitTime::Float64 = cellWidth ./ gliaV
if(MAPKOver)
    tauDiff::Vector{Float64} = [waitTime for i in 1:nCells]
else
    tauDiffBuff::Vector{Float64} = ((cellGliaIntXs .- gliaXInit) ./ gliaV)
    tauDiff::Vector{Float64} = [tauDiffBuff[i] .+ waitTime for i in 1:nCells]
    tauDiff[nCells] = tauDiff[1]
end
println(tauDiff)
tauDiffVaryFlag::Vector{Int} = [0 for i in 1:nCells]
tauDiffVary::Vector{Float64} = [-1. for i in 1:nCells]

# Initialise MAPK activation profile
if(MAPKOver)
    MAPK::Vector{Int} = [1 for i in 1:nCells] # [1 for i in 1:nCells]
else
    MAPK::Vector{Int} = [0 for i in 1:nCells] # [1 for i in 1:nCells]
end

# Initialise Dl and Notch production arrays
dlProd::Vector{Int} = deepcopy(MAPK)
notchProd::Vector{Int} = [1 for i in 1:nCells] # Can be everywhere 1 initially as is also gated by Dl production in neighbours

# Define delay time between Dl production being activated and then switched off
dlProdDelayTime::Float64 = waitTime
currDlProdDelayTime::Vector{Float64} = 1e6 * ones(Float64, nCells)

# Initialise array to store differentiation status of cells
cellDiff::Vector{Int} = [0 for i in 1:nCells]
cellDiffHhProd::Vector{Int} = [0 for i in 1:nCells] # Stores which cells have 'locked in' their differentiation state wrt to producing Hh
cellHhProd::Vector{Int} = [0 for i in 1:nCells]     # Stores which cells actually produce Hh

# Initialise timescales and arrays for saving data during simulations
currSaveTime::Float64 = 0.
dtSaveTime::Float64 = waitTime / 20.
tSaveMax::Float64 = ((2 * cellGliaIntXs[end]) - cellGliaIntXs[end - 1]) ./ gliaV
const dtEps::Float64 = 1e-6
allT::Vector{Float64} = []
allSols::Vector{Vector{Vector{Float64}}} = []
allDiffHh::Vector{Vector{Float64}} = []
allDiffDelta::Vector{Vector{Float64}} = []
allDiffNotch::Vector{Vector{Float64}} = []



####################################
# Import libraries                 #
####################################

using DifferentialEquations, DiffEqCallbacks, LinearAlgebra, Symbolics, Statistics, JLD2

include("Individual_Sims_Functions.jl")



####################################
# Initialise the PR Hh profile     #
####################################

# Define initial conditions and time-scales for solver
u0Hh = zeros(nX)
du0Hh = zeros(nX)
tSpan = (0., Inf) # Always simulate until the system reaches steady-state

# Initialise Jacobian and problem for only Hh equation
jacSolverHh = Symbolics.jacobian_sparsity((du, u) -> solverHh!(du, u, pHh, 0.), du0Hh, u0Hh)
solverHhSparse = ODEFunction(solverHh! ; jac_prototype = float.(jacSolverHh))
probSolverHh = ODEProblem(solverHhSparse, u0Hh, tSpan, pHh, save_everystep = false)

# Initialise combined callback functions for detecting steady state and ensuring positive concentrations
cbTermHh = TerminateSteadyState(absTol, relTol)
cbPosDomain = PositiveDomain()
cbSetHh = CallbackSet(cbTermHh, cbPosDomain)

# Solve for the initial steady-state Hh profile (can use adaptive timestepping here as do not need to update cell state etc. during simulation)
initSol = solve(probSolverHh, Rodas4P(), maxiters = 1e8, callback = cbSetHh, abstol = absTol, reltol = relTol, save_everystep = false, verbose = false)

# Store all possible steady-state Hh profiles
allHh = Vector{Vector{Float64}}(undef, nRegions)
@inbounds for i in 1:nRegions

    # Update indexes corresponding to cell boundary
    if(i == 1)
        global vHhVec = [j < sIndex ? vHh : 0 for j in 1:nX]
    else
        global vHhVec = [j < sIndex ? vHh :
                            (j <= cellMaxInds[i - 1] ? vHhCell : 0) for j in 1:nX]
    end

    # Solve for current Hh profile
    currSol = solve(probSolverHh, Rodas4P(), maxiters = 1e8, callback = cbSetHh, abstol = absTol, reltol = relTol, save_everystep = false, verbose = false)
    allHh[i] = currSol.u[end]

    # Reset production to baseline values
    if(i == nRegions)
        if(!hhOver)
            global vHhVec = [i < sIndex ? vHh : 0 for i in 1:nX]
        else
            global vHhVec = [vHh for i in 1:nX]
        end
    end

end



####################################
# Initialise full system solvers   #
####################################

# Define initial conditions and time-scales for solver
u0 = zeros(nTot)
u0[n2Cells1:nTot] = initSol.u[end] # Expect only small deviations from the analytical profile here due to discretisation of space
du0 = zeros(nTot)
tSpan = (0., Inf) # Always simulate until the system reaches steady-state

# Initialise Jacobian and problem for Hh and Notch-Delta equations
jacSolver = Symbolics.jacobian_sparsity((du, u) -> solver!(du, u, pHND, 0.), du0, u0)
solverSparse = ODEFunction(solver! ; jac_prototype = float.(jacSolver))
probSolver = ODEProblem(solverSparse, u0, tSpan, pHND, save_everystep = false)

# Initialise combined callback functions for detecting steady-state, updating the glia position and MAPK distribution, and saving data
cbTermNew = DiscreteCallback(TerminateSteadyStateCondition, TerminateSteadyStateAffect!)
cbGliaMAPK = DiscreteCallback(gliaMAPKCondition, gliaMAPKAffect!)
cbSave = DiscreteCallback(saveCondition, saveAffect!)
cbSet = CallbackSet(cbTermNew, cbGliaMAPK, cbSave) # Note: positive concentrations are ensured by timestep stability instead of a callback function here



####################################
# Simulate Notch-Delta dynamics    #
####################################

# Loop over solves
@inbounds for i in 1:nSims
    if(nSims > 1)
        println("\nSimulation $(i) / $(nSims)")
    end

    # Change values of parameters using loop ID
    if(nSims > 1)

        # Uncomment options below corresponding to options selected (uncommented) at the top of the code to sweep over different parameter values

        # Vary times of cell differentiation from glia velocities
        #= global waitTime = paraVaryArr[i] * (cellWidth ./ gliaV)
        if(MAPKOver)
            global tauDiff = [waitTime for i in 1:nCells]
        else
            global tauDiffBuff = ((cellGliaIntXs .- gliaXInit) ./ gliaV)
            global tauDiff = [tauDiffBuff[i] .+ waitTime for i in 1:nCells]
            tauDiff[nCells] = tauDiff[1]
        end
        println("Cell fate delay times = $(tauDiff)")
        
        # Can simultaneously vary delay time between Dl production being activated and then switched off
        global dlProdDelayTime = (3. - paraVaryArr[i]) * (cellWidth ./ gliaV)
        println("Delta production delay time = $(dlProdDelayTime)") =#

        # Vary delay time between Dl production being activated and then switched off
        #= global dlProdDelayTime = paraVaryArr[i] * waitTime
        println("Delta production delay time = $(dlProdDelayTime)") =#

        # Vary sigma
        #= pHND[3] = paraVaryArr[i] * sigma
        println("$(paraVaryArr[i]) x sigma = $(pHND[3])") =#

        # Vary epsilon
        #= pHND[7] = paraVaryArr[i] * epsilon
        println("$(paraVaryArr[i]) x epsilon = $(pHND[7])") =#

        # Vary Delta and Notch relaxation times
        #= global betaN = 0.1
        global gammaN = paraVaryArr[i] * 0.008
        global betaD = 0.1
        global gammaD = paraVaryArr[i] * 0.008
        println("$(paraVaryArr[i]) x gammaN,D = $(gammaN)")
        global sigma = 0.0041 .* (betaD / gammaD) # Need to recalculate parameters that depend on these rates
        global epsilon = 0.1319 .* (betaN / gammaN)
        global betaN *= gliaT # Need to re-non-dimensionalise the rates
        global gammaN *= gliaT
        global betaD *= gliaT
        global gammaD *= gliaT
        global pHND = [betaN0, betaN, sigma, m, gammaN, betaD, epsilon, n, gammaD, DHh, kHh]
        println(pHND)
        global gammaD *= gliaT
        global threshNotchS_Low = 0.1075 .* (betaN / gammaN)
        global threshNotchS_High = 0.3770 .* (betaN / gammaN)
        global threshs = [threshHh_Low, threshHh_High, threshNotchS_Low, threshNotchS_High]
        println(threshs) =#

        # Vary glia velocity with the times of cell differentiation and Delta production inactivation
        #= gliaV = paraVaryArr[i] * gliaV0
        println("$(paraVaryArr[i]) x gliaV = $(paraVaryArr[i] * gliaV0)")
        global DHh /= gliaT # Need to re-dimensionalise variables before non-dimensionalising them using new time-scale
        global kHh /= gliaT
        global vHh /= gliaT
        global vHhCell /= gliaT
        global betaN0 /= gliaT
        global betaN /= gliaT
        global gammaN /= gliaT
        global betaD /= gliaT
        global gammaD /= gliaT
        global gliaT = cellWidth / gliaV # Re-dimensionalise variables using new time-scale

        global DHh *= gliaT
        global kHh *= gliaT
        global vHh *= gliaT
        global vHhCell *= gliaT
        global betaN0 *= gliaT
        global betaN *= gliaT
        global gammaN *= gliaT
        global betaD *= gliaT
        global gammaD *= gliaT
        global gliaV *= gliaT

        global pHh = [DHh, kHh]
        global pHND = [betaN0, betaN, sigma, m, gammaN, betaD, epsilon, n, gammaD, DHh, kHh]
        global maxT = 1. / ((2. * DHh) + kHh) # Must re-calculate all other variables that depend on this time-scale
        global vHhVec = [i < sIndex ? vHh : 0 for i in 1:nX]
        global dtSim = dtTarget
        while dtSim > maxT
            println("Time-step too large for stability: $(dtSim) > $(1. / ((2. * DHh) + 1.)), $(dtSim) -> $(dtSim / 2.)")
            global dtSim /= 2.
        end
        global waitTime = cellWidth ./ gliaV
        if(MAPKOver)
            global tauDiff = [waitTime for i in 1:nCells]
        else
            global tauDiffBuff = ((cellGliaIntXs .- gliaXInit) ./ gliaV)
            global tauDiff = [tauDiffBuff[i] .+ waitTime for i in 1:nCells]
            global tauDiff[nCells] = tauDiff[1]
        end
        # println(tauDiff)
        global dlProdDelayTime = waitTime
        global dtSaveTime = waitTime / 20.
        global tSaveMax = ((2 * cellGliaIntXs[end]) - cellGliaIntXs[end - 1]) ./ gliaV

        global allHh = Vector{Vector{Float64}}(undef, nRegions)
        @inbounds for i in 1:nRegions

            # Update indexes corresponding to cell boundary
            if(i == 1)
                global vHhVec = [j < sIndex ? vHh : 0 for j in 1:nX]
            else
                global vHhVec = [j < sIndex ? vHh :
                                    (j <= cellMaxInds[i - 1] ? vHhCell : 0) for j in 1:nX]
            end

            # Solve for current Hh profile
            currSol = solve(probSolverHh, Rodas4P(), maxiters = 1e8, callback = cbSetHh, abstol = absTol, reltol = relTol, save_everystep = false, verbose = false)
            allHh[i] = currSol.u[end]

            # Reset production to baseline values
            if(i == nRegions)
                if(!hhOver)
                    global vHhVec = [i < sIndex ? vHh : 0 for i in 1:nX]
                else
                    global vHhVec = [vHh for i in 1:nX]
                end
            end

        end =#

        # Vary only the glia velocity (without the times of cell differentiation and Delta production inactivation)
        #= gliaV = paraVaryArr[i] * gliaV0
        println("$(paraVaryArr[i]) x gliaV = $(paraVaryArr[i] * gliaV0)")
        global DHh /= gliaT # Need to re-dimensionalise variables before non-dimensionalising them using new time-scale
        global kHh /= gliaT
        global vHh /= gliaT
        global vHhCell /= gliaT
        global betaN0 /= gliaT
        global betaN /= gliaT
        global gammaN /= gliaT
        global betaD /= gliaT
        global gammaD /= gliaT
        global waitTime *= gliaT # These are timescales so re-dimensionalise in reverse to other variables
        global dlProdDelayTime *= gliaT

        global gliaT = cellWidth / gliaV # Re-dimensionalise variables using new time-scale
        global DHh *= gliaT
        global kHh *= gliaT
        global vHh *= gliaT
        global vHhCell *= gliaT
        global betaN0 *= gliaT
        global betaN *= gliaT
        global gammaN *= gliaT
        global betaD *= gliaT
        global gammaD *= gliaT
        global waitTime /= gliaT # These are timescales so re-dimensionalise in reverse to other variables
        global dlProdDelayTime /= gliaT
        global gliaV *= gliaT

        global pHh = [DHh, kHh]
        global pHND = [betaN0, betaN, sigma, m, gammaN, betaD, epsilon, n, gammaD, DHh, kHh]
        global maxT = 1. / ((2. * DHh) + kHh) # Must re-calculate all other variables that depend on this time-scale
        global vHhVec = [i < sIndex ? vHh : 0 for i in 1:nX]
        global dtSim = dtTarget
        while dtSim > maxT
            println("Time-step too large for stability: $(dtSim) > $(1. / ((2. * DHh) + 1.)), $(dtSim) -> $(dtSim / 2.)")
            global dtSim /= 2.
        end
        if(MAPKOver)
            global tauDiff = [waitTime for i in 1:nCells]
        else
            global tauDiffBuff = ((cellGliaIntXs .- gliaXInit) ./ gliaV)
            global tauDiff = [tauDiffBuff[i] .+ waitTime for i in 1:nCells]
            global tauDiff[nCells] = tauDiff[1]
        end
        println(tauDiff)
        global dtSaveTime = (cellWidth ./ gliaV) / 20. # Still save same number of points between glia reaching new cells
        global tSaveMax = ((2 * cellGliaIntXs[end]) - cellGliaIntXs[end - 1]) ./ gliaV

        global allHh = Vector{Vector{Float64}}(undef, nRegions)
        @inbounds for i in 1:nRegions

            # Update indexes corresponding to cell boundary
            if(i == 1)
                global vHhVec = [j < sIndex ? vHh : 0 for j in 1:nX]
            else
                global vHhVec = [j < sIndex ? vHh :
                                    (j <= cellMaxInds[i - 1] ? vHhCell : 0) for j in 1:nX]
            end

            # Solve for current Hh profile
            currSol = solve(probSolverHh, Rodas4P(), maxiters = 1e8, callback = cbSetHh, abstol = absTol, reltol = relTol, save_everystep = false, verbose = false)
            allHh[i] = currSol.u[end]

            # Reset production to baseline values
            if(i == nRegions)
                if(!hhOver)
                    global vHhVec = [i < sIndex ? vHh : 0 for i in 1:nX]
                else
                    global vHhVec = [vHh for i in 1:nX]
                end
            end

        end =#

    end

    # Simulate system
    data::SimData = solveSys(pHND, probSolver, cbSet)
    println("$(reduce(vcat, [threshs[1:2] ./ (vHh / kHh), threshs[3:4] ./ (betaN / gammaN)]))")
    println("$(data.cellDiff)")
    println("$([allDiffHh[i][i] / (vHh / kHh) for i in 1:6])")
    println("$([allDiffNotch[i][i] / (betaN / gammaN) for i in 1:6])")
    println("$([allDiffDelta[i][i] / (betaD / gammaD) for i in 1:6])\n\n")

    # DEBUG: Print results of individual simulations
    #= @inbounds for i in 1:nCells
        println(allDiffNotch[i] ./ (betaN / gammaN))
    end
    println(cellDiffHhProd)
    println(cellHhProd) =#

    # Save data
    if(nSims == 1)
        save_object(fileName_Start * "Struct.jld2", data)
        save_object(fileName_Start * "TData.jld2", [allT, allSols, allDiffDelta, allDiffNotch, allDiffHh])
    else
        save_object(fileName_Start * "$(i)_Struct.jld2", data)
        save_object(fileName_Start * "$(i)_TData.jld2", [allT, allSols, allDiffDelta, allDiffNotch, allDiffHh])
    end

end
