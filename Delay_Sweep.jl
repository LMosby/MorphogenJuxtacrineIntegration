####################################
# Initialise simulation parameters #
####################################

# Check whether input argument(s) have been passed to code
inputs = ( if isempty(ARGS) ; ["1", "3", "10000"] ; else ARGS ; end )

# Read in parameter(s) from input arguments
const nSims = parse(Int, inputs[1])
const simType = parse(Int, inputs[2])
const ID = parse(UInt, inputs[3])

# Initialise the file names for saving data
const folder = @__DIR__
fileName = [folder * "/Data/$(Int(round(((ID - 1) * nSims))) + i)_DelaySweepData.jld2" for i in 1:nSims]

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



####################################
# Dimensional parameters           #
####################################

# Define Hh parameters
DHh::Float64 = 0.2844         # Diffusivity
DHh /= dx2                    # Non-dimensionalise diffusivity
kHh::Float64 = 0.008          # Degradation rate
vHh::Float64 = 0.1            # Production rate in the source region
vHhCell::Float64 = 0.178      # Production rate of the L2 cell

# Define Notch-Delta parameters
# In order: betaN0, betaN, alpha, threshN, m, gammaN, betaD, threshD, n, gammaD
betaN0::Float64 = 0.          # Basline Notch production rate
betaN::Float64 = 0.1          # Notch production rate
sigma::Float64 = 0.0509       # Dissociation constant in Notch production
m::Int = 2                    # Hill coefficient in Notch production
gammaN::Float64 = 0.008       # Notch degradation rate

betaD::Float64 = 0.1          # Delta production rate
epsilon::Float64 = 1.6487     # Dissociation constant in Delta production
n::Int = 3                    # Hill coefficient in Delta production
gammaD::Float64 = 0.008       # Delta degradation rate

# Initialise glial growth speed and position
gliaV::Float64 = 7.0833 * (10 ^ (-4)) # 1.7708 * (10 ^ (-4))

# Non-dimensionalise time by a characteristic glia time-scale
const gliaT::Float64 = cellWidth / gliaV
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

# Initialise Hh and Notch signalling thresholds
const threshHh_Low::Float64 = 1.2710
const threshHh_High::Float64 = 4.2559
threshNotchS_Low::Float64 = 1.3441
threshNotchS_High::Float64 = 4.7129
threshs::Vector{Float64} = [threshHh_Low, threshHh_High, threshNotchS_Low, threshNotchS_High]

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
const dtTarget::Float64 = 2.5e-6 # 1e-5
dtSim::Float64 = dtTarget
while dtSim > maxT
    println("Time-step too large for stability: $(dtSim) > $(1. / ((2. * DHh) + 1.)), $(dtSim) -> $(dtSim / 2.)")
    global dtSim /= 2.
end

# Initialise glial position
gliaXInit::Float64 = -4.

# Initialise MAPK activation profile
MAPK::Vector{Int} = [0 for i in 1:nCells]

# Initialise Dl production array
dlProd::Vector{Int} = [0 for i in 1:nCells]
notchProd::Vector{Int} = [1 for i in 1:nCells] # Can be everywhere 1 initially as is also gated by Dl production in neighbours

# Calculate times of cell differentiation from glia velocities
const waitTime::Float64 = cellWidth ./ gliaV

# Define delay time between Dl production being activated and then switched off
const dlProdDelayTime::Float64 = waitTime
currDlProdDelayTime::Vector{Float64} = 1e6 * ones(Float64, nCells)

# Initialise array to store differentiation status of cells
cellDiff::Vector{Int} = [0 for i in 1:nCells]
cellDiffHhProd::Vector{Int} = [0 for i in 1:nCells] # Stores which cells have 'locked in' their differentiation state wrt to producing Hh
cellHhProd::Vector{Int} = [0 for i in 1:nCells]     # Stores which cells actually produce Hh
cellDiffHh::Vector{Float64} = [-1. for i in 1:nCells]
cellDiffNotch::Vector{Float64} = [-1. for i in 1:nCells]

# Initialise storage of differentiate times for each cell
tauDiff = Vector{Int}(undef, 6)
tauMax::Float64 = 0.
tauDiffVaryFlag::Vector{Int} = [0 for i in 1:nCells]
tauDiffVary::Vector{Float64} = [-1. for i in 1:nCells]

# Calculate average rates of different processes in wild-type case
allGliaT::Vector{Float64} = ((cellGliaIntXs .- gliaXInit) ./ gliaV)
allGliaT[6] = allGliaT[1]
const avGliaT::Float64 = sum(allGliaT) / nCells # Also the average time for initiating Delta and Hh production
const avGliaK::Float64 = 1 / avGliaT
const avFateK::Float64 = 1 / (avGliaT + waitTime)
avKs::Vector{Float64} = [avGliaK, avGliaK, avFateK] # Hh production, Delta production, fate specification
# println(1 ./ avKs)

# Set wild-type conditions for cell fate decisions
hhT::Vector{Float64} = allGliaT
deltaT::Vector{Float64} = allGliaT
fateT::Vector{Float64} = allGliaT .+ waitTime
maxT::Float64 = maximum(fateT)



####################################
# Import libraries                 #
####################################

using DifferentialEquations, DiffEqCallbacks, LinearAlgebra, Symbolics, Statistics, Random, JLD2

include("Delay_Sweep_Functions.jl")



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
cbSet = CallbackSet(cbTermNew, cbGliaMAPK) # Note: positive concentrations are ensured by timestep stability instead of a callback function here



####################################
# Simulate Notch-Delta dynamics    #
####################################

# Loop over all simulations
@inbounds for i in 1:nSims
    println("\nSimulation $(i) / $(nSims)")

    # Reset key variables
    global vHhVec = [i < sIndex ? vHh : 0 for i in 1:nX]
    global MAPK = [0 for i in 1:nCells]
    global dlProd = [0 for i in 1:nCells]
    global currDlProdDelayTime = 1e6 * ones(Float64, nCells)
    global notchProd = [1 for i in 1:nCells]
    global cellDiff = [0 for j in 1:nCells]
    global cellDiffHhProd = [0 for j in 1:nCells]
    global cellHhProd = [0 for j in 1:nCells]
    global cellDiffHh = [-1. for j in 1:nCells]
    global cellDiffNotch = [-1. for j in 1:nCells]
    global tauDiffVaryFlag = [0 for i in 1:nCells]
    global tauDiffVary = [-1. for i in 1:nCells]

	# Seed random numbers for remaining input parameters
	nSeed::UInt = UInt(round(time())) + ID + (1e6i)
	Random.seed!(nSeed)

    # Extract event times using inverse transform sampling
    global ranTau = ranTimes(avKs[simType])
    # println(ranTau)

    # Define times of important events
    if(simType == 1)
        global hhT = deepcopy(ranTau)
    elseif(simType == 2)
        global deltaT = deepcopy(ranTau)
    elseif(simType == 3)
        global fateT = deepcopy(ranTau)
        global fateT[findall(x -> x > 10., fateT)] .= 10.
        global maxT = maximum(fateT)
    end

    # Save initial simulation parameters
    save_object(fileName[i], [nSeed, [hhT, deltaT, fateT], cellDiff, cellHhProd, cellDiffHh, cellDiffNotch])

    # Simulate system
    solFlag = solveSys(pHND, probSolver, cbSet)

    # Check if simulation completed successfully
    if(solFlag == 0)

        # Save simulation results
        println("$(cellDiff)")
        save_object(fileName[i], [nSeed, [hhT, deltaT, fateT], cellDiff, cellHhProd, cellDiffHh, cellDiffNotch])

    else

        # Terminate simulation
        println("Simulation error...")
        break

    end

end
