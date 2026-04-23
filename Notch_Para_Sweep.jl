####################################
# Initialise simulation parameters #
####################################

# Check whether input argument(s) have been passed to code
inputs = ( if isempty(ARGS) ; ["1"] ; else ARGS ; end )

# Read in parameter(s) from input arguments
const ID = parse(UInt, inputs[1])

# Initialise the file name for saving data
const folder = @__DIR__
fileName = folder * "/Data/Notch_Sweep_Data.jld2"

# Initialise number of systems to sweep over
const nSims::Int = 2e6

# Initialise widths of source region and cells (7) total regions
const nCells::Int = 6
const nCells1::Int = nCells + 1
const nCells2::Int = nCells + 2
const nCells_1::Int = nCells - 1
const n2Cells::Int = nCells * 2
const n2Cells1::Int = n2Cells + 1
const n2Cells2::Int = n2Cells + 2
const n2Cells_1::Int = n2Cells - 1
const L::Float64 = 24. # Microns
const cellWidth::Float64 = L / nCells
cellGliaIntXs::Vector{Float64} = [(i - 0.5) * cellWidth for i in 1:nCells]

# Inititalise tolerances for steady-state detection
const absTol::Float64 = 1e-6
const relTol::Float64 = 1e-4

# Initialise ranges and minimum values for non-dimensional Notch-Delta parameters
# In order: sigma, epsilon, thresh 1, thresh 2 (threshs do not necessarily haev to be in order)
range::Vector{Float64} = [6., 6., 0.9, 0.9]
minVal::Vector{Float64} = [-3., -3., 0.1, 0.1]
prefactors::Vector{Float64} = [1., 1., 1., 1.]
gammaN::Float64 = 8e-3
gammaD::Float64 = 8e-3

# Set values for constant variables not being swept over
const m::Int = 2            # Hill coefficient in Notch production
const n::Int = 3            # Hill coefficient in Delta production

# Initialise number of thresholds
const nThreshs::Int = 2

# Initialise MAPK activation profile
MAPK::Vector{Int} = zeros(nCells)

# Initialise Dl production array
dlProd::Vector{Int} = deepcopy(MAPK)

# Initialise glial growth speed and position
gliaV::Float64 = 7.0833 * (10 ^ (-4)) # 1.7708 * (10 ^ (-4))

# Non-dimensionalise time by a characteristic glia time-scale
const gliaT::Float64 = cellWidth / gliaV
gammaN *= gliaT
gammaD *= gliaT
gliaV *= gliaT

# Initialise glial position
gliaXInit::Float64 = -4.
currGliaCell1::Int = 1

# Initialise temporal variables
dtSim::Float64 = 1e-2 # 1 / (4. * gammaN)

# Calculate times of cell differentiation from glia velocities
waitTime::Float64 = cellWidth / gliaV
tauDiffBuff::Vector{Float64} = ((cellGliaIntXs .- gliaXInit) ./ gliaV)
tauDiff::Vector{Float64} = [tauDiffBuff[i] .+ waitTime for i in 1:nCells]
tauDiff[nCells] = tauDiff[1]
# println(tauDiff)
tauDiffVaryFlag::Vector{Int} = [0 for i in 1:nCells]
tauDiffVary::Vector{Float64} = [-1. for i in 1:nCells]

# Initialise differentiation and Dl and Notch production arrays
cellDiff::Vector{Int} = [0 for i in 1:nCells]
dlProd::Vector{Int} = [0 for i in 1:nCells]
notchProd::Vector{Int} = [1 for i in 1:nCells] # Can be everywhere 1 initially as is also gated by Dl production in neighbours

# Define delay time between Dl production being activated and then switched off
dlProdDelayTime::Float64 = waitTime
currDlProdDelayTime::Vector{Float64} = 1e6 * ones(Float64, nCells)

# Initialise arrays to store Notch concentrations at each step of MAPK activation
finNotch::Vector{Vector{Float64}} = [-ones(Float64, 6) for i in 1:6]




####################################
# Import libraries                 #
####################################

using DifferentialEquations, DiffEqCallbacks, LinearAlgebra, Symbolics, Statistics, Random, JLD2

include("Notch_Para_Sweep_Functions.jl")

# Define struct for storing Hh Sweep data
mutable struct NotchData

    seed::Vector{UInt}                           # Store the seeds for each set of random number generation
	paras::Vector{Vector{Float64}}               # Store input parameters for each system
    finNotch::Vector{Vector{Vector{Float64}}}    # Store final Notch concentrations in each cell for each MAPK activation level
    expConstraints::Vector{Int}                  # Store whether or not the Hh profiles reproduce experimental findings
    cellRobs::Vector{Vector{Float64}}            # Store the robustnesses of each cell for each condition

end



####################################
# Initialise full system solvers   #
####################################

# Define initial conditions and time-scales for solver
u0 = zeros(n2Cells)
du0 = zeros(n2Cells)
tSpan = (0., Inf) # Always simulate until the system reaches steady-state

# Initialise Jacobian and problem for Hh and Notch-Delta equations
jacSolver = Symbolics.jacobian_sparsity((du, u) -> solverNotch!(du, u, [1., 1.], 0.), du0, u0)
solverSparse = ODEFunction(solverNotch! ; jac_prototype = float.(jacSolver))
probSolver = ODEProblem(solverSparse, u0, tSpan, [1., 1.], save_everystep = false)

# Initialise combined callback functions for detecting steady-state, updating the glia position and MAPK distribution, and saving data
cbTermNew = DiscreteCallback(TerminateSteadyStateCondition, TerminateSteadyStateAffect!)
cbGlia = DiscreteCallback(gliaCondition, gliaAffect!)
cbSet = CallbackSet(cbTermNew, cbGlia) # Note: positive concentrations are ensured by timestep stability instead of a callback function here



####################################
# Simulate Notch-Delta dynamics    #
####################################

# Initialise vector of output data
data = NotchData(Vector{UInt}(undef, nSims), 
                 Vector{Vector{Float64}}(undef, nSims), 
                 Vector{Vector{Vector{Float64}}}(undef, nSims), 
                 Vector{Int}(undef, nSims), 
                 Vector{Vector{Vector{Float64}}}(undef, nSims))

# Loop over parameters
@inbounds for i in 1:nSims
    if(mod(i, 1e3) == 0)
	    println("\nParameter set $(i) / $(nSims):")
    end

	# Seed random numbers for remaining input parameters
	nSeed::UInt = UInt(round(time())) + ID + (1e6i)
    # println(nSeed)
	
	# Use random number seed to generate random input parameters
	Random.seed!(nSeed)
	r = rand(4)
	r[1:2] = prefactors[1:2] .* (10 .^ ((range[1:2] .* r[1:2]) .+ minVal[1:2]))

    # Ensure that r[5] is the lower Notch threshold
    r[3:4] = minVal[3:4] .+ (range[3:4] .* r[3:4])
    if(r[4] < r[3])
        buff = copy(r[3])
        r[3] = r[4]
        r[4] = buff
    end

    # Simulate system
    solFlag = solveSys(r[1:2], probSolver, cbSet)

    # DEBUG: Print Notch output
    #= @inbounds for i in 1:nCells
        println(finNotch[i])
    end =#

    # Check if simulation completed successfully
    if(solFlag == 0)

        # Only interested in the Notch concentrations at the times of cell fate decisions being made
        cellFateConcs::Vector{Float64} = [finNotch[2][2], finNotch[3][3], finNotch[4][4], finNotch[1][2], finNotch[2][3], finNotch[3][4]]

        # Check whether profile satsifies experimental constraints
        expConstraints::Int = 0
        if((cellFateConcs[1] > r[4]) && (r[3] < cellFateConcs[2] < r[4]) && (cellFateConcs[3] > r[4])
               && (cellFateConcs[4] > r[3]) && (cellFateConcs[5] > r[3]) && (cellFateConcs[6] > r[3]))
            # println("Success 2")
            expConstraints = 2
        elseif((cellFateConcs[1] > r[4]) && (r[3] < cellFateConcs[2] < r[4]) && (cellFateConcs[3] > r[4]))
            # println("Success 1")
            expConstraints = 1
        end

        # Calculate robustness of Hh profiles
        cellRobs = robustCalculation(cellFateConcs, r[3:4])

        # Store data
        data.seed[i] = nSeed
        data.paras[i] = r
        data.finNotch[i] = finNotch
        data.expConstraints[i] = expConstraints
        data.cellRobs[i] = cellRobs

    else

        # Terminate simulation
        println("Simulation error...")
        break

    end

end



####################################
# Output data                      #
####################################

# Save data
save_object(fileName, data)
