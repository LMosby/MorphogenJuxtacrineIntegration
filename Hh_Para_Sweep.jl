####################################
# Import libraries                 #
####################################

using DifferentialEquations, DiffEqCallbacks, LinearAlgebra, Symbolics, Statistics, Random, JLD2



####################################
# Initialise simulation parameters #
####################################

# Check whether input argument(s) have been passed to code
inputs = ( if isempty(ARGS) ; ["1"] ; else ARGS ; end )

# Read in parameter(s) from input arguments
const ID = parse(UInt, inputs[1])

# Initialise the file name for saving data
const folder = @__DIR__
fileName = folder * "/Data/Hh_Sweep_Data.jld2"

# Initialise number of systems to sweep over
const nSims::Int = 2e6

# Initialise spatial variables
const L::Float64 = 28. # Microns
const nRegions::Int = 7
const nCells::Int = nRegions - 1
const cellWidth::Float64 = L / nRegions
const cellWidth2::Float64 = 2 * cellWidth

# Initialise ranges and minimum values for diffusivity, degradation and production parameters
range::Vector{Float64} = [4., 3., 0.9, 0.9]
minVal::Vector{Float64} = [0., -1., 0.1, 0.1]
prefactors::Vector{Float64} = [2.5, 1.]

# Initialise number of thresholds
const nThreshs::Int = 2



####################################
# Initialise functions             #
####################################

# Define struct for storing Hh Sweep data
mutable struct HhData

    seed::Vector{UInt}                           # Store the seeds for each set of random number generation
	paras::Vector{Vector{Float64}}               # Store input parameters for each system
    cellMeanHh::Vector{Vector{Vector{Float64}}}  # Store final mean Hh concentrations at each cell for each condition
    cellStdHh::Vector{Vector{Vector{Float64}}}   # Store final errors in Hh concentrations (due to averaging) at each cell for each condition
    expConstraints::Vector{Vector{Int}}          # Store whether or not the Hh profiles reproduce experimental findings
    cellRobs::Vector{Vector{Vector{Float64}}}    # Store the robustnesses of each cell for each condition

end

# Define function to calculate the Hh profiles for production by only the photo-receptor cells and after Hh production is initiated in cell 1
function HhCalculation(lam::Float64, vHhCell::Float64)

    # Calculate Hh profile for production by only photo-receptor cells
    # meanHh = [(sinh(cellWidth / lam) / (2 * sinh(L / lam))) * lam * ((exp(L / lam) * (exp(-(i * cellWidth) / lam) - exp(-((i + 1) * cellWidth) / lam))) + (exp(-L / lam) * (exp(((i + 1) * cellWidth) / lam) - exp((i * cellWidth) / lam)))) for i in 1:nCells]
    meanHh::Vector{Float64} = [(sinh(cellWidth / lam) / sinh(L / lam)) * lam * (sinh((L - (i * cellWidth)) / lam) - sinh((L - ((i + 1) * cellWidth)) / lam)) for i in 1:nCells]
    meanHh ./= cellWidth
    allMeanHh = copy(meanHh)

    # Calculate corresponding variances
    stdHh::Vector{Float64} = [((sinh(cellWidth / lam) / sinh(L / lam)) ^ 2) * ((cellWidth / 2.) - ((lam / 4.) * (sinh((2. * (L - ((i + 1) * cellWidth))) / lam) - sinh((2. * (L - (i * cellWidth))) / lam)))) for i in 1:nCells]
    stdHh = sqrt.((stdHh ./ cellWidth) .- (meanHh .^ 2))
    allStdHh = copy(stdHh)

    # Calculate Hh profile for production by photo-receptor cells and first lamina cell
    c::Float64 = ((vHhCell * tanh((L - cellWidth2) / lam)) + ((1. - vHhCell) * exp(cellWidth / lam) * tanh(cellWidth / lam) * ((tanh((L - cellWidth2) / lam) + 1.) / (tanh(cellWidth / lam) - 1.)))) / ((exp(-cellWidth2 / lam) * (1. - tanh((L - cellWidth2) / lam))) + ((1 + tanh(cellWidth / lam)) * ((tanh((L - cellWidth2) / lam) + 1.) / (tanh(cellWidth / lam) - 1.))))
    # a::Float64 = ((2. * c * exp(-cellWidth / lam)) - (1. - vHhCell)) / (2. * (cosh(cellWidth / lam) - sinh(cellWidth / lam)))
    d::Float64 = (((1. - vHhCell) * exp(-cellWidth / lam) * tanh(cellWidth / lam)) - (c * exp(-cellWidth2 / lam) * (tanh(cellWidth / lam) + 1.))) / (tanh(cellWidth / lam) - 1.)
    f::Float64 = ((c * exp(-(L + cellWidth2) / lam)) + (d * exp((cellWidth2 - L) / lam)) + (vHhCell * exp(-L / lam))) / (2. * cosh((L - cellWidth2) / lam))
    e::Float64 = f * exp((2. * L) / lam)
    meanHh[1] = (vHhCell * cellWidth) + (c * lam * (exp(-cellWidth / lam) - exp(-cellWidth2 / lam))) + (d * lam * (exp(cellWidth2 / lam) - exp(cellWidth / lam)))
    meanHh[2:nCells] = [(e * lam * (exp(-(i * cellWidth) / lam) - exp(-((i + 1) * cellWidth) / lam))) + (f * lam * (exp(((i + 1) * cellWidth) / lam) - exp((i * cellWidth) / lam))) for i in 2:nCells]
    meanHh ./= cellWidth
    allMeanHh = [allMeanHh, meanHh]

    # Calculate corresponding variances
    stdHh[1] = (((2. * c * d) + (vHhCell ^ 2)) * cellWidth) + ((((c ^ 2) * lam) / 2.) * (exp(-(2. * cellWidth) / lam) - exp(-(2. * cellWidth2) / lam))) + ((2. * c * vHhCell * lam) * (exp(-cellWidth / lam) - exp(-cellWidth2 / lam))) + ((2. * d * vHhCell * lam) * (exp(cellWidth2 / lam) - exp(cellWidth / lam))) + ((((d ^ 2) * lam) / 2.) * (exp((2. * cellWidth2) / lam) - exp((2. * cellWidth) / lam)))
    stdHh[2:nCells] = [((2. * e * f) * cellWidth) + ((((e ^ 2) * lam) / 2.) * (exp(-(2. * (i * cellWidth)) / lam) - exp(-(2. * ((i + 1) * cellWidth)) / lam))) + ((((f ^ 2) * lam) / 2.) * (exp((2. * ((i + 1) * cellWidth)) / lam) - exp((2. * (i * cellWidth)) / lam))) for i in 2:nCells]
    stdHh = sqrt.((stdHh ./ cellWidth) .- (meanHh .^ 2))
    allStdHh = [allStdHh, stdHh]

    return [allMeanHh, allStdHh]

end

# Define function to calculate the robustness of the two Hh profiles calculated above
function robustCalculation(allHh::Vector{Vector{Float64}}, threshs::Vector{Float64})

    # Calculate distance between Hh concentration and nearest threshold for each cell
    rob = Vector{Vector{Float64}}(undef, 2)
    @inbounds for i in 1:2
        rob[i] = Vector{Float64}(undef, nCells)

        @inbounds for j in 1:nCells

            rob[i][j] = minimum([abs(allHh[i][j] - threshs[k]) / allHh[i][j] for k in 1:nThreshs])

        end

    end

    return rob

end



####################################
# Calculate Hh profiles            #
####################################

# Initialise vector of output data
data = HhData(Vector{UInt}(undef, nSims), 
              Vector{Vector{Float64}}(undef, nSims), 
              Vector{Vector{Vector{Float64}}}(undef, nSims), 
              Vector{Vector{Vector{Float64}}}(undef, nSims), 
              Vector{Vector{Int}}(undef, nSims), 
              Vector{Vector{Vector{Float64}}}(undef, nSims))

# Loop over parameters
@inbounds for i in 1:nSims
    if(mod(i, 1e4) == 0)
	    println("\nParameter set $(i) / $(nSims):")
    end

	# Seed random numbers for remaining input parameters
	nSeed::UInt = UInt(round(time())) + ID + (1e6i)
	
	# Use random number seed to generate random input parameters
	Random.seed!(nSeed)
	r = rand(4)
	r[1:2] = prefactors .* (10 .^ ((range[1:2] .* r[1:2]) .+ minVal[1:2]))

    # Ensure that r[3] is the lower Hh threshold
    r[3:4] = minVal[3:4] .+ (range[3:4] .* r[3:4])
    if(r[4] < r[3])
        buff = copy(r[3])
        r[3] = r[4]
        r[4] = buff
    end

    # Calculate mean Hh across cells
    cellMeanHh, cellStdHh = HhCalculation(sqrt(r[1]), r[2])

    # Check whether profile satsifies experimental constraints
    expConstraints::Vector{Int} = zeros(2)
    if((r[3] < cellMeanHh[1][1] < r[4]) && (r[3] < cellMeanHh[1][2] < r[4]) &&
        (cellMeanHh[1][3] < r[3]) && (cellMeanHh[1][4] < r[3]) &&
        (cellMeanHh[1][5] < r[3]) && (cellMeanHh[1][6] < r[3]))
        # println("Success 1")
        expConstraints[1] = 1
    end
    if((cellMeanHh[2][1] > r[4]) && (cellMeanHh[2][2] > r[4]) && 
        (r[3] < cellMeanHh[2][3] < r[4]) && (r[3] < cellMeanHh[2][4] < r[4]) &&
        (cellMeanHh[2][5] < r[3]) && (cellMeanHh[2][6] < r[3]))
        # println("Success 2")
        expConstraints[2] = 1
    end

    # Calculate robustness of Hh profiles
    cellRobs = robustCalculation(cellMeanHh, r[3:4])

    # Store data
    data.seed[i] = nSeed
    data.paras[i] = r
    data.cellMeanHh[i] = cellMeanHh
    data.cellStdHh[i] = cellStdHh
    data.expConstraints[i] = expConstraints
    data.cellRobs[i] = cellRobs

end



####################################
# Output data                      #
####################################

# Save data
save_object(fileName, data)
